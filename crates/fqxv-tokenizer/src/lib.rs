//! Positional read-name tokenizer for FASTQ headers.
//!
//! Each name is split into maximal runs of digits / non-digits, and every token
//! is compared to the previous record's token at the same position:
//!
//! - identical bytes -> `MATCH` (constant instrument/run/lane/tile prefixes),
//! - numeric vs. a previous numeric -> `DELTA` (incrementing x/y coordinates),
//! - otherwise a literal string or number.
//!
//! Tokens are split into separate role streams (per-record token counts,
//! per-column ops, string lengths, string bytes, numeric literals, numeric
//! paddings, numeric deltas), each entropy-coded with [`fqxv_rans`] so every
//! stream models a clean distribution. The ops and numeric deltas are coded
//! **one rANS stream per token column**: the op and delta at a given column are
//! near-constant across records (column 3 is always DELTA, the x-coordinate
//! delta is always small), so each column collapses to almost nothing — far
//! better than one mixed stream, where the entropy coder can't see the
//! record-periodic structure. Every stream is compressed at both rANS orders
//! and the smaller kept (the order is self-describing in the stream header, so
//! the decoder needs no side channel). Round-trips are byte-exact.
//!
//! ```
//! use fqxv_tokenizer::{encode, decode};
//! let names: Vec<&[u8]> = vec![
//!     b"INST:1:FC:1:1101:1000:2000",
//!     b"INST:1:FC:1:1101:1005:2050",
//! ];
//! let enc = encode(&names).unwrap();
//! let out = decode(&enc).unwrap();
//! assert_eq!(out, names);
//! ```

use std::borrow::Cow;

use fqxv_rans::Order;
use thiserror::Error;

/// Errors returned by the tokenizer codec.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// The compressed stream was malformed or truncated.
    #[error("malformed name stream: {0}")]
    Malformed(&'static str),
    /// Underlying rANS coder failure.
    #[error(transparent)]
    Rans(#[from] fqxv_rans::Error),
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;

const FORMAT_VERSION: u8 = 1;
// Op codes. Ops are bucketed per token column, and a per-record token count
// delimits records, so no in-band end-of-record marker is needed.
const MATCH: u8 = 0;
const STR: u8 = 1;
const NUM: u8 = 2;
const DELTA: u8 = 3;
// Numeric runs longer than this don't fit i64; encode them as string literals.
const MAX_NUM_DIGITS: usize = 18;
// Ops and deltas are coded per token column (position within the name) up to
// this cap, so each column (tile / x / y …) is modeled on its own distribution.
const MAX_COL: usize = 63;

/// One name token. On the encode path `bytes` borrows a slice of the input name
/// (no per-token allocation for millions of names); on decode it owns rebuilt
/// bytes. `Cow` lets one type serve both without copying on encode.
#[derive(Clone)]
struct Tok<'a> {
    is_num: bool,
    bytes: Cow<'a, [u8]>,
    value: i64,
}

/// Compress `stream` at both rANS orders and return the smaller encoding. The
/// order is recorded in the rANS stream header, so [`fqxv_rans::decode`]
/// recovers it without any side channel — picking the best order per stream is
/// therefore free on the decode side.
fn encode_best(stream: &[u8]) -> Result<Vec<u8>> {
    let o0 = fqxv_rans::encode(stream, Order::Zero)?;
    let o1 = fqxv_rans::encode(stream, Order::One)?;
    Ok(if o1.len() < o0.len() { o1 } else { o0 })
}

/// Append a flat stream as `[varint comp_len, comp_bytes]`, best order.
fn push_stream(out: &mut Vec<u8>, stream: &[u8]) -> Result<()> {
    let c = encode_best(stream)?;
    write_varint(out, c.len() as u64);
    out.extend_from_slice(&c);
    Ok(())
}

/// Read one flat `[varint comp_len, comp_bytes]` stream and rANS-decode it.
fn read_stream(r: &mut Cursor<'_>) -> Result<Vec<u8>> {
    let len = r.varint()? as usize;
    Ok(fqxv_rans::decode(r.take(len)?)?)
}

/// Append the per-column buckets, each as its own `[varint comp_len, bytes]`
/// stream (empty columns collapse to a single `0` byte). Coding each column
/// separately lets its near-constant run compress to almost nothing.
fn push_cols(out: &mut Vec<u8>, cols: &[Vec<u8>]) -> Result<()> {
    for col in cols {
        if col.is_empty() {
            write_varint(out, 0);
        } else {
            let c = encode_best(col)?;
            write_varint(out, c.len() as u64);
            out.extend_from_slice(&c);
        }
    }
    Ok(())
}

/// Inverse of [`push_cols`]: decode `MAX_COL + 1` per-column buckets.
fn read_cols(r: &mut Cursor<'_>) -> Result<Vec<Vec<u8>>> {
    let mut cols = Vec::with_capacity(MAX_COL + 1);
    for _ in 0..=MAX_COL {
        let len = r.varint()? as usize;
        if len == 0 {
            cols.push(Vec::new());
        } else {
            cols.push(fqxv_rans::decode(r.take(len)?)?);
        }
    }
    Ok(cols)
}

/// Number of decimal digits in a non-negative numeric token's canonical form,
/// used to derive a token's width from its value (only genuine leading-zero
/// padding then needs to be stored).
fn natural_digits(value: i64) -> usize {
    value.to_string().len()
}

/// Encode a list of read names.
///
/// Tokens are split into separate role streams — per-column ops, string lengths,
/// string bytes, numeric literals, leading-zero paddings, and per-column numeric
/// deltas — each entropy coded on its own so rANS models a clean distribution
/// per stream (in particular the incrementing x/y-coordinate deltas compress far
/// better apart from string bytes, and per column rather than mixed together).
pub fn encode(names: &[&[u8]]) -> Result<Vec<u8>> {
    let mut counts: Vec<u8> = Vec::new();
    let mut op_cols: Vec<Vec<u8>> = vec![Vec::new(); MAX_COL + 1];
    let mut str_lens: Vec<u8> = Vec::new();
    let mut str_data: Vec<u8> = Vec::new();
    let mut num_vals: Vec<u8> = Vec::new();
    // Leading-zero padding per numeric token: width minus the token's natural
    // digit count. Almost always 0 (names rarely zero-pad), so this stream is
    // near-degenerate — far cheaper than storing the full width.
    let mut pads: Vec<u8> = Vec::new();
    let mut delta_cols: Vec<Vec<u8>> = vec![Vec::new(); MAX_COL + 1];
    let mut prev: Vec<Tok> = Vec::new();

    for name in names {
        let toks = tokenize(name);
        for (i, t) in toks.iter().enumerate() {
            // Ops go in a per-column bucket: the op at a given column is almost
            // always the same across records (col 3 is always DELTA, …), so each
            // column collapses to near-nothing once coded on its own.
            let col = i.min(MAX_COL);
            let p = prev.get(i);
            if p.is_some_and(|p| p.bytes == t.bytes) {
                op_cols[col].push(MATCH);
            } else if t.is_num {
                let pad = t.bytes.len() - natural_digits(t.value);
                if let Some(p) = p.filter(|p| p.is_num) {
                    op_cols[col].push(DELTA);
                    write_varint(&mut delta_cols[col], zigzag(t.value - p.value));
                } else {
                    op_cols[col].push(NUM);
                    write_varint(&mut num_vals, t.value as u64);
                }
                pads.push(pad as u8);
            } else {
                op_cols[col].push(STR);
                write_varint(&mut str_lens, t.bytes.len() as u64);
                str_data.extend_from_slice(&t.bytes);
            }
        }
        write_varint(&mut counts, toks.len() as u64);
        prev = toks;
    }

    let mut out = Vec::new();
    out.push(FORMAT_VERSION);
    write_varint(&mut out, names.len() as u64);
    // Flat streams first; order must match `decode`.
    push_stream(&mut out, &counts)?;
    push_stream(&mut out, &str_lens)?;
    push_stream(&mut out, &str_data)?;
    push_stream(&mut out, &num_vals)?;
    push_stream(&mut out, &pads)?;
    // Then the per-column op and delta buckets.
    push_cols(&mut out, &op_cols)?;
    push_cols(&mut out, &delta_cols)?;
    Ok(out)
}

/// Decode a stream produced by [`encode`], returning the read names.
pub fn decode(src: &[u8]) -> Result<Vec<Vec<u8>>> {
    let mut r = Cursor::new(src);
    if r.u8()? != FORMAT_VERSION {
        return Err(Error::Malformed("unsupported version"));
    }
    let n_records = r.varint()? as usize;

    let counts = read_stream(&mut r)?;
    let str_lens = read_stream(&mut r)?;
    let str_data = read_stream(&mut r)?;
    let num_vals = read_stream(&mut r)?;
    let pads = read_stream(&mut r)?;
    let op_data = read_cols(&mut r)?;
    let delta_data = read_cols(&mut r)?;

    let mut op_cursors: Vec<Cursor> = op_data.iter().map(|v| Cursor::new(v)).collect();
    let mut delta_cursors: Vec<Cursor> = delta_data.iter().map(|v| Cursor::new(v)).collect();
    let mut c_counts = Cursor::new(&counts);
    let mut c_str_lens = Cursor::new(&str_lens);
    let mut c_num = Cursor::new(&num_vals);
    let (mut str_pos, mut pad_pos) = (0usize, 0usize);
    let next_pad = |pad_pos: &mut usize| -> Result<usize> {
        let p = *pads.get(*pad_pos).ok_or(Error::Malformed("pad underrun"))?;
        *pad_pos += 1;
        Ok(p as usize)
    };

    let mut names = Vec::with_capacity(n_records.min(1 << 20));
    let mut prev: Vec<Tok<'static>> = Vec::new();
    let mut cur: Vec<Tok<'static>> = Vec::new();

    for _ in 0..n_records {
        let n_toks = c_counts.varint()? as usize;
        for _ in 0..n_toks {
            let col = cur.len().min(MAX_COL);
            let op = op_cursors[col].u8()?;
            match op {
                MATCH => {
                    let t = prev
                        .get(cur.len())
                        .ok_or(Error::Malformed("MATCH without prior token"))?
                        .clone();
                    cur.push(t);
                }
                STR => {
                    let len = c_str_lens.varint()? as usize;
                    let bytes = str_data
                        .get(str_pos..str_pos + len)
                        .ok_or(Error::Malformed("string data underrun"))?
                        .to_vec();
                    str_pos += len;
                    cur.push(Tok {
                        is_num: false,
                        bytes: Cow::Owned(bytes),
                        value: 0,
                    });
                }
                NUM => {
                    let value = c_num.varint()? as i64;
                    let width = natural_digits(value) + next_pad(&mut pad_pos)?;
                    cur.push(num_tok(value, width));
                }
                DELTA => {
                    let d = unzigzag(delta_cursors[col].varint()?);
                    let p = prev
                        .get(cur.len())
                        .filter(|p| p.is_num)
                        .ok_or(Error::Malformed("DELTA without numeric prior"))?;
                    let value = p.value + d;
                    let width = natural_digits(value) + next_pad(&mut pad_pos)?;
                    cur.push(num_tok(value, width));
                }
                _ => return Err(Error::Malformed("unknown op")),
            }
        }
        let name = cur.iter().flat_map(|t| t.bytes.iter().copied()).collect();
        names.push(name);
        prev = std::mem::take(&mut cur);
    }
    Ok(names)
}

/// Split a name into maximal digit / non-digit runs. Tokens borrow slices of
/// `name`, so the common encode path allocates nothing per token.
fn tokenize(name: &[u8]) -> Vec<Tok<'_>> {
    let mut toks = Vec::new();
    let mut i = 0;
    while i < name.len() {
        let is_digit = name[i].is_ascii_digit();
        let start = i;
        while i < name.len() && name[i].is_ascii_digit() == is_digit {
            i += 1;
        }
        let bytes = &name[start..i];
        // Only treat as numeric (for delta) if it fits i64.
        let is_num = is_digit && bytes.len() <= MAX_NUM_DIGITS;
        let value = if is_num {
            std::str::from_utf8(bytes)
                .ok()
                .and_then(|s| s.parse::<i64>().ok())
                .unwrap_or(0)
        } else {
            0
        };
        toks.push(Tok {
            is_num,
            bytes: Cow::Borrowed(bytes),
            value,
        });
    }
    toks
}

/// Reconstruct a numeric token's bytes, zero-padded to its original width.
fn num_tok(value: i64, width: usize) -> Tok<'static> {
    let s = value.to_string();
    let bytes = if s.len() >= width {
        s.into_bytes()
    } else {
        let mut v = vec![b'0'; width - s.len()];
        v.extend_from_slice(s.as_bytes());
        v
    };
    Tok {
        is_num: true,
        bytes: Cow::Owned(bytes),
        value,
    }
}

fn zigzag(d: i64) -> u64 {
    ((d << 1) ^ (d >> 63)) as u64
}

fn unzigzag(z: u64) -> i64 {
    ((z >> 1) as i64) ^ -((z & 1) as i64)
}

fn write_varint(out: &mut Vec<u8>, mut v: u64) {
    loop {
        let byte = (v & 0x7f) as u8;
        v >>= 7;
        if v == 0 {
            out.push(byte);
            break;
        }
        out.push(byte | 0x80);
    }
}

struct Cursor<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn new(buf: &'a [u8]) -> Self {
        Cursor { buf, pos: 0 }
    }
    fn u8(&mut self) -> Result<u8> {
        let b = *self
            .buf
            .get(self.pos)
            .ok_or(Error::Malformed("truncated"))?;
        self.pos += 1;
        Ok(b)
    }
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self.pos + n;
        let s = self
            .buf
            .get(self.pos..end)
            .ok_or(Error::Malformed("truncated slice"))?;
        self.pos = end;
        Ok(s)
    }
    fn varint(&mut self) -> Result<u64> {
        let mut v = 0u64;
        let mut shift = 0u32;
        loop {
            let byte = self.u8()?;
            v |= u64::from(byte & 0x7f) << shift;
            if byte & 0x80 == 0 {
                return Ok(v);
            }
            shift += 7;
            if shift >= 64 {
                return Err(Error::Malformed("varint too long"));
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(names: &[&[u8]]) {
        let enc = encode(names).expect("encode");
        let out = decode(&enc).expect("decode");
        let expect: Vec<Vec<u8>> = names.iter().map(|n| n.to_vec()).collect();
        assert_eq!(out, expect, "name round-trip mismatch");
    }

    #[test]
    fn roundtrip_empty() {
        roundtrip(&[]);
    }

    #[test]
    fn roundtrip_illumina_incrementing() {
        roundtrip(&[
            b"INST:1:FC:1:1101:1000:2000",
            b"INST:1:FC:1:1101:1005:2050",
            b"INST:1:FC:1:1101:1010:2100",
        ]);
    }

    #[test]
    fn roundtrip_leading_zeros_and_punct() {
        roundtrip(&[
            b"SRR453566.1 HWI-ST167:4:1101:0042:1986 length=101",
            b"SRR453566.2 HWI-ST167:4:1101:0043:1990 length=101",
        ]);
    }

    #[test]
    fn roundtrip_leading_zero_widths() {
        // Zero-padded numerics of differing widths must survive exactly, since
        // width is derived from value + stored padding.
        roundtrip(&[
            b"X:007:0000:5",
            b"X:008:0001:5",
            b"X:009:0002:5",
            b"X:010:1000:5",
        ]);
    }

    #[test]
    fn roundtrip_varying_structure() {
        roundtrip(&[b"a", b"bb99", b"", b"12", b"x1y2z3"]);
    }

    #[test]
    fn incrementing_names_compress_well() {
        // Enough names that the per-stream rANS table overhead amortizes.
        let names: Vec<Vec<u8>> = (0..20000)
            .map(|i| format!("INST:1:FC:1:1101:{}:{}", 1000 + i, 2000 + i * 2).into_bytes())
            .collect();
        let refs: Vec<&[u8]> = names.iter().map(|n| n.as_slice()).collect();
        let enc = encode(&refs).expect("encode");
        let raw: usize = names.iter().map(|n| n.len() + 1).sum();
        assert!(
            enc.len() < raw / 10,
            "expected >10x on incrementing names, got {raw} -> {}",
            enc.len()
        );
        assert_eq!(decode(&enc).unwrap(), names);
    }

    proptest::proptest! {
        #[test]
        fn roundtrip_arbitrary(
            names in proptest::collection::vec(
                proptest::collection::vec(
                    proptest::sample::select(b"AB:._-0129 ".to_vec()), 0..40),
                0..60)
        ) {
            let refs: Vec<&[u8]> = names.iter().map(|n| n.as_slice()).collect();
            let enc = encode(&refs).expect("encode");
            let out = decode(&enc).expect("decode");
            proptest::prop_assert_eq!(out, names);
        }
    }
}
