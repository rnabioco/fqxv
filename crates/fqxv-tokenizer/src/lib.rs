//! Positional read-name tokenizer for FASTQ headers.
//!
//! Each name is split into maximal runs of digits / non-digits, and every token
//! is compared to the previous record's token at the same position:
//!
//! - identical bytes -> `MATCH` (constant instrument/run/lane/tile prefixes),
//! - numeric vs. a previous numeric -> `DELTA` (incrementing x/y coordinates),
//! - otherwise a literal string or number.
//!
//! Tokens are split into separate role streams (ops, string lengths, string
//! bytes, numeric literals, token widths, numeric deltas), each entropy-coded
//! with [`fqxv_rans`] so every stream models a clean distribution — the
//! incrementing coordinate deltas in particular compress far better on their
//! own than mixed with string bytes. Round-trips are byte-exact.
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

const FORMAT_VERSION: u8 = 0;
// Op codes.
const MATCH: u8 = 0;
const STR: u8 = 1;
const NUM: u8 = 2;
const DELTA: u8 = 3;
const REC_END: u8 = 4;
// Numeric runs longer than this don't fit i64; encode them as string literals.
const MAX_NUM_DIGITS: usize = 18;

#[derive(Clone)]
struct Tok {
    is_num: bool,
    bytes: Vec<u8>,
    value: i64,
}

/// Encode a list of read names.
///
/// Tokens are split into separate role streams — ops, string lengths, string
/// bytes, numeric literals, token widths, and numeric deltas — each entropy
/// coded on its own so rANS models a clean distribution per stream (in
/// particular the incrementing x/y-coordinate deltas compress far better apart
/// from string bytes).
pub fn encode(names: &[&[u8]]) -> Result<Vec<u8>> {
    let mut ops: Vec<u8> = Vec::new();
    let mut str_lens: Vec<u8> = Vec::new();
    let mut str_data: Vec<u8> = Vec::new();
    let mut num_vals: Vec<u8> = Vec::new();
    let mut widths: Vec<u8> = Vec::new();
    let mut delta_vals: Vec<u8> = Vec::new();
    let mut prev: Vec<Tok> = Vec::new();

    for name in names {
        let toks = tokenize(name);
        for (i, t) in toks.iter().enumerate() {
            let p = prev.get(i);
            if p.is_some_and(|p| p.bytes == t.bytes) {
                ops.push(MATCH);
            } else if t.is_num {
                if let Some(p) = p.filter(|p| p.is_num) {
                    ops.push(DELTA);
                    write_varint(&mut delta_vals, zigzag(t.value - p.value));
                    widths.push(t.bytes.len() as u8);
                } else {
                    ops.push(NUM);
                    write_varint(&mut num_vals, t.value as u64);
                    widths.push(t.bytes.len() as u8);
                }
            } else {
                ops.push(STR);
                write_varint(&mut str_lens, t.bytes.len() as u64);
                str_data.extend_from_slice(&t.bytes);
            }
        }
        ops.push(REC_END);
        prev = toks;
    }

    let mut out = Vec::new();
    out.push(FORMAT_VERSION);
    write_varint(&mut out, names.len() as u64);
    // Stream order must match `decode`; text-like streams use order-1.
    for (stream, order) in [
        (&ops, Order::One),
        (&str_lens, Order::Zero),
        (&str_data, Order::One),
        (&num_vals, Order::Zero),
        (&widths, Order::Zero),
        (&delta_vals, Order::Zero),
    ] {
        let c = fqxv_rans::encode(stream, order)?;
        write_varint(&mut out, c.len() as u64);
        out.extend_from_slice(&c);
    }
    Ok(out)
}

/// Decode a stream produced by [`encode`], returning the read names.
pub fn decode(src: &[u8]) -> Result<Vec<Vec<u8>>> {
    let mut r = Cursor::new(src);
    if r.u8()? != FORMAT_VERSION {
        return Err(Error::Malformed("unsupported version"));
    }
    let n_records = r.varint()? as usize;
    let mut stream = || -> Result<Vec<u8>> {
        let len = r.varint()? as usize;
        Ok(fqxv_rans::decode(r.take(len)?)?)
    };
    let ops = stream()?;
    let str_lens = stream()?;
    let str_data = stream()?;
    let num_vals = stream()?;
    let widths = stream()?;
    let delta_vals = stream()?;

    let mut c_str_lens = Cursor::new(&str_lens);
    let mut c_num = Cursor::new(&num_vals);
    let mut c_delta = Cursor::new(&delta_vals);
    let mut str_pos = 0usize;
    let mut w_pos = 0usize;
    let next_width = |w_pos: &mut usize| -> Result<usize> {
        let w = *widths
            .get(*w_pos)
            .ok_or(Error::Malformed("width underrun"))?;
        *w_pos += 1;
        Ok(w as usize)
    };

    let mut names = Vec::with_capacity(n_records.min(1 << 20));
    let mut prev: Vec<Tok> = Vec::new();
    let mut cur: Vec<Tok> = Vec::new();

    for &op in &ops {
        match op {
            REC_END => {
                let name = cur.iter().flat_map(|t| t.bytes.iter().copied()).collect();
                names.push(name);
                prev = std::mem::take(&mut cur);
            }
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
                    bytes,
                    value: 0,
                });
            }
            NUM => {
                let value = c_num.varint()? as i64;
                let width = next_width(&mut w_pos)?;
                cur.push(num_tok(value, width));
            }
            DELTA => {
                let d = unzigzag(c_delta.varint()?);
                let width = next_width(&mut w_pos)?;
                let p = prev
                    .get(cur.len())
                    .filter(|p| p.is_num)
                    .ok_or(Error::Malformed("DELTA without numeric prior"))?;
                cur.push(num_tok(p.value + d, width));
            }
            _ => return Err(Error::Malformed("unknown op")),
        }
    }
    Ok(names)
}

/// Split a name into maximal digit / non-digit runs.
fn tokenize(name: &[u8]) -> Vec<Tok> {
    let mut toks = Vec::new();
    let mut i = 0;
    while i < name.len() {
        let is_digit = name[i].is_ascii_digit();
        let start = i;
        while i < name.len() && name[i].is_ascii_digit() == is_digit {
            i += 1;
        }
        let bytes = name[start..i].to_vec();
        // Only treat as numeric (for delta) if it fits i64.
        let is_num = is_digit && bytes.len() <= MAX_NUM_DIGITS;
        let value = if is_num {
            std::str::from_utf8(&bytes)
                .ok()
                .and_then(|s| s.parse::<i64>().ok())
                .unwrap_or(0)
        } else {
            0
        };
        toks.push(Tok {
            is_num,
            bytes,
            value,
        });
    }
    toks
}

/// Reconstruct a numeric token's bytes, zero-padded to its original width.
fn num_tok(value: i64, width: usize) -> Tok {
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
        bytes,
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
