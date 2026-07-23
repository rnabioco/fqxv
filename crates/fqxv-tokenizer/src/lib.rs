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

use fqxv_bytes::{ReaderError, unzigzag, write_varint, zigzag};
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

// Cap on any single rANS sub-stream's decoded length, as a multiple of the whole
// compressed names stream (with a floor for tiny inputs). Names decode one block
// at a time, and no individual stream legitimately expands past this — the large
// expansion in repetitive data comes from token reconstruction (MATCH/DELTA
// replays), not from one stream. On the untrusted direct-decode path this stops
// a corrupt output-length header from allocating gigabytes; see
// [`fqxv_rans::decode_bounded`].
const MAX_STREAM_EXPANSION: usize = 1 << 16;
const MIN_STREAM_BOUND: usize = 1 << 20;

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

/// Read one flat `[varint comp_len, comp_bytes]` stream and rANS-decode it,
/// rejecting a decoded length above `max_out` before it is allocated.
fn read_stream(r: &mut Cursor<'_>, max_out: usize) -> Result<Vec<u8>> {
    let len = r.varint()? as usize;
    Ok(fqxv_rans::decode_bounded(r.take(len)?, max_out)?)
}

/// Append the per-column numeric deltas as fixed-width little-endian byte
/// planes: for each column, `[varint width]` then one best-order stream per
/// byte plane. Splitting a value's low/high bytes into separate streams lets
/// each plane model its own distribution — the high planes of x/y coordinate
/// deltas are near-zero and compress away, which coding the packed varints in
/// one model cannot reach. `width == 0` marks an empty column.
fn push_delta_planes(out: &mut Vec<u8>, cols: &[Vec<u64>]) -> Result<()> {
    for col in cols {
        if col.is_empty() {
            write_varint(out, 0);
            continue;
        }
        let maxv = col.iter().copied().max().unwrap();
        let width = (64 - maxv.leading_zeros()).div_ceil(8).max(1) as usize;
        write_varint(out, width as u64);
        for plane in 0..width {
            let s: Vec<u8> = col.iter().map(|&v| (v >> (8 * plane)) as u8).collect();
            let c = encode_best(&s)?;
            write_varint(out, c.len() as u64);
            out.extend_from_slice(&c);
        }
    }
    Ok(())
}

/// Inverse of [`push_delta_planes`]. `n_delta[col]` (derived from the op stream:
/// one delta per `DELTA` op in that column) gives each column's value count, so
/// no per-column length is stored; it also bounds the per-column allocation.
fn read_delta_planes(r: &mut Cursor<'_>, n_delta: &[usize]) -> Result<Vec<Vec<u64>>> {
    let mut cols = Vec::with_capacity(MAX_COL + 1);
    for &n in n_delta.iter().take(MAX_COL + 1) {
        let width = r.varint()? as usize;
        if width == 0 {
            cols.push(Vec::new());
            continue;
        }
        // `n` derives from the untrusted `n_records` header, so reserve fallibly:
        // a corrupt count must error, not abort on a huge infallible allocation.
        let mut vals: Vec<u64> = Vec::new();
        vals.try_reserve_exact(n)
            .map_err(|_| Error::Malformed("delta plane too large to allocate"))?;
        vals.resize(n, 0);
        for plane in 0..width {
            let len = r.varint()? as usize;
            // Each plane must decode to exactly `n` bytes, so bound the decode by
            // `n` up front — a corrupt length header can't over-allocate.
            let bytes = fqxv_rans::decode_bounded(r.take(len)?, n)?;
            if bytes.len() != n {
                return Err(Error::Malformed("delta plane length mismatch"));
            }
            for (v, &b) in vals.iter_mut().zip(bytes.iter()) {
                *v |= u64::from(b) << (8 * plane);
            }
        }
        cols.push(vals);
    }
    Ok(cols)
}

/// RLE-encode the op columns into (per-column run counts, run symbols, run
/// lengths), each its own best-order stream. Op columns are near-periodic —
/// paired mates share x/y coordinates, so a coordinate column alternates
/// MATCH/DELTA — so a handful of run symbols and small run lengths carry the
/// whole thing, an order of magnitude below coding the raw ops.
fn push_op_rle(out: &mut Vec<u8>, cols: &[Vec<u8>]) -> Result<()> {
    let mut counts = Vec::new();
    let mut syms = Vec::new();
    let mut lens = Vec::new();
    for col in cols {
        let mut runs = 0u64;
        let mut i = 0;
        while i < col.len() {
            let s = col[i];
            let mut run = 1u64;
            while i + (run as usize) < col.len() && col[i + run as usize] == s {
                run += 1;
            }
            syms.push(s);
            write_varint(&mut lens, run);
            runs += 1;
            i += run as usize;
        }
        write_varint(&mut counts, runs);
    }
    push_stream(out, &counts)?;
    push_stream(out, &syms)?;
    push_stream(out, &lens)?;
    Ok(())
}

/// Inverse of [`push_op_rle`]. `n_records` bounds each column's length (a column
/// holds at most one op per record), guarding a corrupt run length from
/// allocating unboundedly.
fn read_op_rle(r: &mut Cursor<'_>, n_records: usize, max_out: usize) -> Result<Vec<Vec<u8>>> {
    let counts = read_stream(r, max_out)?;
    let syms = read_stream(r, max_out)?;
    let lens = read_stream(r, max_out)?;
    let mut c_counts = Cursor::new(&counts);
    let mut c_lens = Cursor::new(&lens);
    let mut sym_pos = 0usize;
    let mut cols = Vec::with_capacity(MAX_COL + 1);
    // Aggregate op bytes across every column. The per-column `n_records` guard
    // below bounds each of the 65 columns alone, but 65 columns can still reach
    // `65 * n_records` — and `n_records` itself can be as large as `max_out` — so
    // the columns together (and the delta planes they later size) could allocate
    // ~65x past the per-stream cap every other stream obeys (#142). The op
    // content is one logical stream; hold its total to `max_out` too.
    let mut total_ops = 0usize;
    for _ in 0..=MAX_COL {
        let n_runs = c_counts.varint()?;
        let mut col = Vec::new();
        for _ in 0..n_runs {
            let s = *syms
                .get(sym_pos)
                .ok_or(Error::Malformed("op sym underrun"))?;
            sym_pos += 1;
            let run = c_lens.varint()? as usize;
            // `n_records` is bounded by the token stream (see `decode`), so this
            // is a tight cap on the column length. `checked_add` also stops a
            // near-`usize::MAX` run from overflowing the sum past the check.
            if col
                .len()
                .checked_add(run)
                .is_none_or(|total| total > n_records)
            {
                return Err(Error::Malformed("op run too long"));
            }
            total_ops = total_ops
                .checked_add(run)
                .filter(|&t| t <= max_out)
                .ok_or(Error::Malformed("op columns exceed the stream bound"))?;
            // Grow fallibly so even a within-bound run can't abort on OOM.
            col.try_reserve(run)
                .map_err(|_| Error::Malformed("op run too large to allocate"))?;
            col.resize(col.len() + run, s);
        }
        cols.push(col);
    }
    Ok(cols)
}

/// Number of decimal digits in a non-negative numeric token's canonical form,
/// used to derive a token's width from its value (only genuine leading-zero
/// padding then needs to be stored).
fn natural_digits(value: i64) -> usize {
    // Arithmetic digit count avoids a per-token `String` allocation in the
    // hot encode/decode paths. Tokens are parsed from digit runs, so `value`
    // is non-negative; `ilog10` is undefined at 0, hence the special case.
    if value <= 0 {
        1
    } else {
        value.ilog10() as usize + 1
    }
}

/// One column of a regenerable read-name template.
#[derive(Debug, Clone, PartialEq, Eq)]
enum TemplateColumn {
    /// Constant bytes (separators, prefixes, and constant numeric fields).
    Const(Vec<u8>),
    /// A per-read counter: the value at position `j` is `start + j`, rendered in
    /// decimal and left-zero-padded to `pad` (0 = the natural, unpadded width).
    Counter { start: i64, pad: usize },
}

/// A read-name template that regenerates each name purely from its position — a
/// name is the concatenation of its columns with counters evaluated at the
/// position. It is detected (see [`detect_template`]) only when every name shares
/// one digit/non-digit structure and each numeric column is either constant or a
/// `start + row_index` counter, i.e. the names carry nothing but read order (as
/// with SRA `@RUN.N N …` headers). Regenerating them makes the names stream
/// essentially free — at the cost of renumbering reads, so it is reorder-lossy
/// and must be opt-in.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NameTemplate {
    columns: Vec<TemplateColumn>,
}

/// Upper bound on template columns (bounds header size; names with more
/// digit/non-digit runs are treated as not regenerable).
const MAX_TEMPLATE_COLS: usize = 64;

/// Upper bound on a counter column's zero-pad width. A real read-name numeric
/// field is a handful of digits; this cap (comfortably above any genuine value)
/// stops a corrupt template from driving an unbounded fill in [`NameTemplate::regenerate`].
const MAX_COUNTER_PAD: usize = 1024;

impl NameTemplate {
    /// A synthetic renumbering template: the name at output position `j` is the
    /// bare decimal counter `j + 1`. Used by the reorder-lossy renumber path when
    /// [`detect_template`] finds no regenerable structure in the original names, so
    /// the name stream AND the order permutation are both dropped (SPRING-style
    /// renumbering) rather than stored losslessly. The reads themselves (sequence
    /// and quality) are preserved exactly; only the original names are discarded.
    #[must_use]
    pub fn renumber() -> Self {
        NameTemplate {
            columns: vec![TemplateColumn::Counter { start: 1, pad: 0 }],
        }
    }

    /// Regenerate the name at position `index` (counters evaluated at `index`).
    #[must_use]
    pub fn regenerate(&self, index: usize) -> Vec<u8> {
        let mut out = Vec::new();
        for col in &self.columns {
            match col {
                TemplateColumn::Const(b) => out.extend_from_slice(b),
                TemplateColumn::Counter { start, pad } => {
                    let s = (start + index as i64).to_string();
                    if s.len() < *pad {
                        out.extend(std::iter::repeat_n(b'0', pad - s.len()));
                    }
                    out.extend_from_slice(s.as_bytes());
                }
            }
        }
        out
    }

    /// Serialize the template for storage in a container header.
    #[must_use]
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut out = Vec::new();
        write_varint(&mut out, self.columns.len() as u64);
        for col in &self.columns {
            match col {
                TemplateColumn::Const(b) => {
                    out.push(0);
                    write_varint(&mut out, b.len() as u64);
                    out.extend_from_slice(b);
                }
                TemplateColumn::Counter { start, pad } => {
                    out.push(1);
                    write_varint(&mut out, zigzag(*start));
                    write_varint(&mut out, *pad as u64);
                }
            }
        }
        out
    }

    /// Deserialize a template written by [`NameTemplate::to_bytes`].
    pub fn from_bytes(src: &[u8]) -> Result<Self> {
        let mut r = Cursor::new(src);
        let ncol = r.varint()? as usize;
        if ncol > MAX_TEMPLATE_COLS {
            return Err(Error::Malformed("name template too large"));
        }
        let mut columns = Vec::with_capacity(ncol);
        for _ in 0..ncol {
            match r.u8()? {
                0 => {
                    let len = r.varint()? as usize;
                    columns.push(TemplateColumn::Const(r.take(len)?.to_vec()));
                }
                1 => {
                    let start = unzigzag(r.varint()?);
                    let pad = r.varint()? as usize;
                    if pad > MAX_COUNTER_PAD {
                        return Err(Error::Malformed("name-template counter pad too large"));
                    }
                    columns.push(TemplateColumn::Counter { start, pad });
                }
                _ => return Err(Error::Malformed("bad name-template column tag")),
            }
        }
        Ok(NameTemplate { columns })
    }
}

/// Detect a regenerable counter template over `names` (given in ORIGINAL order),
/// or `None` if the names carry more than their position. Every name must share
/// the same digit/non-digit column structure; each column must be either
/// identical across all names, or a numeric `first + row_index` counter (natural
/// or fixed-zero-padded width). This is the common SRA case (`@RUN.N N length=L`)
/// where the number is just the read index; other schemes (Illumina tile/x/y)
/// return `None`.
#[must_use]
pub fn detect_template(names: &[&[u8]]) -> Option<NameTemplate> {
    let first = tokenize(names.first()?);
    let ncol = first.len();
    if ncol == 0 || ncol > MAX_TEMPLATE_COLS {
        return None;
    }
    struct Col {
        first_bytes: Vec<u8>,
        first_val: i64,
        is_num: bool,
        const_ok: bool,
        nat_ok: bool,   // natural-width counter still possible
        fixed_ok: bool, // fixed-width (leading-zero) counter still possible
        width: usize,
    }
    let mut cols: Vec<Col> = first
        .iter()
        .map(|t| Col {
            first_bytes: t.bytes.to_vec(),
            first_val: t.value,
            is_num: t.is_num,
            const_ok: true,
            nat_ok: t.is_num,
            fixed_ok: t.is_num,
            width: t.bytes.len(),
        })
        .collect();

    for (row, name) in names.iter().enumerate() {
        let toks = tokenize(name);
        if toks.len() != ncol {
            return None;
        }
        for (c, t) in cols.iter_mut().zip(&toks) {
            if t.is_num != c.is_num {
                return None; // structure must be stable across all names
            }
            if c.const_ok && t.bytes.as_ref() != c.first_bytes.as_slice() {
                c.const_ok = false;
            }
            if c.is_num && (c.nat_ok || c.fixed_ok) {
                if c.first_val.checked_add(row as i64) != Some(t.value) {
                    c.nat_ok = false;
                    c.fixed_ok = false;
                    continue;
                }
                let s = t.value.to_string();
                if t.bytes.as_ref() != s.as_bytes() {
                    c.nat_ok = false;
                }
                if t.bytes.len() != c.width || s.len() > c.width {
                    c.fixed_ok = false;
                } else if c.fixed_ok {
                    let mut b = vec![b'0'; c.width - s.len()];
                    b.extend_from_slice(s.as_bytes());
                    if t.bytes.as_ref() != b.as_slice() {
                        c.fixed_ok = false;
                    }
                }
            }
        }
    }

    let mut columns = Vec::with_capacity(ncol);
    for c in cols {
        if c.const_ok {
            columns.push(TemplateColumn::Const(c.first_bytes));
        } else if c.nat_ok {
            columns.push(TemplateColumn::Counter {
                start: c.first_val,
                pad: 0,
            });
        } else if c.fixed_ok {
            columns.push(TemplateColumn::Counter {
                start: c.first_val,
                pad: c.width,
            });
        } else {
            return None; // this column carries more than position
        }
    }
    Some(NameTemplate { columns })
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
    // Zigzag numeric deltas per column, coded later as byte planes.
    let mut delta_vals: Vec<Vec<u64>> = vec![Vec::new(); MAX_COL + 1];
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
                    delta_vals[col].push(zigzag(t.value - p.value));
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
    // Ops are RLE'd (near-periodic runs); deltas keep their raw per-column
    // coding (they carry the genuine coordinate entropy, not runs).
    push_op_rle(&mut out, &op_cols)?;
    push_delta_planes(&mut out, &delta_vals)?;
    Ok(out)
}

/// Decode a stream produced by [`encode`], returning the read names.
pub fn decode(src: &[u8]) -> Result<Vec<Vec<u8>>> {
    let mut r = Cursor::new(src);
    if r.u8()? != FORMAT_VERSION {
        return Err(Error::Malformed("unsupported version"));
    }
    let n_records = r.varint()? as usize;

    // Per-stream decode cap for the untrusted direct-decode path (see
    // [`MAX_STREAM_EXPANSION`]). Well above any legitimate single stream, but
    // finite, so a corrupt length header cannot allocate gigabytes.
    let max_out = src
        .len()
        .saturating_mul(MAX_STREAM_EXPANSION)
        .max(MIN_STREAM_BOUND);

    let counts = read_stream(&mut r, max_out)?;
    // `counts` holds exactly one token-count varint per record (>= 1 byte each),
    // so a valid `n_records` can never exceed `counts.len()`. Enforcing that here
    // turns the `col.len() + run > n_records` check in `read_op_rle` and the
    // delta-plane sizing into tight bounds: without it a corrupt `n_records` lets
    // a crafted run drive a multi-gigabyte zero-fill that still fits in memory
    // (so `try_reserve` never trips) and stalls decode for minutes.
    if n_records > counts.len() {
        return Err(Error::Malformed("record count exceeds token stream"));
    }
    let str_lens = read_stream(&mut r, max_out)?;
    let str_data = read_stream(&mut r, max_out)?;
    let num_vals = read_stream(&mut r, max_out)?;
    let pads = read_stream(&mut r, max_out)?;
    let op_data = read_op_rle(&mut r, n_records, max_out)?;
    // Each column holds one delta per DELTA op in it; that count sizes the byte
    // planes, so it need not be stored.
    let n_delta: Vec<usize> = op_data
        .iter()
        .map(|col| col.iter().filter(|&&o| o == DELTA).count())
        .collect();
    let delta_vals = read_delta_planes(&mut r, &n_delta)?;

    let mut op_cursors: Vec<Cursor> = op_data.iter().map(|v| Cursor::new(v)).collect();
    let mut delta_pos = vec![0usize; MAX_COL + 1];
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
                    let z = *delta_vals[col]
                        .get(delta_pos[col])
                        .ok_or(Error::Malformed("delta underrun"))?;
                    delta_pos[col] += 1;
                    let d = unzigzag(z);
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

/// Shared byte cursor specialized to this crate's [`Error`].
type Cursor<'a> = fqxv_bytes::Reader<'a, Error>;

impl ReaderError for Error {
    fn truncated() -> Self {
        Error::Malformed("truncated")
    }
    fn bad_varint() -> Self {
        Error::Malformed("varint too long")
    }
    fn oversized() -> Self {
        Error::Malformed("length count too large to allocate")
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

    /// A corrupt stream must not be able to drive decode into a huge zero-fill.
    ///
    /// Hand-build a well-formed stream whose op-RLE claims a ~1 GiB run while the
    /// header inflates `n_records` so the `col.len() + run > n_records` check
    /// would pass. The record-count bound in `decode` must reject it up front
    /// rather than `resize`-ing gigabytes (which fits in memory, so `try_reserve`
    /// would not catch it).
    #[test]
    fn corrupt_run_length_cannot_drive_huge_allocation() {
        let mut out = Vec::new();
        out.push(FORMAT_VERSION);
        write_varint(&mut out, u64::MAX / 2); // inflated n_records

        // Flat streams the decoder reads before the op RLE. `counts` is short —
        // that is what now bounds `n_records`; the rest are empty.
        let mut counts = Vec::new();
        write_varint(&mut counts, 1); // record 0 has one token
        push_stream(&mut out, &counts).unwrap();
        push_stream(&mut out, &[]).unwrap(); // str_lens
        push_stream(&mut out, &[]).unwrap(); // str_data
        push_stream(&mut out, &[]).unwrap(); // num_vals
        push_stream(&mut out, &[]).unwrap(); // pads

        // Op RLE: column 0 has one run of ~1 GiB; the other columns are empty.
        let mut op_counts = Vec::new();
        for c in 0..=MAX_COL {
            write_varint(&mut op_counts, u64::from(c == 0));
        }
        let syms = vec![MATCH];
        let mut op_lens = Vec::new();
        write_varint(&mut op_lens, 1 << 30);
        push_stream(&mut out, &op_counts).unwrap();
        push_stream(&mut out, &syms).unwrap();
        push_stream(&mut out, &op_lens).unwrap();

        // Must fail at the record-count bound, before the op RLE resizes: any
        // other `Malformed` means decode walked into the huge run first.
        let err = decode(&out).expect_err("inflated n_records must be rejected");
        assert!(
            matches!(err, Error::Malformed("record count exceeds token stream")),
            "expected the record-count guard to fire, got {err:?}"
        );
    }

    /// A detected template must regenerate the names in ORIGINAL order exactly.
    fn assert_template_regenerates(names: &[&[u8]]) -> NameTemplate {
        let t = detect_template(names).expect("expected a regenerable template");
        // Round-trips through serialization.
        let t2 = NameTemplate::from_bytes(&t.to_bytes()).expect("template deser");
        assert_eq!(t, t2, "template serialization round-trip");
        for (i, name) in names.iter().enumerate() {
            assert_eq!(t.regenerate(i), *name, "regenerate[{i}] mismatch");
        }
        t
    }

    #[test]
    fn detect_sra_counter_names() {
        // The common SRA pattern: prefix + counter, repeated, + constant suffix.
        let names: Vec<Vec<u8>> = (1..=2000)
            .map(|i| format!("@DRR174812.{i} {i} length=150").into_bytes())
            .collect();
        let refs: Vec<&[u8]> = names.iter().map(|n| n.as_slice()).collect();
        let t = assert_template_regenerates(&refs);
        // Renumbering (regenerate at a different position) stays well-formed.
        assert_eq!(t.regenerate(0), b"@DRR174812.1 1 length=150");
        assert!(t.to_bytes().len() < 64, "template must be tiny");
    }

    #[test]
    fn renumber_template_is_a_bare_counter() {
        // The synthetic renumber template: position j -> "j+1", and it survives a
        // storage round-trip (it is written to / read from the container header).
        let t = NameTemplate::renumber();
        assert_eq!(t.regenerate(0), b"1");
        assert_eq!(t.regenerate(41), b"42");
        let rt = NameTemplate::from_bytes(&t.to_bytes()).unwrap();
        assert_eq!(rt, t);
        assert!(t.to_bytes().len() < 16, "renumber template must be tiny");
    }

    #[test]
    fn detect_zero_padded_counter() {
        let names: Vec<Vec<u8>> = (0..1500)
            .map(|i| format!("read{:06}", i).into_bytes())
            .collect();
        let refs: Vec<&[u8]> = names.iter().map(|n| n.as_slice()).collect();
        assert_template_regenerates(&refs);
    }

    #[test]
    fn reject_illumina_tile_coordinates() {
        // x/y coordinates are not a per-read counter — not regenerable.
        let names: &[&[u8]] = &[
            b"INST:1:FC:1:1101:1000:2000",
            b"INST:1:FC:1:1101:1005:2050",
            b"INST:1:FC:1:1101:1010:2100",
        ];
        assert!(detect_template(names).is_none());
    }

    #[test]
    fn reject_non_counter_number() {
        // A number that isn't first+index (jumps) is not regenerable.
        let names: &[&[u8]] = &[b"r.1", b"r.2", b"r.9"];
        assert!(detect_template(names).is_none());
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

    /// Apply byte substitutions, an optional 0xFF window (to drive a length/count
    /// field high), and an optional truncation to `data`.
    fn corrupt(
        mut data: Vec<u8>,
        subs: &[(usize, u8)],
        wipe: Option<(usize, usize)>,
        trunc: Option<usize>,
    ) -> Vec<u8> {
        for &(pos, val) in subs {
            if !data.is_empty() {
                let i = pos % data.len();
                data[i] = val;
            }
        }
        if let (Some((pos, w)), false) = (wipe, data.is_empty()) {
            let i = pos % data.len();
            for b in data.iter_mut().skip(i).take(w) {
                *b = 0xFF;
            }
        }
        if let Some(t) = trunc
            && !data.is_empty()
        {
            data.truncate(t % data.len());
        }
        data
    }

    proptest::proptest! {
        /// Encode valid names, corrupt the stream, then decode: a mutated name
        /// stream must never panic or abort — only return Ok/Err.
        #[test]
        fn decode_survives_mutation(
            names in proptest::collection::vec(
                proptest::collection::vec(proptest::sample::select(b"AB:._-0129 ".to_vec()), 0..40),
                0..60),
            subs in proptest::collection::vec((proptest::num::usize::ANY, proptest::num::u8::ANY), 0..12),
            wipe in proptest::option::of((proptest::num::usize::ANY, 1usize..9)),
            trunc in proptest::option::of(proptest::num::usize::ANY),
        ) {
            let refs: Vec<&[u8]> = names.iter().map(|n| n.as_slice()).collect();
            let enc = encode(&refs).expect("encode");
            let _ = decode(&corrupt(enc, &subs, wipe, trunc));
        }
    }

    // --- hardening: corrupt-input allocation guards ---

    #[test]
    fn from_bytes_rejects_huge_counter_pad() {
        // A template counter with an absurd zero-pad width would drive an
        // unbounded fill in `regenerate`; deserialization must reject it.
        let t = NameTemplate {
            columns: vec![TemplateColumn::Counter {
                start: 1,
                pad: 1 << 40,
            }],
        };
        assert!(NameTemplate::from_bytes(&t.to_bytes()).is_err());
    }

    #[test]
    fn read_op_rle_rejects_huge_run_without_aborting() {
        // One run of an enormous length in column 0. With an equally corrupt
        // `n_records` the length guard passes, so the fallible reserve is what
        // must reject it — pre-hardening this was an infallible resize that
        // aborted the process.
        let mut counts = Vec::new();
        write_varint(&mut counts, 1); // column 0: one run
        for _ in 0..MAX_COL {
            write_varint(&mut counts, 0); // columns 1..=MAX_COL: no runs
        }
        let mut lens = Vec::new();
        write_varint(&mut lens, 1u64 << 63); // run length past isize::MAX
        let mut data = Vec::new();
        push_stream(&mut data, &counts).unwrap();
        push_stream(&mut data, &[0u8]).unwrap(); // syms: one op byte
        push_stream(&mut data, &lens).unwrap();
        let mut r = Cursor::new(&data);
        let err = read_op_rle(&mut r, usize::MAX, usize::MAX).unwrap_err();
        assert!(matches!(err, Error::Malformed(_)), "got {err:?}");
    }

    #[test]
    fn take_rejects_overflowing_length_without_panicking() {
        // A stream whose length varint overflows `pos + n` must error, not panic
        // (debug) or wrap (release). Found by the cargo-fuzz `tokenizer` target.
        let src = [
            0x01, 0x01, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x37,
        ];
        assert!(decode(&src).is_err());
    }

    #[test]
    fn read_delta_planes_rejects_huge_count_without_aborting() {
        // A per-column value count taken from the untrusted `n_records` must size
        // its plane buffer fallibly rather than aborting.
        let mut data = Vec::new();
        write_varint(&mut data, 1); // width = 1 (non-zero)
        let mut r = Cursor::new(&data);
        let n_delta = [1usize << 61]; // (1<<61)*8 bytes overflows isize
        let err = read_delta_planes(&mut r, &n_delta).unwrap_err();
        assert!(matches!(err, Error::Malformed(_)), "got {err:?}");
    }

    #[test]
    fn read_op_rle_caps_aggregate_across_columns() {
        // #142 F3: each column is bounded by `n_records`, but the 65 columns
        // TOGETHER were not — 65 * n_records could allocate ~65x past the
        // per-stream cap (and then size the delta planes). Two columns, each a
        // single run of exactly `n_records` (so the per-column guard passes),
        // must trip the aggregate `max_out` cap on the second.
        let max_out = 100usize;
        let mut counts = Vec::new();
        write_varint(&mut counts, 1); // column 0: one run
        write_varint(&mut counts, 1); // column 1: one run
        for _ in 2..=MAX_COL {
            write_varint(&mut counts, 0); // remaining columns: empty
        }
        let mut lens = Vec::new();
        write_varint(&mut lens, max_out as u64); // col 0: fills to the cap
        write_varint(&mut lens, max_out as u64); // col 1: pushes the total past it
        let mut data = Vec::new();
        push_stream(&mut data, &counts).unwrap();
        push_stream(&mut data, &[0u8, 0u8]).unwrap(); // syms: two ops
        push_stream(&mut data, &lens).unwrap();
        let mut r = Cursor::new(&data);
        // n_records == max_out, so neither run trips the per-column guard; only
        // the aggregate cap can reject this.
        let err = read_op_rle(&mut r, max_out, max_out).unwrap_err();
        assert!(
            matches!(err, Error::Malformed("op columns exceed the stream bound")),
            "expected the aggregate cap, got {err:?}"
        );
    }

    proptest::proptest! {
        /// Arbitrary bytes fed straight to `decode` must never panic or abort —
        /// only return Ok/Err. Unlike `decode_survives_mutation` (which perturbs
        /// VALID encodings), this feeds pure garbage, which reaches header/length
        /// paths a mutation of a valid stream rarely does. Guards htscodecs-style
        /// crashes on adversarial name streams.
        #[test]
        fn decode_never_aborts_on_garbage(bytes in proptest::collection::vec(0u8..=255, 0..256)) {
            let _ = decode(&bytes);
        }
    }

    /// Names containing bytes >= 0x80 must round-trip byte-exactly. Guards
    /// htscodecs #105, a signed-char crash on the input `[0x80, 0x0a]`: the
    /// existing arbitrary-name generators only sample bytes < 0x80, so the
    /// high-bit path is otherwise untested.
    #[test]
    fn roundtrip_high_bit_name_bytes() {
        roundtrip(&[
            &[0x80u8, 0x0a],                   // the htscodecs #105 crasher
            &[0xffu8],                         // all bits set
            b"read\x80\xa5\xffname",           // high-bit bytes among ASCII
            &[0x00u8, 0x80, 0x7f, 0xff, 0x0a], // mix of low, boundary, and high
            b"plain",                          // an ordinary name alongside them
        ]);
    }

    /// Token-count boundary at `MAX_COL`. Guards the htscodecs 1.5.2 overflow when
    /// writing the token past the column cap. A name of single-character
    /// alternating letter/digit runs tokenizes to one column per character (each
    /// character flips digit/non-digit, so no runs merge), letting us hit exactly
    /// `MAX_COL`, `MAX_COL + 1`, and one past. Encode must not panic/overflow at
    /// the boundary, and each name must recover byte-exactly.
    #[test]
    fn roundtrip_at_max_col_boundary() {
        // `n` alternating single-char tokens: "a1a1a1…".
        let name_of = |n: usize| -> Vec<u8> {
            (0..n)
                .map(|i| if i % 2 == 0 { b'a' } else { b'1' })
                .collect()
        };
        // Sanity: the construction really yields `n` tokens (columns).
        assert_eq!(tokenize(&name_of(MAX_COL)).len(), MAX_COL);
        assert_eq!(tokenize(&name_of(MAX_COL + 1)).len(), MAX_COL + 1);

        for &n in &[MAX_COL - 1, MAX_COL, MAX_COL + 1] {
            let name = name_of(n);
            roundtrip(&[name.as_slice()]);
        }
    }
}
