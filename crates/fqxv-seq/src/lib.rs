//! Nucleotide sequence coding via an order-k adaptive context model.
//!
//! Each base is one of A/C/G/T (2 bits); the model conditions on the previous
//! `k` bases and is range-coded ([`fqxv_range`]). Non-ACGT bytes (`N`, IUPAC
//! codes, lowercase) are recorded in an exception list and restored verbatim on
//! decode, so the codec is byte-exact. Context resets at every read boundary.
//!
//! This is the sequence path for reads that are *not* reordered; the reordered
//! path lives in `fqxv-reorder`. Because the model is adaptive it uses the range
//! coder, not rANS (whose reverse encode can't carry adaptive state).
//!
//! ```
//! use fqxv_seq::{encode, decode};
//! let lens = [4u32, 5];
//! let seq = b"ACGTACGTN"; // note the trailing N
//! let enc = encode(&lens, seq, 6).unwrap();
//! let (out_lens, out_seq) = decode(&enc).unwrap();
//! assert_eq!(out_lens, lens);
//! assert_eq!(out_seq, seq);
//! ```

use fqxv_range::{Decoder, Encoder, SimpleModel};
use thiserror::Error;

/// Errors returned by the sequence codec.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// The compressed stream was malformed or truncated.
    #[error("malformed sequence stream: {0}")]
    Malformed(&'static str),
    /// The provided lengths do not sum to the sequence-buffer size.
    #[error("read lengths ({lens}) do not match sequence bytes ({seq})")]
    LengthMismatch {
        /// Sum of the provided read lengths.
        lens: usize,
        /// Number of sequence bytes provided.
        seq: usize,
    },
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;

const FORMAT_VERSION: u8 = 0;
/// Largest context order (4^11 contexts ≈ 4.2M models).
const MAX_ORDER: usize = 11;

/// byte -> 2-bit symbol, 255 for non-ACGT (an exception).
const BASE_LUT: [u8; 256] = base_lut();
const SYM2BASE: [u8; 4] = *b"ACGT";

const fn base_lut() -> [u8; 256] {
    let mut t = [255u8; 256];
    t[b'A' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'T' as usize] = 3;
    t
}

/// Encode per-read sequences with an order-`order` context model.
pub fn encode(lens: &[u32], seq: &[u8], order: usize) -> Result<Vec<u8>> {
    let total: usize = lens.iter().map(|&l| l as usize).sum();
    if total != seq.len() {
        return Err(Error::LengthMismatch {
            lens: total,
            seq: seq.len(),
        });
    }
    let k = order.clamp(1, MAX_ORDER);
    let ctx_mask = (1usize << (2 * k)) - 1;

    let mut models = vec![SimpleModel::<4>::new(); ctx_mask + 1];
    let mut enc = Encoder::new();
    let mut exceptions: Vec<(usize, u8)> = Vec::new();
    let mut idx = 0usize;
    for &l in lens {
        let mut ctx = 0usize;
        for _ in 0..l {
            let byte = seq[idx];
            let raw = BASE_LUT[byte as usize];
            let sym = if raw == 255 {
                exceptions.push((idx, byte));
                0
            } else {
                raw as usize
            };
            models[ctx].encode(&mut enc, sym);
            ctx = ((ctx << 2) | sym) & ctx_mask;
            idx += 1;
        }
    }
    let payload = enc.finish();

    let mut out = Vec::with_capacity(16 + lens.len() + exceptions.len() * 2 + payload.len());
    out.push(FORMAT_VERSION);
    out.push(k as u8);
    write_lens(&mut out, lens);
    write_exceptions(&mut out, &exceptions);
    out.extend_from_slice(&payload);
    Ok(out)
}

/// Decode a stream produced by [`encode`], returning `(lengths, sequence)`.
pub fn decode(src: &[u8]) -> Result<(Vec<u32>, Vec<u8>)> {
    let mut r = ByteReader::new(src);
    if r.u8()? != FORMAT_VERSION {
        return Err(Error::Malformed("unsupported version"));
    }
    let k = r.u8()? as usize;
    if !(1..=MAX_ORDER).contains(&k) {
        return Err(Error::Malformed("order out of range"));
    }
    let ctx_mask = (1usize << (2 * k)) - 1;
    let lens = read_lens(&mut r)?;
    let exceptions = read_exceptions(&mut r)?;

    let mut models = vec![SimpleModel::<4>::new(); ctx_mask + 1];
    let mut dec = Decoder::new(r.rest());
    let total: usize = lens.iter().map(|&l| l as usize).sum();
    let mut seq = Vec::with_capacity(total);
    for &l in &lens {
        let mut ctx = 0usize;
        for _ in 0..l {
            let sym = models[ctx].decode(&mut dec);
            seq.push(SYM2BASE[sym]);
            ctx = ((ctx << 2) | sym) & ctx_mask;
        }
    }
    // Restore non-ACGT bytes verbatim.
    for (pos, byte) in exceptions {
        *seq.get_mut(pos)
            .ok_or(Error::Malformed("exception out of range"))? = byte;
    }
    Ok((lens, seq))
}

// --- length + exception side streams (LEB128 varints) ------------------------

fn write_lens(out: &mut Vec<u8>, lens: &[u32]) {
    write_varint(out, lens.len() as u64);
    let fixed = lens.first().is_some_and(|&f| lens.iter().all(|&l| l == f));
    out.push(u8::from(fixed));
    if fixed {
        if let Some(&f) = lens.first() {
            write_varint(out, u64::from(f));
        }
    } else {
        for &l in lens {
            write_varint(out, u64::from(l));
        }
    }
}

fn read_lens(r: &mut ByteReader<'_>) -> Result<Vec<u32>> {
    let n = r.varint()? as usize;
    let fixed = r.u8()? != 0;
    let mut lens = Vec::with_capacity(n.min(1 << 20));
    if fixed {
        if n > 0 {
            let f = r.varint()? as u32;
            lens.resize(n, f);
        }
    } else {
        for _ in 0..n {
            lens.push(r.varint()? as u32);
        }
    }
    Ok(lens)
}

fn write_exceptions(out: &mut Vec<u8>, exceptions: &[(usize, u8)]) {
    write_varint(out, exceptions.len() as u64);
    let mut prev = 0usize;
    for &(pos, byte) in exceptions {
        write_varint(out, (pos - prev) as u64);
        out.push(byte);
        prev = pos;
    }
}

fn read_exceptions(r: &mut ByteReader<'_>) -> Result<Vec<(usize, u8)>> {
    let n = r.varint()? as usize;
    let mut v = Vec::with_capacity(n.min(1 << 20));
    let mut pos = 0usize;
    for _ in 0..n {
        pos += r.varint()? as usize;
        let byte = r.u8()?;
        v.push((pos, byte));
    }
    Ok(v)
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

struct ByteReader<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> ByteReader<'a> {
    fn new(buf: &'a [u8]) -> Self {
        ByteReader { buf, pos: 0 }
    }
    fn u8(&mut self) -> Result<u8> {
        let b = *self
            .buf
            .get(self.pos)
            .ok_or(Error::Malformed("truncated stream"))?;
        self.pos += 1;
        Ok(b)
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
    fn rest(&self) -> &'a [u8] {
        &self.buf[self.pos..]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn roundtrip(lens: &[u32], seq: &[u8], order: usize) {
        let enc = encode(lens, seq, order).expect("encode");
        let (out_lens, out_seq) = decode(&enc).expect("decode");
        assert_eq!(out_lens, lens, "lengths mismatch");
        assert_eq!(out_seq, seq, "sequence mismatch (order {order})");
    }

    #[test]
    fn roundtrip_empty() {
        roundtrip(&[], b"", 4);
    }

    #[test]
    fn roundtrip_basic() {
        roundtrip(&[4, 5], b"ACGTACGTN", 6);
    }

    #[test]
    fn roundtrip_with_exceptions() {
        roundtrip(&[10], b"ACGTNNRYAC", 4);
    }

    #[test]
    fn roundtrip_orders() {
        let seq: Vec<u8> = (0..1000).map(|i| SYM2BASE[(i * 7 + i / 3) % 4]).collect();
        for order in [1usize, 2, 4, 8, 11] {
            roundtrip(&[100; 10], &seq, order);
        }
    }

    #[test]
    fn repetitive_compresses() {
        // A repeated motif should compress well at a modest order.
        let unit = b"ACGTACGGTA";
        let seq: Vec<u8> = unit.iter().cycle().take(50_000).copied().collect();
        let enc = encode(&[50_000], &seq, 8).expect("encode");
        assert!(
            enc.len() < seq.len() / 4,
            "expected strong compression on repeats, got {} -> {}",
            seq.len(),
            enc.len()
        );
    }

    proptest::proptest! {
        #[test]
        fn roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(
                    proptest::sample::select(b"ACGTNacgtRYKM".to_vec()), 0..60),
                0..40),
            order in 1usize..=11,
        ) {
            let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
            let seq: Vec<u8> = reads.concat();
            roundtrip(&lens, &seq, order);
        }
    }
}
