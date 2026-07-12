//! fqzcomp-style quality-score context model.
//!
//! Each quality symbol is range-coded ([`fqxv_range`]) under a context built
//! from the two previous quality values and the position within the read — the
//! dominant signals in Illumina quality streams. One adaptive model per context.
//! Context resets at every read boundary, so [`encode`] takes per-read lengths.
//!
//! Lossy quality binning ([`QualityBinning`], Illumina 2/4/8-level) is applied
//! before modeling; the default is lossless.
//!
//! ```
//! use fqxv_fqzcomp::{encode, decode, QualityBinning};
//! let lens = [5u32, 3];
//! let quals = b"IIIII##F"; // two reads
//! let enc = encode(&lens, quals, QualityBinning::Lossless).unwrap();
//! let (out_lens, out_quals) = decode(&enc).unwrap();
//! assert_eq!(out_lens, lens);
//! assert_eq!(out_quals, quals);
//! ```

use fqxv_range::{Decoder, Encoder, SimpleModel};
use thiserror::Error;

/// Optional lossy quantization applied to quality scores before modeling.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum QualityBinning {
    /// No quantization — fully lossless (default).
    #[default]
    Lossless,
    /// Illumina 8-level binning.
    Bin8,
    /// Illumina 4-level binning.
    Bin4,
    /// 2-level (binary) binning.
    Bin2,
}

impl QualityBinning {
    fn tag(self) -> u8 {
        match self {
            QualityBinning::Lossless => 0,
            QualityBinning::Bin8 => 1,
            QualityBinning::Bin4 => 2,
            QualityBinning::Bin2 => 3,
        }
    }

    fn from_tag(t: u8) -> Result<Self> {
        Ok(match t {
            0 => QualityBinning::Lossless,
            1 => QualityBinning::Bin8,
            2 => QualityBinning::Bin4,
            3 => QualityBinning::Bin2,
            _ => return Err(Error::Malformed("unknown quality-binning tag")),
        })
    }

    /// Map a Phred+33 quality byte through the (possibly lossy) bin table.
    #[must_use]
    pub fn apply(self, byte: u8) -> u8 {
        if self == QualityBinning::Lossless {
            return byte;
        }
        let q = byte.saturating_sub(33);
        let b = match self {
            QualityBinning::Bin8 => match q {
                0..=1 => q,
                2..=9 => 6,
                10..=19 => 15,
                20..=24 => 22,
                25..=29 => 27,
                30..=34 => 33,
                35..=39 => 37,
                _ => 40,
            },
            QualityBinning::Bin4 => match q {
                0..=9 => 6,
                10..=24 => 18,
                25..=34 => 30,
                _ => 37,
            },
            QualityBinning::Bin2 => match q {
                0..=24 => 15,
                _ => 37,
            },
            QualityBinning::Lossless => q,
        };
        33 + b
    }
}

/// Errors returned by the quality codec.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// The compressed stream was malformed or truncated.
    #[error("malformed fqzcomp stream: {0}")]
    Malformed(&'static str),
    /// The quality alphabet exceeds what this codec models (64 symbols).
    #[error("quality alphabet too large ({0} > 64 symbols)")]
    AlphabetTooLarge(usize),
    /// The provided lengths do not sum to the quality-buffer size.
    #[error("read lengths ({lens}) do not match quality bytes ({quals})")]
    LengthMismatch {
        /// Sum of the provided read lengths.
        lens: usize,
        /// Number of quality bytes provided.
        quals: usize,
    },
    /// A code path that is not yet implemented in this scaffold.
    #[error("not yet implemented: {0}")]
    NotImplemented(&'static str),
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Max quality alphabet the model handles.
const QMAX: usize = 64;
/// Number of contexts: q1(6) | q2(6) | position-bucket(4) = 16 bits.
const N_CTX: usize = 1 << 16;
const FORMAT_VERSION: u8 = 0;

/// Build the context index from the two previous symbols and the position.
#[inline]
fn context(q1: u8, q2: u8, pos: usize) -> usize {
    let pb = (pos >> 3).min(15);
    (q1 as usize) | ((q2 as usize) << 6) | (pb << 12)
}

/// Encode per-read quality strings.
///
/// `lens` gives each read's quality length; `quals` is their concatenation.
/// `binning` optionally quantizes qualities before modeling (lossy).
pub fn encode(lens: &[u32], quals: &[u8], binning: QualityBinning) -> Result<Vec<u8>> {
    let total: usize = lens.iter().map(|&l| l as usize).sum();
    if total != quals.len() {
        return Err(Error::LengthMismatch {
            lens: total,
            quals: quals.len(),
        });
    }

    // Apply (optional) lossy binning, then map to a dense 0-based alphabet.
    let binned: Vec<u8> = quals.iter().map(|&b| binning.apply(b)).collect();
    let (qmin, qsize) = alphabet(&binned)?;

    let mut models = vec![SimpleModel::<QMAX>::new(); N_CTX];
    let mut enc = Encoder::new();
    let mut idx = 0usize;
    for &l in lens {
        let (mut q1, mut q2) = (0u8, 0u8);
        for pos in 0..l as usize {
            let sym = binned[idx] - qmin;
            idx += 1;
            models[context(q1, q2, pos)].encode(&mut enc, sym as usize);
            q2 = q1;
            q1 = sym;
        }
    }
    let payload = enc.finish();

    let mut out = Vec::with_capacity(16 + lens.len() + payload.len());
    out.push(FORMAT_VERSION);
    out.push(binning.tag());
    out.push(qmin);
    out.push(qsize as u8);
    write_lens(&mut out, lens);
    out.extend_from_slice(&payload);
    Ok(out)
}

/// Decode a stream produced by [`encode`], returning `(lengths, qualities)`.
/// In lossy modes the qualities are the binned values, not the originals.
pub fn decode(src: &[u8]) -> Result<(Vec<u32>, Vec<u8>)> {
    let mut r = ByteReader::new(src);
    if r.u8()? != FORMAT_VERSION {
        return Err(Error::Malformed("unsupported version"));
    }
    let _binning = QualityBinning::from_tag(r.u8()?)?;
    let qmin = r.u8()?;
    let qsize = r.u8()? as usize;
    if qsize > QMAX {
        return Err(Error::AlphabetTooLarge(qsize));
    }
    let lens = read_lens(&mut r)?;

    let mut models = vec![SimpleModel::<QMAX>::new(); N_CTX];
    let mut dec = Decoder::new(r.rest());
    let total: usize = lens.iter().map(|&l| l as usize).sum();
    let mut quals = Vec::with_capacity(total);
    for &l in &lens {
        let (mut q1, mut q2) = (0u8, 0u8);
        for pos in 0..l as usize {
            let sym = models[context(q1, q2, pos)].decode(&mut dec) as u8;
            if sym as usize >= qsize {
                return Err(Error::Malformed("decoded symbol outside alphabet"));
            }
            quals.push(sym + qmin);
            q2 = q1;
            q1 = sym;
        }
    }
    Ok((lens, quals))
}

/// Determine `(min_byte, alphabet_size)` over the quality bytes.
fn alphabet(quals: &[u8]) -> Result<(u8, usize)> {
    if quals.is_empty() {
        return Ok((0, 1));
    }
    let mut lo = u8::MAX;
    let mut hi = 0u8;
    for &b in quals {
        lo = lo.min(b);
        hi = hi.max(b);
    }
    let size = (hi - lo) as usize + 1;
    if size > QMAX {
        return Err(Error::AlphabetTooLarge(size));
    }
    Ok((lo, size))
}

// --- length stream (LEB128 varints, with a fixed-length fast path) -----------

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
    let mut lens = Vec::with_capacity(n);
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
            .ok_or(Error::Malformed("truncated header"))?;
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

    fn roundtrip(lens: &[u32], quals: &[u8], binning: QualityBinning) {
        let enc = encode(lens, quals, binning).expect("encode");
        let (out_lens, out_quals) = decode(&enc).expect("decode");
        assert_eq!(out_lens, lens, "lengths mismatch");
        let expect: Vec<u8> = quals.iter().map(|&b| binning.apply(b)).collect();
        assert_eq!(out_quals, expect, "qualities mismatch");
    }

    #[test]
    fn roundtrip_empty() {
        roundtrip(&[], b"", QualityBinning::Lossless);
    }

    #[test]
    fn roundtrip_two_reads() {
        roundtrip(&[5, 3], b"IIIII##F", QualityBinning::Lossless);
    }

    #[test]
    fn roundtrip_variable_lengths() {
        roundtrip(
            &[3, 1, 4, 1, 5],
            b"ABCDEFGHIJKLMN",
            QualityBinning::Lossless,
        );
    }

    #[test]
    fn roundtrip_binned() {
        let quals: Vec<u8> = (0..300).map(|i| b'!' + (i % 42) as u8).collect();
        for b in [
            QualityBinning::Bin8,
            QualityBinning::Bin4,
            QualityBinning::Bin2,
        ] {
            roundtrip(&[100, 100, 100], &quals, b);
        }
    }

    #[test]
    fn beats_raw_on_correlated_quality() {
        // Slowly drifting quality (like a real read): should compress well.
        let mut quals = Vec::new();
        let mut q = 30i32;
        let mut state = 0x2545_f491u32;
        for _ in 0..50_000 {
            state ^= state << 13;
            state ^= state >> 17;
            state ^= state << 5;
            q = (q + (state % 5) as i32 - 2).clamp(2, 40);
            quals.push(b'!' + q as u8);
        }
        let lens = vec![100u32; 500];
        let enc = encode(&lens, &quals, QualityBinning::Lossless).expect("encode");
        assert!(
            enc.len() < quals.len() / 2,
            "expected >2x on correlated quality, got {} -> {}",
            quals.len(),
            enc.len()
        );
    }

    proptest::proptest! {
        #[test]
        fn roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(33u8..=74, 0..50), 0..40)
        ) {
            let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
            let quals: Vec<u8> = reads.concat();
            roundtrip(&lens, &quals, QualityBinning::Lossless);
        }
    }
}
