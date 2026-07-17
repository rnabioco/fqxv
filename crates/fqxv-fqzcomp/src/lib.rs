//! fqzcomp-style quality-score context model.
//!
//! Each quality symbol is range-coded ([`fqxv_range`]) under a context built
//! from the previous three quality values (`q3` coarsely quantized), a running
//! "how noisy has this read been so far" delta counter, and the position within
//! the read — the dominant signals in Illumina quality streams (the same
//! features fqz_comp conditions on). One adaptive model per context. Context
//! resets at every read boundary, so [`encode`] takes per-read lengths.
//!
//! Lossy quality binning ([`QualityBinning`]) is applied before modeling; the
//! default is lossless. Three quantization tables are offered (exact ranges in
//! [`QualityBinning::apply`]):
//!
//! - **`Bin8`** — the standard Illumina 8-level scheme (HiSeq 2500/4000 and the
//!   Illumina "Reducing Whole-Genome Data Storage Footprint" whitepaper).
//!   Representatives `{6, 15, 22, 27, 33, 37, 40}` with Q0/Q1 preserved.
//! - **`Bin4`** — Illumina's current *documented* 4-level scheme (NovaSeq X /
//!   RTA4 control software v1.2): raw `0–2 → 2`, `3–17 → 12`, `18–29 → 24`,
//!   `30+ → 40`. This is deliberately the RTA4 table; Illumina does not publish
//!   the older NovaSeq 6000 / RTA3 cut points (whose representatives were
//!   `{2, 12, 23, 37}`), so Bin4 is *not* a no-op on RTA3-binned NovaSeq 6000
//!   data — it re-bins `23 → 24` and `37 → 40`.
//! - **`Bin2`** — a *custom* binary split with no Illumina equivalent: below
//!   Q25 → Q15, Q25+ → Q37.
//!
//! Binning is irreversible: only the binned values are entropy-coded, so decode
//! returns the binned qualities, never the originals.
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

use std::borrow::Cow;

use fqxv_bytes::{read_lens, write_lens, ReaderError};
use fqxv_range::{Decoder, Encoder, SimpleModel};
use thiserror::Error;

/// Optional lossy quantization applied to quality scores before modeling.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum QualityBinning {
    /// No quantization — fully lossless (default).
    #[default]
    Lossless,
    /// Standard Illumina 8-level binning (HiSeq); representatives
    /// `{6, 15, 22, 27, 33, 37, 40}`, Q0/Q1 preserved.
    Bin8,
    /// Illumina documented 4-level binning (NovaSeq X / RTA4);
    /// representatives `{2, 12, 24, 40}`.
    Bin4,
    /// Custom 2-level (binary) binning — no Illumina equivalent.
    Bin2,
    /// Long-read **Oxford Nanopore** 4-level binning; representatives
    /// `{3, 10, 18, 35}`. Cutpoints match CoLoRd's ONT table (validated to
    /// preserve downstream analysis), not Illumina's cycle-quality profile.
    BinOnt,
    /// Long-read **PacBio HiFi** 5-level binning; representatives
    /// `{3, 10, 18, 35}` plus Q93 kept exact. HiFi packs most bases near the top
    /// of the scale and encodes the max-quality symbol (Q93) with application
    /// meaning, so it is preserved as its own level (CoLoRd's HiFi table).
    BinHifi,
}

impl QualityBinning {
    fn tag(self) -> u8 {
        match self {
            QualityBinning::Lossless => 0,
            QualityBinning::Bin8 => 1,
            QualityBinning::Bin4 => 2,
            QualityBinning::Bin2 => 3,
            QualityBinning::BinOnt => 4,
            QualityBinning::BinHifi => 5,
        }
    }

    fn from_tag(t: u8) -> Result<Self> {
        Ok(match t {
            0 => QualityBinning::Lossless,
            1 => QualityBinning::Bin8,
            2 => QualityBinning::Bin4,
            3 => QualityBinning::Bin2,
            4 => QualityBinning::BinOnt,
            5 => QualityBinning::BinHifi,
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
            // NovaSeq X / RTA4 control software v1.2 documented 4-bin table.
            QualityBinning::Bin4 => match q {
                0..=2 => 2,
                3..=17 => 12,
                18..=29 => 24,
                _ => 40,
            },
            QualityBinning::Bin2 => match q {
                0..=24 => 15,
                _ => 37,
            },
            // CoLoRd ONT 4-level table (representatives 3/10/18/35).
            QualityBinning::BinOnt => match q {
                0..=6 => 3,
                7..=13 => 10,
                14..=25 => 18,
                _ => 35,
            },
            // CoLoRd HiFi 5-level table: as ONT, but the top Q93 symbol carries
            // application meaning and is preserved exactly rather than folded
            // into the 26+ bin.
            QualityBinning::BinHifi => match q {
                0..=6 => 3,
                7..=13 => 10,
                14..=25 => 18,
                26..=92 => 35,
                _ => 93,
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

/// Max quality alphabet the model handles. Covers the full Sanger FASTQ range
/// (Phred 0..=93, ASCII `!`..=`~`), so long-read (nanopore) data — whose quality
/// range routinely exceeds Illumina's ~40 levels and can span the whole scale —
/// is accepted rather than rejected. `context` masks its fields (below) so a
/// symbol beyond the old 64-cap can never index past `N_CTX`.
const QMAX: usize = 94;
/// Number of contexts: q1(6) | q2>>2(4) | delta(2) | q3>>4(2) | pos-bucket(4) =
/// 18 bits. `q3` (the third-previous quality, coarse) captures more of the local
/// quality trajectory. Keeping `q2` at 4 bits and growing to 18 bits beat the
/// 16-bit rebalance (coarsening `q2` to 2 bits lost more than `q3` added on
/// full-range data), at 4x the quality-model memory (~34 MB/block).
const N_CTX: usize = 1 << 18;
/// Saturating cap on the running delta counter (2 bits).
const DELTA_MAX: u8 = 3;
const FORMAT_VERSION: u8 = 1;
/// Absolute ceiling on a single decode's total quality bytes, guarding the
/// [`decode`] allocation/loop against a corrupt length header on memory-overcommit
/// systems (where `try_reserve` can't). Sized far above any real decode — quality
/// is one byte per base and a container row group caps at 256 MiB of sequence — so
/// it never rejects legitimate data; 16 GiB leaves a wide margin.
const MAX_DECODED_QUALS: usize = 1 << 34;

/// Position bucket: fine near the read start, then 32-wide buckets so the
/// low-quality tail of long reads keeps positional resolution. The old `pos>>3`
/// collapsed every position >= 120 into bucket 15; this saturates near 224.
#[inline]
fn pos_bucket(pos: usize) -> usize {
    if pos < 16 {
        pos >> 1 // 0..7, two positions per bucket
    } else {
        (8 + (pos >> 5)).min(15) // 8..15, 32 positions per bucket
    }
}

/// Build the context index from the previous three symbols, the running delta
/// counter, and the position.
///
/// Each field is masked to its bit width so the packed context stays within the
/// 18-bit `N_CTX` bound even when a symbol exceeds 63 (possible now that `QMAX`
/// spans the full Phred scale). For alphabets that fit the old 64-symbol cap the
/// masks are no-ops, so short-read output is byte-identical.
#[inline]
fn context(q1: u8, q2: u8, q3: u8, delta: u8, pos: usize) -> usize {
    ((q1 as usize) & 0x3F)                    // bits 0..5   (previous symbol)
        | (((q2 as usize >> 2) & 0xF) << 6)   // bits 6..9   (q2 coarsened)
        | (((delta as usize) & 0x3) << 10)    // bits 10..11
        | (((q3 as usize >> 4) & 0x3) << 12)  // bits 12..13 (q3 coarse)
        | ((pos_bucket(pos) & 0xF) << 14) // bits 14..17
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

    // Apply (optional) lossy binning, then map to a dense 0-based alphabet. On
    // the lossless default `apply` is the identity, so borrow `quals` directly
    // instead of allocating and copying a block-sized duplicate.
    let binned: Cow<[u8]> = if binning == QualityBinning::Lossless {
        Cow::Borrowed(quals)
    } else {
        Cow::Owned(quals.iter().map(|&b| binning.apply(b)).collect())
    };
    let (syms, dense) = dense_alphabet(&binned)?;
    let qmin = syms[0];
    let k = syms.len();

    // Models are sized to the symbols that actually occur (`k`), not the 0..QMAX
    // capacity — see `SimpleModel::with_active`. Context features stay on the
    // original Phred scale (`cv = b - qmin`); only the coded symbol is the dense
    // index (`dv`).
    let mut models = vec![SimpleModel::<QMAX>::with_active(k); N_CTX];
    let mut enc = Encoder::new();
    // Walk `binned` as a cursor rather than a running `idx`, so the per-quality
    // access iterates a slice directly and carries no bounds check. `total`
    // (checked above) equals `binned.len()`, so every `split_at` is in range.
    let mut rest: &[u8] = &binned;
    for &l in lens {
        let (read, tail) = rest.split_at(l as usize);
        rest = tail;
        let (mut q1, mut q2, mut q3) = (0u8, 0u8, 0u8);
        let mut delta = 0u8;
        for (pos, &b) in read.iter().enumerate() {
            let cv = b - qmin; // context value (original Phred scale)
            let dv = dense[b as usize]; // dense coding symbol
            let c = context(q1, q2, q3, delta, pos);
            debug_assert!(c < N_CTX);
            // SAFETY: `context` packs into 18 bits, so `c < N_CTX == models.len()`.
            unsafe { models.get_unchecked_mut(c) }.encode(&mut enc, dv as usize);
            if pos > 0 && cv != q1 {
                delta = (delta + 1).min(DELTA_MAX);
            }
            q3 = q2;
            q2 = q1;
            q1 = cv;
        }
    }
    let payload = enc.finish();

    let mut out = Vec::with_capacity(16 + k + lens.len() + payload.len());
    out.push(FORMAT_VERSION);
    out.push(binning.tag());
    out.push(k as u8);
    out.extend_from_slice(&syms);
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
    let k = r.u8()? as usize;
    if k == 0 || k > QMAX {
        return Err(Error::AlphabetTooLarge(k));
    }
    let syms = r.take(k)?.to_vec();
    let qmin = syms[0];
    let lens = read_lens(&mut r)?;

    let mut models = vec![SimpleModel::<QMAX>::with_active(k); N_CTX];
    // Checked sum: a malformed stream can declare lengths whose total wraps
    // `usize`, which would under-allocate `quals` and then over-push.
    let total = lens
        .iter()
        .try_fold(0usize, |acc, &l| acc.checked_add(l as usize))
        .ok_or(Error::Malformed("total length overflows usize"))?;
    let mut dec = Decoder::new(r.rest());
    // Decompression-bomb guard: bound the *allocation*, not the ratio. A single
    // repeated quality symbol codes to almost nothing (a 1-symbol alphabet costs
    // ~0 bits/symbol), so there is no finite "symbols per compressed byte" bound —
    // any ratio cap would reject legitimately compressible constant/low-entropy
    // quality (including the lossy `--quality-bin` modes). The container caps a
    // block's read count and cross-checks it against the decoded content digest,
    // so a lying `total` can't turn into wrong output.
    //
    // Reserving `total` fallibly rejects a hostile length on systems that refuse
    // the allocation — but memory overcommit (macOS always, Linux by default)
    // accepts a multi-terabyte reservation and then stalls in the decode loop
    // below. So first reject any `total` past an absolute ceiling. This is not a
    // ratio (which would lose data): no real single decode approaches it — quality
    // is one byte per base and a container block caps at `MAX_BLOCK_SEQ_BYTES`.
    if total > MAX_DECODED_QUALS {
        return Err(Error::Malformed("declared total length exceeds maximum"));
    }
    let mut quals = Vec::new();
    quals
        .try_reserve(total)
        .map_err(|_| Error::Malformed("declared total length too large to allocate"))?;
    for &l in &lens {
        let (mut q1, mut q2, mut q3) = (0u8, 0u8, 0u8);
        let mut delta = 0u8;
        for pos in 0..l as usize {
            let c = context(q1, q2, q3, delta, pos);
            debug_assert!(c < N_CTX);
            // SAFETY: `context` packs into 18 bits, so `c < N_CTX == models.len()`.
            let dv = unsafe { models.get_unchecked_mut(c) }.decode(&mut dec);
            // `with_active(k)` only ever assigns probability to symbols `0..k`, so
            // a well-formed stream decodes within `syms`; guard anyway.
            let b = *syms
                .get(dv)
                .ok_or(Error::Malformed("decoded symbol outside alphabet"))?;
            let cv = b - qmin;
            quals.push(b);
            if pos > 0 && cv != q1 {
                delta = (delta + 1).min(DELTA_MAX);
            }
            q3 = q2;
            q2 = q1;
            q1 = cv;
        }
    }
    Ok((lens, quals))
}

/// The set of quality values actually present, as a compact coding alphabet.
///
/// Returns the sorted distinct bytes (`syms`) and a 256-entry map from byte to
/// its dense index in `syms`. The coder sizes its per-context models to
/// `syms.len()` — only the values that occur — so a stream using few of the
/// possible Phred levels (e.g. NovaSeq's 4) pays nothing for the ones it never
/// uses. `syms[0]` is the minimum byte; it doubles as the context origin so the
/// context features stay on the original (spread-out) Phred scale rather than the
/// dense indices, which would collapse together under the context's bit-shifts.
fn dense_alphabet(quals: &[u8]) -> Result<(Vec<u8>, [u8; 256])> {
    let mut present = [false; 256];
    for &b in quals {
        present[b as usize] = true;
    }
    let syms: Vec<u8> = (0..=255u8).filter(|&b| present[b as usize]).collect();
    if syms.is_empty() {
        // No symbols (empty input): a single dummy so models are well-formed.
        return Ok((vec![0], [0u8; 256]));
    }
    if syms.len() > QMAX {
        return Err(Error::AlphabetTooLarge(syms.len()));
    }
    let mut map = [0u8; 256];
    for (i, &b) in syms.iter().enumerate() {
        map[b as usize] = i as u8;
    }
    Ok((syms, map))
}

// --- length stream (LEB128 varints, with a fixed-length fast path) -----------

/// Shared byte cursor specialized to this crate's [`Error`].
type ByteReader<'a> = fqxv_bytes::Reader<'a, Error>;

impl ReaderError for Error {
    fn truncated() -> Self {
        Error::Malformed("truncated header")
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
    // Only the corruption tests hand-build streams; `write_lens` moved to
    // fqxv-bytes, so this varint helper is no longer used outside tests.
    use fqxv_bytes::write_varint;

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
    fn roundtrip_large_constant_quality() {
        // Constant quality compresses to almost nothing, so total / payload_len
        // far exceeds any fixed ratio — the old "declared length exceeds payload
        // capacity" guard wrongly rejected this on decode, so `compress` produced
        // an archive `decompress` refused. Regression: it must round-trip.
        let n = 100_000usize;
        let read_len = 100u32;
        let lens = vec![read_len; n];
        let quals = vec![b'I'; n * read_len as usize]; // 10M identical symbols
        let enc = encode(&lens, &quals, QualityBinning::Lossless).expect("encode");
        assert!(
            enc.len() < 100_000,
            "constant quality must compress tiny (got {} bytes for 10M symbols)",
            enc.len()
        );
        let (out_lens, out_quals) = decode(&enc).expect("high compression must not be rejected");
        assert_eq!(out_lens, lens);
        assert_eq!(out_quals, quals);
    }

    #[test]
    fn decode_rejects_length_bomb_without_aborting() {
        // A tiny stream declaring an enormous fixed-length total must fail the
        // fallible reserve, not abort (the removed ratio guard used to catch this;
        // `try_reserve` on `total` still does).
        let mut src = vec![FORMAT_VERSION, 0, 1, b'I']; // version, lossless, k=1, symbol
        write_varint(&mut src, 1000); // n = 1000 reads
        src.push(1); // fixed-length flag
        write_varint(&mut src, u64::from(u32::MAX)); // f -> total ~4.3e12 bytes
        assert!(matches!(decode(&src), Err(Error::Malformed(_))));
    }

    #[test]
    fn roundtrip_binned() {
        let quals: Vec<u8> = (0..300).map(|i| b'!' + (i % 42) as u8).collect();
        for b in [
            QualityBinning::Bin8,
            QualityBinning::Bin4,
            QualityBinning::Bin2,
            QualityBinning::BinOnt,
            QualityBinning::BinHifi,
        ] {
            roundtrip(&[100, 100, 100], &quals, b);
        }
    }

    #[test]
    fn roundtrip_binned_hifi_high_q() {
        // HiFi packs quality near the top of the scale and keeps Q93 exact;
        // exercise the wide-alphabet path with values through Q93.
        let quals: Vec<u8> = (0..300).map(|i| b'!' + (i % 94) as u8).collect();
        for b in [QualityBinning::BinHifi, QualityBinning::BinOnt] {
            roundtrip(&[150, 150], &quals, b);
        }
    }

    #[test]
    fn bin_tables_map_expected_values() {
        // Phred value q -> Phred+33 byte.
        let ch = |q: u8| q + 33;
        // Bin8 (standard Illumina 8-level): Q0/Q1 preserved, then bands.
        for (q, want) in [
            (0, 0),
            (1, 1),
            (5, 6),
            (9, 6),
            (12, 15),
            (22, 22),
            (27, 27),
            (33, 33),
            (37, 37),
            (41, 40),
        ] {
            assert_eq!(QualityBinning::Bin8.apply(ch(q)), ch(want), "Bin8 q={q}");
        }
        // Bin4 (NovaSeq X / RTA4 v1.2): 0-2->2, 3-17->12, 18-29->24, 30+->40.
        for (q, want) in [
            (0, 2),
            (2, 2),
            (3, 12),
            (17, 12),
            (18, 24),
            (29, 24),
            (30, 40),
            (41, 40),
        ] {
            assert_eq!(QualityBinning::Bin4.apply(ch(q)), ch(want), "Bin4 q={q}");
        }
        // The four RTA4 representatives are fixed points of Bin4.
        for q in [2, 12, 24, 40] {
            assert_eq!(
                QualityBinning::Bin4.apply(ch(q)),
                ch(q),
                "Bin4 fixed point q={q}"
            );
        }
        // Bin2 (custom binary): <Q25 -> Q15, Q25+ -> Q37.
        for (q, want) in [(0, 15), (24, 15), (25, 37), (41, 37)] {
            assert_eq!(QualityBinning::Bin2.apply(ch(q)), ch(want), "Bin2 q={q}");
        }
        // BinOnt (CoLoRd ONT 4-level): 0-6->3, 7-13->10, 14-25->18, 26+->35.
        for (q, want) in [
            (0, 3),
            (6, 3),
            (7, 10),
            (13, 10),
            (14, 18),
            (25, 18),
            (26, 35),
            (93, 35),
        ] {
            assert_eq!(
                QualityBinning::BinOnt.apply(ch(q)),
                ch(want),
                "BinOnt q={q}"
            );
        }
        // BinHifi (CoLoRd HiFi 5-level): as ONT below 26, 26-92->35, Q93 exact.
        for (q, want) in [(0, 3), (6, 3), (14, 18), (26, 35), (92, 35), (93, 93)] {
            assert_eq!(
                QualityBinning::BinHifi.apply(ch(q)),
                ch(want),
                "BinHifi q={q}"
            );
        }
        // Lossless is the identity.
        for q in 0..=42 {
            assert_eq!(
                QualityBinning::Lossless.apply(ch(q)),
                ch(q),
                "Lossless q={q}"
            );
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

    #[test]
    fn roundtrip_wide_nanopore_alphabet() {
        // Nanopore quality spans far more than Illumina's ~40 levels. Exercise
        // the whole '!'..='~' alphabet (Phred 0..=93, 94 symbols) across long
        // reads — the old 64-symbol cap rejected this outright.
        let quals: Vec<u8> = (0..12_000u32).map(|i| b'!' + (i % 94) as u8).collect();
        roundtrip(&[3000, 3000, 3000, 3000], &quals, QualityBinning::Lossless);
    }

    #[test]
    fn accepts_full_sanger_range_rejects_beyond() {
        // '!' (33) .. '~' (126) is exactly 94 symbols — the widest valid FASTQ
        // quality alphabet — and must encode.
        let full: Vec<u8> = (b'!'..=b'~').collect();
        assert_eq!(full.len(), 94);
        encode(&[full.len() as u32], &full, QualityBinning::Lossless)
            .expect("full Sanger range encodes");
        // A contiguous span of 95 distinct bytes exceeds the model and must be a
        // clean AlphabetTooLarge error, not a panic or silent corruption.
        let over: Vec<u8> = (33u8..33 + 95).collect();
        assert!(matches!(
            encode(&[over.len() as u32], &over, QualityBinning::Lossless),
            Err(Error::AlphabetTooLarge(95))
        ));
    }

    fn push_varint(out: &mut Vec<u8>, mut v: u64) {
        loop {
            let b = (v & 0x7f) as u8;
            v >>= 7;
            if v == 0 {
                out.push(b);
                break;
            }
            out.push(b | 0x80);
        }
    }

    // A hostile length header must produce a clean `Err`, never abort the
    // process on an impossible allocation. Regression for the DoS where a
    // ~13-byte stream requested hundreds of petabytes via `Vec::with_capacity`.
    #[test]
    fn rejects_huge_length_count() {
        let mut buf = vec![FORMAT_VERSION, 0, 1, 0]; // version, binning, k=1, syms=[0]
        push_varint(&mut buf, u64::MAX >> 8); // n: absurd length count
        buf.push(1); // fixed = true -> resize(n, f) path
        push_varint(&mut buf, 100); // f
        assert!(matches!(decode(&buf), Err(Error::Malformed(_))));
    }

    #[test]
    fn rejects_huge_total_length() {
        let mut buf = vec![FORMAT_VERSION, 0, 1, 0]; // version, binning, k=1, syms=[0]
        push_varint(&mut buf, 1000); // n reads
        buf.push(1); // fixed = true
        push_varint(&mut buf, u32::MAX as u64); // each read u32::MAX long
        assert!(matches!(decode(&buf), Err(Error::Malformed(_))));
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

        // Full Sanger quality range ('!'..='~'), as long-read basecallers emit.
        #[test]
        fn roundtrip_wide_alphabet_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(33u8..=126, 0..80), 0..30)
        ) {
            let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
            let quals: Vec<u8> = reads.concat();
            roundtrip(&lens, &quals, QualityBinning::Lossless);
        }

        // Arbitrary bytes must never panic or abort the decoder — only Ok/Err.
        #[test]
        fn decode_never_aborts_on_garbage(bytes in proptest::collection::vec(0u8..=255, 0..256)) {
            let _ = decode(&bytes);
        }
    }
}
