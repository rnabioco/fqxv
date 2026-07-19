//! fqzcomp-style quality-score context model.
//!
//! Each quality symbol is range-coded ([`fqxv_range`]) under a per-context
//! adaptive model; the context resets at every read boundary, so [`encode`] takes
//! per-read lengths. Two context modes are chosen automatically by mean read
//! length (a self-describing byte in the stream header records which):
//!
//! - **Position context** (short reads): the previous three quality values (`q3`
//!   coarsely quantized), a running "how noisy has this read been so far" delta
//!   counter, and the position within the read — the dominant signals in Illumina
//!   quality streams (the same features fqz_comp conditions on).
//! - **Sequence context** (long reads, via [`encode_seq`]): the two previous
//!   qualities plus the current base, the next base, and the homopolymer
//!   run-length — where HiFi/ONT quality actually lives. This drops the
//!   position/delta features (useless on long reads) and requires the decoded
//!   sequence at decode time (see [`decode_seq`]/[`needs_sequence`]).
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

use fqxv_bytes::{ReaderError, read_lens, write_lens};
use fqxv_range::{Decoder, Encoder, SimpleModel};
use thiserror::Error;

mod binmix;

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
/// Number of contexts, sized to the wider of the two modes. Position context:
/// q1(6) | q2>>2(4) | delta(2) | q3>>4(2) | pos-bucket(4) = 18 bits. `q3` (the
/// third-previous quality, coarse) captures more of the local quality trajectory;
/// keeping `q2` at 4 bits and growing to 18 bits beat the 16-bit rebalance
/// (coarsening `q2` to 2 bits lost more than `q3` added on full-range data), at 4x
/// the quality-model memory (~34 MB/block). Sequence context ([`context_lr`]) packs
/// into 16 bits, well inside the same table.
const N_CTX: usize = 1 << 18;
/// Saturating cap on the running delta counter (2 bits).
const DELTA_MAX: u8 = 3;
/// Stream format version. Bumped to 2 when the header gained its context-mode byte
/// (see `MODE_POS`/`MODE_SEQ`); a v1 stream has no such byte, so rejecting it
/// here is cleaner than silently misreading the old layout.
const FORMAT_VERSION: u8 = 2;

/// Context-model selector, stored as a byte in the stream header (after the
/// binning tag). The decoder reads it to know which context to reconstruct — and,
/// for `MODE_SEQ`, that it must be handed the block's decoded sequence.
///
/// - `MODE_POS`: the original sequence-blind context (q1, q2, q3, delta, pos).
///   Tuned for short reads, whose quality tracks position; needs no sequence, so
///   the container decodes it in parallel with the sequence stream.
/// - `MODE_SEQ`: a long-read context (q1, q2, current base, next base,
///   homopolymer run-length) that drops the position/delta features — useless on
///   long reads — for the base identity that HiFi/ONT quality actually follows.
///   Requires the decoded sequence, so the container serializes seq → qual.
const MODE_POS: u8 = 0;
/// Single-context long-read quality model (base + next base + homopolymer run).
/// An early long-read design, superseded by `MODE_SEQ_BINMIX` and no longer
/// emitted; still decodable. (Mode 2 was a k-way softmax mixer, since retired.)
const MODE_SEQ: u8 = 1;
/// Long-read **binary-decomposition** logistic-mixing quality model ([`binmix`]):
/// codes each quality as `ceil(log2 k)` bit-tree decisions with binary logistic
/// mixing across coarse/mid/rich context tiers. Beats a single packed context on
/// ratio — a per-block adaptive model can't exploit a richer *single* context but
/// *can* mix a well-trained coarse model with a sparse rich one — and runs faster
/// than the single-context coder. The **default** long-read quality mode; needs the
/// decoded bases like `MODE_SEQ`.
const MODE_SEQ_BINMIX: u8 = 3;

/// Mean read length (bases) above which [`encode_seq`] selects `MODE_SEQ`. Long
/// reads (HiFi ~15 kb, ONT ~10 kb) clear this comfortably; Illumina (≤250 bp)
/// never does, so short-read archives keep the position context and the parallel
/// decode. A middle ground here would only ever misclassify synthetic data.
const SEQ_MODE_MIN_MEAN_LEN: usize = 500;
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

/// 2-bit code for a base (A/C/G/T → 0/1/2/3; anything else, including `N` and the
/// end-of-read sentinel, folds to 0). Only feeds the quality context, so the fold
/// is harmless — encode and decode compute it identically.
#[inline]
fn base_code(b: u8) -> usize {
    match b {
        b'A' | b'a' => 0,
        b'C' | b'c' => 1,
        b'G' | b'g' => 2,
        b'T' | b't' => 3,
        _ => 0,
    }
}

/// Long-read (`MODE_SEQ`) context: condition each quality byte on the two
/// previous qualities, the current and next base, and the homopolymer run-length
/// ending at this position. This is where HiFi/ONT quality lives — base identity
/// and homopolymer runs, not read position — so we spend the 18 bits there and
/// drop the position/delta features [`context`] uses for short reads.
///
/// Packs into 16 bits (`< N_CTX`): q1 (coarse /2, 6 bits) | q2 (coarse /8, 3) |
/// base (2) | next base (2) | run-length capped at 7 (3).
#[inline]
fn context_lr(q1: u8, q2: u8, base: usize, next: usize, hp_run: usize) -> usize {
    ((q1 as usize >> 1) & 0x3F)              // bits 0..5   previous quality (coarse)
        | (((q2 as usize >> 3) & 0x7) << 6)  // bits 6..8   second-previous quality (coarse)
        | ((base & 0x3) << 9)                // bits 9..10  current base
        | ((next & 0x3) << 11)               // bits 11..12 next base
        | ((hp_run.min(7)) << 13) // bits 13..15 homopolymer run-length
}

/// Size-class dispatch for the per-context quality models.
///
/// `SimpleModel<N>` is `[u16; N] + u32`, so the `N_CTX`-entry table costs
/// `N_CTX * (2N + 4)` bytes — 48 MiB at `N = QMAX = 94`. Only the first `k`
/// symbols are ever active (`with_active`), and the coder never touches a slot
/// past `k`: `encode` reads `freq[..=sym]` with `sym < k`, and `decode`'s scan
/// stops at `k-1` because `cum + freq[k-1] == tot > target`. So any `N >= k`
/// produces a byte-identical stream; picking the smallest class that fits `k`
/// shrinks the table (and the cache footprint) without touching the format.
macro_rules! by_size_class {
    ($k:expr, $f:ident, $($arg:expr),* $(,)?) => {
        match $k {
            0..=4 => $f::<4>($($arg),*),
            5..=8 => $f::<8>($($arg),*),
            9..=16 => $f::<16>($($arg),*),
            17..=32 => $f::<32>($($arg),*),
            33..=64 => $f::<64>($($arg),*),
            _ => $f::<QMAX>($($arg),*),
        }
    };
}

fn encode_payload<const NM: usize>(
    lens: &[u32],
    binned: &[u8],
    seq: Option<&[u8]>,
    dense: &[u8; 256],
    qmin: u8,
    k: usize,
) -> Vec<u8> {
    let mut models = vec![SimpleModel::<NM>::with_active(k); N_CTX];
    let mut enc = Encoder::new();
    let mut rest: &[u8] = binned;
    let seq_mode = seq.is_some();
    let mut srest: &[u8] = seq.unwrap_or(&[]);
    for &l in lens {
        let (read, tail) = rest.split_at(l as usize);
        rest = tail;
        // In `MODE_SEQ` the bases run in lockstep with the qualities; slice this
        // read's bases off the parallel stream. `MODE_POS` never touches `sread`.
        let sread: &[u8] = if seq_mode {
            let (sr, st) = srest.split_at(l as usize);
            srest = st;
            sr
        } else {
            &[]
        };
        let (mut q1, mut q2, mut q3) = (0u8, 0u8, 0u8);
        let mut delta = 0u8;
        let mut prev_base = u8::MAX;
        let mut run = 0usize;
        for (pos, &b) in read.iter().enumerate() {
            let cv = b - qmin;
            let dv = dense[b as usize];
            let c = if seq_mode {
                let base = sread[pos];
                let next = sread.get(pos + 1).copied().unwrap_or(u8::MAX);
                run = if base == prev_base { run + 1 } else { 1 };
                prev_base = base;
                context_lr(q1, q2, base_code(base), base_code(next), run)
            } else {
                context(q1, q2, q3, delta, pos)
            };
            debug_assert!(c < N_CTX);
            // SAFETY: both contexts pack into ≤18 bits, so `c < N_CTX == models.len()`.
            unsafe { models.get_unchecked_mut(c) }.encode(&mut enc, dv as usize);
            if pos > 0 && cv != q1 {
                delta = (delta + 1).min(DELTA_MAX);
            }
            q3 = q2;
            q2 = q1;
            q1 = cv;
        }
    }
    enc.finish()
}

fn dispatch_encode(
    lens: &[u32],
    binned: &[u8],
    seq: Option<&[u8]>,
    dense: &[u8; 256],
    qmin: u8,
    k: usize,
) -> Vec<u8> {
    by_size_class!(k, encode_payload, lens, binned, seq, dense, qmin, k)
}

fn decode_payload<const NM: usize>(
    lens: &[u32],
    syms: &[u8],
    seq: Option<&[u8]>,
    qmin: u8,
    k: usize,
    dec: &mut Decoder<'_>,
    quals: &mut Vec<u8>,
) -> Result<()> {
    let mut models = vec![SimpleModel::<NM>::with_active(k); N_CTX];
    let seq_mode = seq.is_some();
    let mut srest: &[u8] = seq.unwrap_or(&[]);
    for &l in lens {
        // `MODE_SEQ` reconstructs the same base/next/run features the encoder
        // used, from the already-decoded sequence handed in by the caller.
        let sread: &[u8] = if seq_mode {
            if srest.len() < l as usize {
                return Err(Error::Malformed("sequence shorter than quality lengths"));
            }
            let (sr, st) = srest.split_at(l as usize);
            srest = st;
            sr
        } else {
            &[]
        };
        let (mut q1, mut q2, mut q3) = (0u8, 0u8, 0u8);
        let mut delta = 0u8;
        let mut prev_base = u8::MAX;
        let mut run = 0usize;
        for pos in 0..l as usize {
            let c = if seq_mode {
                let base = sread[pos];
                let next = sread.get(pos + 1).copied().unwrap_or(u8::MAX);
                run = if base == prev_base { run + 1 } else { 1 };
                prev_base = base;
                context_lr(q1, q2, base_code(base), base_code(next), run)
            } else {
                context(q1, q2, q3, delta, pos)
            };
            debug_assert!(c < N_CTX);
            // SAFETY: both contexts pack into ≤18 bits, so `c < N_CTX == models.len()`.
            let dv = unsafe { models.get_unchecked_mut(c) }.decode(dec);
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
    Ok(())
}

fn dispatch_decode(
    lens: &[u32],
    syms: &[u8],
    seq: Option<&[u8]>,
    qmin: u8,
    k: usize,
    dec: &mut Decoder<'_>,
    quals: &mut Vec<u8>,
) -> Result<()> {
    by_size_class!(k, decode_payload, lens, syms, seq, qmin, k, dec, quals)
}

/// Encode per-read quality strings (sequence-blind).
///
/// `lens` gives each read's quality length; `quals` is their concatenation.
/// `binning` optionally quantizes qualities before modeling (lossy). Always
/// selects the position context (`MODE_POS`); see [`encode_seq`] to let long
/// reads condition quality on their bases.
pub fn encode(lens: &[u32], quals: &[u8], binning: QualityBinning) -> Result<Vec<u8>> {
    encode_seq(lens, quals, &[], binning)
}

/// Encode per-read quality strings, optionally conditioning long reads on their
/// bases.
///
/// `seq` is the reads' concatenated bases in the same order and per-read lengths
/// as `quals`. When it is present, non-empty, and the mean read length exceeds
/// [`SEQ_MODE_MIN_MEAN_LEN`], the stream is coded in `MODE_SEQ` (base + next
/// base + homopolymer run-length context); otherwise it falls back to the
/// sequence-blind `MODE_POS` and `seq` is ignored. Pass `&[]` for `seq` to
/// force `MODE_POS` — that is exactly what [`encode`] does. A `MODE_SEQ` stream
/// requires the decoded sequence at [`decode_seq`] time.
pub fn encode_seq(
    lens: &[u32],
    quals: &[u8],
    seq: &[u8],
    binning: QualityBinning,
) -> Result<Vec<u8>> {
    let total: usize = lens.iter().map(|&l| l as usize).sum();
    if total != quals.len() {
        return Err(Error::LengthMismatch {
            lens: total,
            quals: quals.len(),
        });
    }

    // Sequence context needs the bases to run in exact lockstep with the
    // qualities, and only pays off on long reads. If either fails, code
    // sequence-blind — never a correctness risk, just the old behaviour.
    let seq_mode = !seq.is_empty()
        && seq.len() == total
        && !lens.is_empty()
        && total / lens.len() >= SEQ_MODE_MIN_MEAN_LEN;

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
    //
    // Long reads take the binary-decomposition mixing coder (`binmix`,
    // `MODE_SEQ_BINMIX`); short reads keep the position context (`MODE_POS`). The
    // retired single-context `MODE_SEQ` stays decodable but is no longer emitted.
    let (payload, mode) = if seq_mode {
        // Long reads take the binary-decomposition mixing coder ([`binmix`]): it
        // beats the single-context model on ratio and runs faster than it.
        (
            binmix::encode(lens, &binned, seq, &dense, qmin, k),
            MODE_SEQ_BINMIX,
        )
    } else {
        (
            dispatch_encode(lens, &binned, None, &dense, qmin, k),
            MODE_POS,
        )
    };

    let mut out = Vec::with_capacity(16 + k + lens.len() + payload.len());
    out.push(FORMAT_VERSION);
    out.push(binning.tag());
    out.push(mode);
    out.push(k as u8);
    out.extend_from_slice(&syms);
    write_lens(&mut out, lens);
    out.extend_from_slice(&payload);
    Ok(out)
}

/// Whether a stream was coded in `MODE_SEQ` and so needs the block's decoded
/// sequence at [`decode_seq`] time. Lets the container peek the header cheaply and
/// decide whether to serialize seq → qual or keep decoding them in parallel,
/// without paying that serialization on the short-read common case. A truncated or
/// foreign stream reads as `false`; [`decode`]/[`decode_seq`] then reject it.
pub fn needs_sequence(src: &[u8]) -> bool {
    // Header layout: version(0), binning tag(1), mode(2), ...
    src.first() == Some(&FORMAT_VERSION)
        && matches!(src.get(2), Some(&MODE_SEQ) | Some(&MODE_SEQ_BINMIX))
}

/// Decode a sequence-blind stream produced by [`encode`], returning
/// `(lengths, qualities)`. A `MODE_SEQ` stream (from [`encode_seq`] on long
/// reads) is rejected here — use [`decode_seq`] with its decoded sequence.
pub fn decode(src: &[u8]) -> Result<(Vec<u32>, Vec<u8>)> {
    decode_seq(src, &[])
}

/// Decode a stream produced by [`encode_seq`], returning `(lengths, qualities)`.
/// In lossy modes the qualities are the binned values, not the originals.
///
/// For a `MODE_SEQ` stream, `seq` must be the block's decoded sequence (same
/// order and per-read lengths as the qualities) — see [`needs_sequence`]. For a
/// `MODE_POS` stream `seq` is ignored and may be `&[]`.
pub fn decode_seq(src: &[u8], seq: &[u8]) -> Result<(Vec<u32>, Vec<u8>)> {
    let mut r = ByteReader::new(src);
    if r.u8()? != FORMAT_VERSION {
        return Err(Error::Malformed("unsupported version"));
    }
    let _binning = QualityBinning::from_tag(r.u8()?)?;
    let mode = r.u8()?;
    // Both sequence modes need the decoded bases; `MODE_POS` does not.
    let seq_mode = match mode {
        MODE_POS => false,
        MODE_SEQ | MODE_SEQ_BINMIX => true,
        _ => return Err(Error::Malformed("unknown quality context mode")),
    };
    let k = r.u8()? as usize;
    if k == 0 || k > QMAX {
        return Err(Error::AlphabetTooLarge(k));
    }
    let syms = r.take(k)?.to_vec();
    let qmin = syms[0];
    let lens = read_lens(&mut r)?;

    // Checked sum: a malformed stream can declare lengths whose total wraps
    // `usize`, which would under-allocate `quals` and then over-push.
    let total = lens
        .iter()
        .try_fold(0usize, |acc, &l| acc.checked_add(l as usize))
        .ok_or(Error::Malformed("total length overflows usize"))?;
    let payload = r.rest();
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
    // A `MODE_SEQ` stream conditions each quality on its base, so the caller must
    // hand back exactly the decoded sequence — same total length as the qualities.
    // Reject a mismatch up front rather than panicking on a short slice mid-decode.
    if seq_mode && seq.len() != total {
        return Err(Error::Malformed(
            "sequence context requires the decoded sequence",
        ));
    }
    let mut quals = Vec::new();
    quals
        .try_reserve(total)
        .map_err(|_| Error::Malformed("declared total length too large to allocate"))?;
    match mode {
        MODE_SEQ_BINMIX => binmix::decode(&lens, payload, seq, &syms, qmin, k, &mut quals)?,
        _ => {
            let mut dec = Decoder::new(payload);
            dispatch_decode(
                &lens,
                &syms,
                seq_mode.then_some(seq),
                qmin,
                k,
                &mut dec,
                &mut quals,
            )?;
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

    /// The size-class dispatch is only legitimate if every `N >= k` codes to the
    /// SAME bytes — otherwise picking a class by `k` silently rewrites archives.
    ///
    /// That rests on a property of `SimpleModel`, which lives in a crate BELOW
    /// this one: phantom slots (`freq[k..]`) must stay at frequency 0 forever, so
    /// they never enter `tot` and never shift a code. They do, because the
    /// rescale is `(f + 1) >> 1` and `(0 + 1) >> 1 == 0`. If that ever becomes
    /// `max(1, f / 2)` — which its own comment, "halve, keep >= 1", already
    /// describes and which is the obvious "fix" for a coder that must not emit
    /// zero-probability symbols — every phantom goes live, `tot` diverges by
    /// class, and this dispatch changes the output of every archive it touches
    /// while every round-trip test still passes.
    ///
    /// So pin the property here rather than trusting a comment in another crate.
    #[test]
    fn every_size_class_codes_identically() {
        // A 3-symbol alphabet, so `k = 3` and the dispatch picks `N = 4`, while
        // the classes above it are all phantom-heavy.
        let mut quals = Vec::new();
        let mut x: u32 = 12345;
        for _ in 0..20_000 {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            quals.push(b'#' + ((x >> 16) % 3) as u8);
        }
        let lens = vec![100u32; quals.len() / 100];
        let (syms, dense) = dense_alphabet(&quals).expect("alphabet");
        let (qmin, k) = (syms[0], syms.len());
        assert_eq!(k, 3, "fixture must land in the smallest size class");

        let want = encode_payload::<QMAX>(&lens, &quals, None, &dense, qmin, k);
        for got in [
            encode_payload::<4>(&lens, &quals, None, &dense, qmin, k),
            encode_payload::<8>(&lens, &quals, None, &dense, qmin, k),
            encode_payload::<16>(&lens, &quals, None, &dense, qmin, k),
            encode_payload::<32>(&lens, &quals, None, &dense, qmin, k),
            encode_payload::<64>(&lens, &quals, None, &dense, qmin, k),
        ] {
            assert_eq!(
                got, want,
                "a size class changed the coded bytes: the model's phantom slots \
                 are no longer inert, and the dispatch is now a format change"
            );
        }
    }
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

    /// Deterministic pseudo-random long-read fixture: `n` reads of `len` bases,
    /// each quality correlated with its base so `MODE_SEQ` has real signal to fit.
    fn longread_fixture(n: usize, len: usize) -> (Vec<u32>, Vec<u8>, Vec<u8>) {
        let bases = b"ACGT";
        let mut x: u32 = 0x1234_5678;
        let mut rng = move || {
            x ^= x << 13;
            x ^= x >> 17;
            x ^= x << 5;
            x
        };
        let (mut seq, mut qual) = (Vec::new(), Vec::new());
        for _ in 0..n {
            for _ in 0..len {
                let base = bases[(rng() % 4) as usize];
                seq.push(base);
                // Quality leans on base identity plus a little noise — exactly the
                // structure the sequence context is meant to exploit.
                let bias = match base {
                    b'A' => 40,
                    b'C' => 35,
                    b'G' => 30,
                    _ => 25,
                };
                let noise = (rng() % 8) as u8;
                qual.push(33 + bias + noise);
            }
        }
        (vec![len as u32; n], seq, qual)
    }

    #[test]
    fn roundtrip_seq_context_long_reads() {
        let (lens, seq, quals) = longread_fixture(40, 2000);
        // Long reads: encode_seq must pick MODE_SEQ and round-trip through decode_seq.
        let enc = encode_seq(&lens, &quals, &seq, QualityBinning::Lossless).expect("encode");
        assert!(needs_sequence(&enc), "long reads must select MODE_SEQ");
        let (out_lens, out_quals) = decode_seq(&enc, &seq).expect("decode");
        assert_eq!(out_lens, lens);
        assert_eq!(out_quals, quals);
    }

    #[test]
    fn seq_context_beats_blind_on_correlated_quality() {
        // The whole point: when quality tracks the base, conditioning on sequence
        // must code smaller than the sequence-blind position context.
        let (lens, seq, quals) = longread_fixture(40, 2000);
        let with_seq = encode_seq(&lens, &quals, &seq, QualityBinning::Lossless).expect("seq");
        let blind = encode(&lens, &quals, QualityBinning::Lossless).expect("blind");
        assert!(
            with_seq.len() < blind.len(),
            "sequence context ({} B) should beat blind ({} B) on base-correlated quality",
            with_seq.len(),
            blind.len()
        );
    }

    #[test]
    fn seq_context_decode_requires_sequence() {
        // A MODE_SEQ stream cannot be decoded blind: the plain `decode` (empty
        // sequence) must reject it cleanly rather than produce wrong output.
        let (lens, seq, quals) = longread_fixture(20, 1000);
        let enc = encode_seq(&lens, &quals, &seq, QualityBinning::Lossless).expect("encode");
        assert!(matches!(decode(&enc), Err(Error::Malformed(_))));
        // A wrong-length sequence is rejected too.
        assert!(matches!(
            decode_seq(&enc, &seq[..seq.len() - 1]),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn short_reads_stay_sequence_blind() {
        // Below the mean-length gate, encode_seq must fall back to MODE_POS and
        // produce byte-identical output to the sequence-blind encode.
        let lens = vec![100u32; 50];
        let mut quals = Vec::new();
        let mut seq = Vec::new();
        let mut x: u32 = 99;
        for _ in 0..5000 {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            quals.push(b'#' + ((x >> 16) % 4) as u8);
            seq.push(b"ACGT"[((x >> 20) % 4) as usize]);
        }
        let with_seq = encode_seq(&lens, &quals, &seq, QualityBinning::Lossless).expect("seq");
        let blind = encode(&lens, &quals, QualityBinning::Lossless).expect("blind");
        assert!(!needs_sequence(&with_seq), "short reads must stay MODE_POS");
        assert_eq!(with_seq, blind, "short-read output must be byte-identical");
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
        let mut src = vec![FORMAT_VERSION, 0, MODE_POS, 1, b'I']; // version, lossless, mode, k=1, symbol
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
        let mut buf = vec![FORMAT_VERSION, 0, MODE_POS, 1, 0]; // version, binning, mode, k=1, syms=[0]
        push_varint(&mut buf, u64::MAX >> 8); // n: absurd length count
        buf.push(1); // fixed = true -> resize(n, f) path
        push_varint(&mut buf, 100); // f
        assert!(matches!(decode(&buf), Err(Error::Malformed(_))));
    }

    #[test]
    fn rejects_huge_total_length() {
        let mut buf = vec![FORMAT_VERSION, 0, MODE_POS, 1, 0]; // version, binning, mode, k=1, syms=[0]
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
