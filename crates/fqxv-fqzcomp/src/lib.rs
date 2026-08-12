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
//! - **Sequence context** (long reads, via [`encode_seq`]): each quality is coded
//!   as `ceil(log2 k)` binary bit-tree decisions, every bit predicted by several
//!   context tiers (coarse/mid/rich) whose probabilities are mixed in the logit
//!   domain — the binary-decomposition logistic mixer in `binmix`. The contexts
//!   are built from the neighbouring qualities and the sequence (current base,
//!   next base, homopolymer run-length), which is where HiFi/ONT quality actually
//!   lives. This drops the position/delta features (useless on long reads) and
//!   requires the decoded sequence at decode time (see
//!   [`decode_seq`]/[`needs_sequence`]).
//!
//! Long-read streams can additionally be **chunked** for within-block parallel
//! decode and encode ([`encode_seq_chunked`] with a [`QualityChunking`]): a
//! serially-coded warmup prefix trains the model, and the remaining reads
//! split into whole-read chunks coded from snapshots of the warmed state — see
//! `MODE_SEQ_BINMIX_CHUNKED` for the wire layout and
//! `docs/design/parallel-decode.md` for the measured costs.
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

use fqxv_bytes::{ReaderError, read_lens, write_lens, write_varint};
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

/// Warmup policy for chunked quality coding ([`QualityChunking`]): how many
/// leading bases the serially-coded warmup prefix covers before the model is
/// snapshotted and the remaining reads fan into parallel chunks.
///
/// The warmup is rounded up to whole reads, and the exact boundary is
/// transmitted (`warm_bases`), so the decoder never re-derives the policy —
/// only the encoder resolves it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ChunkWarmup {
    /// A fixed number of leading bases (ONT's measured operating point is
    /// 8 MiB: chunk cost becomes nearly K-invariant, +0.36% of archive at
    /// `K = 8`).
    Bases(usize),
    /// `total / K` bases — the warmup is one K-th of the block, like every
    /// other segment (HiFi's measured operating point: +0.19% of archive at
    /// `K = 8`; CCS quality drifts, so a proportional warmup tracks it).
    TotalOverSegments,
}

/// Chunk-model variant for chunked quality coding ([`QualityChunking`]) — how
/// each parallel chunk relates to the model snapshot taken after the warmup.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum ChunkVariant {
    /// Each chunk starts from a **clone** of the snapshot and keeps adapting.
    /// The default and the recommended variant on both long-read platforms;
    /// decode memory scales with concurrently-active chunks (each clone is the
    /// full model, up to ~336 MiB at long-read alphabet depth — the decoder
    /// runs chunk-at-a-time tasks so the rayon pool bounds active clones).
    #[default]
    WarmClone,
    /// Every chunk codes against the single **read-only** snapshot: `O(1)`
    /// model memory at any fan-out, at the price of no adaptation after the
    /// warmup — measured near-free on ONT (+0.61%, K-invariant) and
    /// disqualified on HiFi (+5.5–9.6%), so callers should not select it for
    /// PacBio data.
    WarmFrozen,
}

/// Chunked-quality coding parameters for [`encode_seq_chunked`]
/// (`MODE_SEQ_BINMIX_CHUNKED`): split a long-read block's quality stream into
/// a warmup prefix plus `segments - 1` chunks that encode and decode in
/// parallel.
///
/// Chunking **costs** bytes (~0.2–0.4% of archive at `segments = 8`), so it is
/// opt-in per call rather than a keep-smaller trial the serial mode would
/// always win; the container gates it by platform. Everything here is a pure
/// function of the input and these parameters — never of the thread count — so
/// the emitted stream is byte-identical at any parallelism.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct QualityChunking {
    /// Total segment count `K` (the warmup plus `K - 1` parallel chunks).
    /// Must be at least 2; the measured knee on both long-read platforms is 8.
    pub segments: usize,
    /// Warmup policy (see [`ChunkWarmup`]).
    pub warmup: ChunkWarmup,
    /// Chunk-model variant (see [`ChunkVariant`]).
    pub variant: ChunkVariant,
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
    /// The stream's quality context mode is newer than this reader knows — an
    /// archive written by a future fqxv (a new mode byte is the blessed
    /// evolution path for alternative encodings of this stream). Upgrading fqxv
    /// is the fix, so the message says so instead of reporting corruption.
    #[error("unsupported quality context mode {0}; this archive needs a newer fqxv")]
    UnsupportedMode(u8),
    /// A chunked quality stream's variant tag is one this reader does not
    /// implement (tag `0`, reset, is reserved; a future fqxv may define more).
    /// Like [`Error::UnsupportedMode`], this is an upgrade message, not
    /// corruption: the variant byte is the chunked layout's own extension
    /// point.
    #[error("unsupported quality chunk variant {0}; this archive needs a newer fqxv")]
    UnsupportedChunkVariant(u8),
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
/// Long-read binary-mixing quality model with a **per-block context quantizer**
/// ([`binmix::QCtx`]): identical coder to [`MODE_SEQ_BINMIX`], but the three
/// recent-quality context fields are quantized through a small table built from
/// this block's quality histogram (transmitted in the header after `syms`) rather
/// than the fixed `>>1`/`>>3`/`>>4` shifts. This spends full context resolution
/// where quality actually varies — the fqzcomp `qtab` / CoLoRd platform-quantizer
/// idea, which helps HiFi/Revio whose distinct high-quality values differ by one
/// Phred step (the flat `q1>>1` merges adjacent values, this keeps them apart).
///
/// Chosen only when it codes a block **strictly smaller** than `MODE_SEQ_BINMIX`:
/// [`encode_seq`] trials both and keeps the smaller, so a block can only match or
/// shrink. Needs the decoded bases exactly like `MODE_SEQ_BINMIX`.
const MODE_SEQ_BINMIX_Q: u8 = 4;
/// **Chunked** long-read binary-mixing quality model: the same coder as
/// [`MODE_SEQ_BINMIX`], but the block's reads are split into a serially-coded
/// **warmup** prefix plus `K - 1` whole-read **chunks** that encode and decode
/// in parallel, each from a snapshot of the model taken at the end of the
/// warmup (docs/design/parallel-decode.md; the format change #277 measured).
///
/// Wire layout — the chunk table sits where [`MODE_SEQ_BINMIX_Q`]'s qtable
/// sits, after `syms` and before the length array:
///
/// ```text
/// version(=2) | binning | mode=5 | k | syms[k]
/// [1] variant     (1 warm-clone, 2 warm-frozen; 0 reserved for reset)
/// [v] K           (total segment count, warmup included; 2..=MAX_CHUNK_SEGMENTS)
/// [v] warm_bases  (the warmup prefix's exact base count — a prefix sum of lens)
/// [v] warm_len    (the warmup segment's coded byte length)
/// [v] chunk_len × (K-1)   (each parallel chunk's coded byte length)
/// lens            (the existing read-length array)
/// payload         (segments concatenated: warmup, chunk 1, .., chunk K-1)
/// ```
///
/// Chunk boundaries are **never transmitted**: both sides recompute them as a
/// pure function of `(lens, K)` (`chunk_boundaries`, whole-read,
/// equal-cumulative-bases) over the reads after the warmup, and the warmup
/// read count is recovered from `warm_bases` (`warmup_reads`, self-checked to
/// land exactly on a read boundary). Segment byte lengths ARE transmitted so
/// the decoder can slice the payload and hand each chunk to a worker without
/// parsing it. Under **warm-clone** each chunk starts from a clone of the
/// snapshot and keeps adapting; under **warm-frozen** every chunk codes
/// against the one read-only snapshot (O(1) model memory at any fan-out; near
/// free on ONT, disqualified on HiFi — the variant byte is why the header
/// carries it). Chunking costs a measured ~0.2–0.4% of archive at `K = 8`
/// (the per-chunk range-coder flush plus the clone's lost adaptation), so the
/// encoder emits it only where the caller asks ([`QualityChunking`], platform
/// gated by the container) — never through a keep-smaller trial, which the
/// serial mode would always win. Needs the decoded bases like every binmix
/// mode.
const MODE_SEQ_BINMIX_CHUNKED: u8 = 5;
/// Chunked long-read binary-mixing quality model **with the per-block context
/// quantizer**: [`MODE_SEQ_BINMIX_CHUNKED`]'s layout with
/// [`MODE_SEQ_BINMIX_Q`]'s transmitted qtable (the `2 * k` bytes directly
/// after `syms`, before the chunk table). The chunked twin of the mode-3/4
/// pair: mode is to quantizer as 5 is to 6, so the qtable's presence stays a
/// pure function of the mode byte rather than a second flag mechanism.
/// [`encode_seq_chunked`] trials 5 and 6 under the same segmentation and keeps
/// the smaller, exactly as 3/4 are trialled serially.
const MODE_SEQ_BINMIX_CHUNKED_Q: u8 = 6;

/// Chunk-variant tag: warm-clone — each parallel chunk starts from a clone of
/// the warmed snapshot and keeps adapting. The shipping default on both
/// long-read platforms (ONT +0.36%, HiFi +0.19% of archive at `K = 8`).
const CHUNK_VARIANT_WARM_CLONE: u8 = 1;
/// Chunk-variant tag: warm-frozen — every chunk codes against the single
/// read-only snapshot (`O(1)` model memory at any fan-out; ONT-safe only).
const CHUNK_VARIANT_WARM_FROZEN: u8 = 2;
/// Decode-side ceiling on the transmitted segment count `K`. The encoder emits
/// single-digit K (8 at the measured knee); this bound only stops a hostile
/// chunk table from allocating per-segment state for millions of declared
/// segments before the payload-length cross-check would catch it.
const MAX_CHUNK_SEGMENTS: usize = 4096;
/// Minimum total quality bases for the encoder to emit a chunked stream.
/// Below this the block decodes fast serially and the per-chunk overhead is
/// proportionally larger, so chunking buys nothing; the trailing block of a
/// long-read file falls back to the serial mode. Decode accepts any size —
/// this gates the writer only.
const MIN_CHUNKED_BASES: usize = 16 << 20;

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

/// Code the reads `lens`/`binned` into `enc` under the per-context `models`.
///
/// The chunking twin of `binmix`'s `encode_reads`: every per-read feature (the
/// recent qualities, delta counter, position, and — in `MODE_SEQ` — the base
/// window) resets at read boundaries, so the only cross-read state is the model
/// table plus the range coder. `ADAPT = true` is the shipped path (byte-identical
/// to the pre-refactor coder; the const flag monomorphizes away); `ADAPT = false`
/// reads the models without updating them ([`SimpleModel::encode_frozen`]), the
/// frozen-snapshot variant the chunked-decode diagnostics measure.
fn encode_payload_into<const NM: usize, const ADAPT: bool>(
    models: &mut [SimpleModel<NM>],
    enc: &mut Encoder,
    lens: &[u32],
    binned: &[u8],
    seq: Option<&[u8]>,
    dense: &[u8; 256],
    qmin: u8,
) {
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
            let model = unsafe { models.get_unchecked_mut(c) };
            if ADAPT {
                model.encode(enc, dv as usize);
            } else {
                model.encode_frozen(enc, dv as usize);
            }
            if pos > 0 && cv != q1 {
                delta = (delta + 1).min(DELTA_MAX);
            }
            q3 = q2;
            q2 = q1;
            q1 = cv;
        }
    }
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
    encode_payload_into::<NM, true>(&mut models, &mut enc, lens, binned, seq, dense, qmin);
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
    // Sequence-blind (`seq = &[]`) ⇒ short-read `MODE_POS`, which never reaches the
    // quantizer trial, so the flag is moot here.
    encode_seq(lens, quals, &[], binning, false)
}

/// Encode per-read quality strings, optionally conditioning long reads on their
/// bases.
///
/// `seq` is the reads' concatenated bases in the same order and per-read lengths
/// as `quals`. When it is present, non-empty, and the mean read length exceeds
/// [`SEQ_MODE_MIN_MEAN_LEN`], the stream is coded in `MODE_SEQ_BINMIX` — the
/// binary-decomposition logistic mixer (`binmix`), which conditions on the
/// neighbouring qualities and the sequence (base, next base, homopolymer run);
/// otherwise it falls back to the sequence-blind `MODE_POS` and `seq` is ignored.
/// Pass `&[]` for `seq` to force `MODE_POS` — that is exactly what [`encode`]
/// does. Any sequence-context stream requires the decoded sequence at
/// [`decode_seq`] time. (The earlier single-context `MODE_SEQ` is still decodable
/// but is no longer emitted.)
pub fn encode_seq(
    lens: &[u32],
    quals: &[u8],
    seq: &[u8],
    binning: QualityBinning,
    // Whether to trial the per-block context quantizer (`MODE_SEQ_BINMIX_Q`). The
    // quantizer only wins on skewed quality (PacBio HiFi/Revio); on the flatter
    // Nanopore distribution it never does, so the caller passes `false` there to
    // skip the histogram build + bounded probe entirely (0 bytes saved, pure cost —
    // measured ~+5% ONT compress). The kept stream is byte-identical to the pre-#235
    // `MODE_SEQ_BINMIX` baseline either way, so this only changes speed.
    try_quantizer: bool,
) -> Result<Vec<u8>> {
    encode_seq_chunked(lens, quals, seq, binning, try_quantizer, None)
}

/// [`encode_seq`] with opt-in **chunked** coding of a long-read quality stream
/// (`MODE_SEQ_BINMIX_CHUNKED`): a serially-coded warmup prefix trains the
/// model, then the remaining reads split into `segments - 1` whole-read chunks
/// coded **in parallel**, each from a snapshot of the warmed model — and a
/// decoder fans the same chunks out, un-serializing the block's dominant
/// stream (measured 92–98% of long-read decode time) at a measured ~0.2–0.4%
/// size cost at `segments = 8`.
///
/// `chunking` is honored only where it makes sense: when the stream selects a
/// long-read (sequence-conditioned) mode, the block has at least
/// `MIN_CHUNKED_BASES` quality bytes and `segments` whole reads, and the
/// warmup leaves reads to fan. Otherwise — short reads, tiny trailing blocks,
/// `None` — the output is byte-identical to [`encode_seq`]. The emitted bytes
/// are a pure function of the input and `chunking` (chunk boundaries derive
/// from the read lengths, the parallel encode concatenates by chunk index), so
/// archives stay byte-identical at any thread count.
pub fn encode_seq_chunked(
    lens: &[u32],
    quals: &[u8],
    seq: &[u8],
    binning: QualityBinning,
    try_quantizer: bool,
    chunking: Option<QualityChunking>,
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
    // Short reads keep the position context (`MODE_POS`). Long reads take the
    // binary-decomposition mixing coder (`binmix`) and, on top of that, trial two
    // context quantizers per block and keep the smaller (never-worse by
    // construction — the same rule the sequence path uses):
    //
    // - `MODE_SEQ_BINMIX`   — the fixed `>>1`/`>>3`/`>>4` context shifts (baseline).
    // - `MODE_SEQ_BINMIX_Q` — a per-block quantizer built from this block's quality
    //   histogram (fqzcomp `qtab`). Its small table is transmitted in the header.
    //
    // The retired single-context `MODE_SEQ` stays decodable but is no longer emitted.
    if seq_mode {
        // Chunked coding is opt-in and platform-gated by the caller, never a
        // keep-smaller trial — a chunked stream always costs a little over the
        // serial one (per-chunk flush + the snapshot's lost adaptation), so a
        // size trial would always keep serial and the mode would never ship.
        // Ineligible blocks (too small, too few reads) fall back to the serial
        // modes, byte-identical to a `None` request.
        if let Some(c) = chunking.filter(|c| chunking_eligible(lens, total, c)) {
            return encode_binmix_chunked(
                lens,
                &binned,
                seq,
                &dense,
                qmin,
                k,
                &syms,
                binning,
                try_quantizer,
                c,
            );
        }
        let flat = binmix::QCtx::flat();
        let base_payload = binmix::encode(lens, &binned, seq, &dense, qmin, k, &flat);
        let baseline = assemble_quality_stream(
            binning,
            MODE_SEQ_BINMIX,
            k,
            &syms,
            &[],
            &[],
            lens,
            &base_payload,
        );

        // On regimes where the quantizer never wins (Nanopore), skip the histogram
        // build and the probe outright — the baseline is what would be kept anyway.
        if !try_quantizer {
            if diag_qchunk_enabled() {
                qchunk_probe_binmix(
                    lens,
                    &binned,
                    seq,
                    &dense,
                    qmin,
                    k,
                    &flat,
                    MODE_SEQ_BINMIX,
                    base_payload.len(),
                );
            }
            return Ok(baseline);
        }

        // The per-block context quantizer, built from the full histogram.
        let (qc, qtable) = build_quant_ctx(&binned, &syms, qmin);

        // Bound the trial cost. Coding a second full quality stream to keep the
        // smaller doubles the (dominant) long-read quality-encode work, and the
        // quantizer only wins on **skewed** quality — HiFi/Revio, where the top of
        // the Phred scale is packed — never on the flatter Nanopore distribution
        // (measured: 0 bytes saved, so the second full encode is pure waste there).
        // So probe a bounded prefix under both quantizers first and only pay the
        // second full encode when the quantizer wins the probe. Skipping keeps the
        // baseline, so this stays never-worse *by construction*; the probe can only
        // ever forgo a win it failed to predict, never enlarge a block.
        if !quant_wins_probe(lens, &binned, seq, &dense, qmin, k, &flat, &qc) {
            if diag_qchunk_enabled() {
                qchunk_probe_binmix(
                    lens,
                    &binned,
                    seq,
                    &dense,
                    qmin,
                    k,
                    &flat,
                    MODE_SEQ_BINMIX,
                    base_payload.len(),
                );
            }
            return Ok(baseline);
        }

        let q_payload = binmix::encode(lens, &binned, seq, &dense, qmin, k, &qc);
        let quantized = assemble_quality_stream(
            binning,
            MODE_SEQ_BINMIX_Q,
            k,
            &syms,
            &qtable,
            &[],
            lens,
            &q_payload,
        );

        // Keep the smaller; a tie keeps the baseline (no table, no swap for no
        // gain). Fixed candidate order + this tie-break make the choice
        // thread-independent and byte-deterministic.
        //
        // The FQXV_DIAG_QCHUNK sweep runs at this winner-decision point — under
        // the winning quantizer, against the winning payload — so its deltas
        // are relative to the bytes an archive would actually carry. Reporting
        // only: the returned stream is untouched either way.
        Ok(if quantized.len() < baseline.len() {
            if diag_qchunk_enabled() {
                qchunk_probe_binmix(
                    lens,
                    &binned,
                    seq,
                    &dense,
                    qmin,
                    k,
                    &qc,
                    MODE_SEQ_BINMIX_Q,
                    q_payload.len(),
                );
            }
            quantized
        } else {
            if diag_qchunk_enabled() {
                qchunk_probe_binmix(
                    lens,
                    &binned,
                    seq,
                    &dense,
                    qmin,
                    k,
                    &flat,
                    MODE_SEQ_BINMIX,
                    base_payload.len(),
                );
            }
            baseline
        })
    } else {
        let payload = dispatch_encode(lens, &binned, None, &dense, qmin, k);
        if diag_qchunk_enabled() {
            dispatch_qchunk_probe_pos(lens, &binned, &dense, qmin, k, payload.len());
        }
        Ok(assemble_quality_stream(
            binning,
            MODE_POS,
            k,
            &syms,
            &[],
            &[],
            lens,
            &payload,
        ))
    }
}

/// Assemble a full quality stream: `version | binning | mode | k | syms[k] |
/// qtable | chunk_table | lens | payload`. `qtable` is empty for every mode
/// except [`MODE_SEQ_BINMIX_Q`] / [`MODE_SEQ_BINMIX_CHUNKED_Q`], whose
/// per-block context quantizer is transmitted here (after `syms`, before the
/// length array) so the decoder can rebuild it; `chunk_table` is empty for
/// every mode except the two chunked ones, whose segment table
/// (see [`MODE_SEQ_BINMIX_CHUNKED`]) follows the qtable slot.
#[allow(clippy::too_many_arguments)] // one arg per header field, in wire order
fn assemble_quality_stream(
    binning: QualityBinning,
    mode: u8,
    k: usize,
    syms: &[u8],
    qtable: &[u8],
    chunk_table: &[u8],
    lens: &[u32],
    payload: &[u8],
) -> Vec<u8> {
    let mut out =
        Vec::with_capacity(16 + k + qtable.len() + chunk_table.len() + lens.len() + payload.len());
    out.push(FORMAT_VERSION);
    out.push(binning.tag());
    out.push(mode);
    out.push(k as u8);
    out.extend_from_slice(syms);
    out.extend_from_slice(qtable);
    out.extend_from_slice(chunk_table);
    write_lens(&mut out, lens);
    out.extend_from_slice(payload);
    out
}

/// Whether the encoder should honor a chunking request for this block: the
/// segment count is sane, the block carries at least [`MIN_CHUNKED_BASES`]
/// quality bytes and one whole read per segment, and the warmup prefix leaves
/// reads for the parallel chunks. A `false` falls back to the serial modes —
/// never an error, exactly the `chunking: None` output.
fn chunking_eligible(lens: &[u32], total: usize, c: &QualityChunking) -> bool {
    if c.segments < 2 || c.segments > MAX_CHUNK_SEGMENTS {
        return false;
    }
    if total < MIN_CHUNKED_BASES || lens.len() < c.segments {
        return false;
    }
    warmup_reads(lens, chunk_warm_target(total, c)) < lens.len()
}

/// The warmup prefix's target size in bases for a chunking request — the
/// encoder-side resolution of [`ChunkWarmup`]. The achieved boundary (whole
/// reads covering at least this) is what the stream transmits; the decoder
/// never re-derives the policy.
fn chunk_warm_target(total: usize, c: &QualityChunking) -> usize {
    match c.warmup {
        ChunkWarmup::Bases(b) => b.min(total),
        ChunkWarmup::TotalOverSegments => total.div_ceil(c.segments),
    }
}

/// Code a long-read block's quality as a **chunked** stream
/// ([`MODE_SEQ_BINMIX_CHUNKED`]), trialling the per-block context quantizer
/// ([`MODE_SEQ_BINMIX_CHUNKED_Q`]) under the same segmentation when
/// `try_quantizer` is set — the chunked twin of the serial mode-3/4 trial,
/// with the same probe gate, keep-smaller rule, and tie-break.
#[allow(clippy::too_many_arguments)] // mirrors the coder's input list
fn encode_binmix_chunked(
    lens: &[u32],
    binned: &[u8],
    seq: &[u8],
    dense: &[u8; 256],
    qmin: u8,
    k: usize,
    syms: &[u8],
    binning: QualityBinning,
    try_quantizer: bool,
    c: QualityChunking,
) -> Result<Vec<u8>> {
    let flat = binmix::QCtx::flat();
    let (base_table, base_payload) =
        encode_chunked_payload(lens, binned, seq, dense, qmin, k, &flat, c);
    let baseline = assemble_quality_stream(
        binning,
        MODE_SEQ_BINMIX_CHUNKED,
        k,
        syms,
        &[],
        &base_table,
        lens,
        &base_payload,
    );
    if !try_quantizer {
        return Ok(baseline);
    }
    let (qc, qtable) = build_quant_ctx(binned, syms, qmin);
    if !quant_wins_probe(lens, binned, seq, dense, qmin, k, &flat, &qc) {
        return Ok(baseline);
    }
    let (q_table, q_payload) = encode_chunked_payload(lens, binned, seq, dense, qmin, k, &qc, c);
    let quantized = assemble_quality_stream(
        binning,
        MODE_SEQ_BINMIX_CHUNKED_Q,
        k,
        syms,
        &qtable,
        &q_table,
        lens,
        &q_payload,
    );
    // Keep the smaller; a tie keeps the un-quantized baseline — the same
    // deterministic candidate order and tie-break as the serial trial.
    Ok(if quantized.len() < baseline.len() {
        quantized
    } else {
        baseline
    })
}

/// Code the chunked layout's segments under quantizer `qc`: the warmup prefix
/// serially from a fresh model, then the `segments - 1` chunks — **in
/// parallel** from per-chunk clones of the warmed snapshot (warm-clone), or
/// serially against the shared read-only snapshot (warm-frozen, whose chunks
/// are cheap enough that parallel encode isn't worth a second code path).
/// Returns `(chunk_table, payload)` for [`assemble_quality_stream`].
///
/// Determinism: the segmentation is a pure function of `(lens, c)`, each
/// chunk's bytes are a pure function of (its reads, the snapshot), and the
/// results are concatenated by chunk index — so the emitted bytes are
/// identical whether the chunks are coded serially or by any number of rayon
/// workers. Encode memory: one model clone per *running* chunk task, bounded
/// by the rayon pool's thread count, never by `segments`.
#[allow(clippy::too_many_arguments)] // mirrors the coder's input list
fn encode_chunked_payload(
    lens: &[u32],
    binned: &[u8],
    seq: &[u8],
    dense: &[u8; 256],
    qmin: u8,
    k: usize,
    qc: &binmix::QCtx,
    c: QualityChunking,
) -> (Vec<u8>, Vec<u8>) {
    use rayon::prelude::*;
    let total = binned.len();
    let wreads = warmup_reads(lens, chunk_warm_target(total, &c));
    let offs = prefix_offsets(lens);
    let wbases = offs[wreads];
    let (warm, warm_payload) = binmix::warm_mixer(
        &lens[..wreads],
        &binned[..wbases],
        &seq[..wbases],
        dense,
        qmin,
        k,
        qc,
    );
    let n_chunks = c.segments - 1;
    let bounds = chunk_boundaries(&lens[wreads..], n_chunks);
    let chunks: Vec<Vec<u8>> = match c.variant {
        ChunkVariant::WarmClone => (0..n_chunks)
            .into_par_iter()
            .map(|i| {
                let (a, b) = (wreads + bounds[i], wreads + bounds[i + 1]);
                binmix::encode_chunk::<true>(
                    &mut warm.clone(),
                    &lens[a..b],
                    &binned[offs[a]..offs[b]],
                    &seq[offs[a]..offs[b]],
                    dense,
                    qmin,
                    qc,
                )
            })
            .collect(),
        ChunkVariant::WarmFrozen => {
            // Frozen never adapts (`ADAPT = false` provably leaves the mixer
            // untouched), so one serial pass over the shared snapshot codes
            // every chunk in O(1) model memory — the variant's whole point.
            let mut mx = warm;
            (0..n_chunks)
                .map(|i| {
                    let (a, b) = (wreads + bounds[i], wreads + bounds[i + 1]);
                    binmix::encode_chunk::<false>(
                        &mut mx,
                        &lens[a..b],
                        &binned[offs[a]..offs[b]],
                        &seq[offs[a]..offs[b]],
                        dense,
                        qmin,
                        qc,
                    )
                })
                .collect()
        }
    };
    let mut table = Vec::with_capacity(1 + 10 * (n_chunks + 3));
    table.push(match c.variant {
        ChunkVariant::WarmClone => CHUNK_VARIANT_WARM_CLONE,
        ChunkVariant::WarmFrozen => CHUNK_VARIANT_WARM_FROZEN,
    });
    write_varint(&mut table, c.segments as u64);
    write_varint(&mut table, wbases as u64);
    write_varint(&mut table, warm_payload.len() as u64);
    for ch in &chunks {
        write_varint(&mut table, ch.len() as u64);
    }
    let mut payload = warm_payload;
    for ch in &chunks {
        payload.extend_from_slice(ch);
    }
    (table, payload)
}

/// Per-read prefix offsets into the concatenated quality/sequence buffers:
/// `lens.len() + 1` entries with `offs[i]` = bases before read `i`.
fn prefix_offsets(lens: &[u32]) -> Vec<usize> {
    let mut offs = Vec::with_capacity(lens.len() + 1);
    offs.push(0usize);
    for &l in lens {
        offs.push(offs.last().copied().unwrap_or(0) + l as usize);
    }
    offs
}

/// Per-block context quantizer for the long-read binary-mixing coder
/// ([`MODE_SEQ_BINMIX_Q`]).
///
/// Builds the three per-`cv` bucket tables ([`binmix::QCtx`]) that replace the
/// fixed `>>1`/`>>3`/`>>4` context shifts, plus the compact table transmitted in
/// the header so the decoder reconstructs the identical quantizer.
///
/// The `q1` field (6 bits, coarse tier) gets **full resolution**: for the small
/// alphabets long-read quality uses (Revio's ~7 levels, HiFi's clustered high-Q
/// values), each distinct value maps to its own bucket — the win the flat `q1>>1`
/// throws away by merging adjacent Phred values. Alphabets past 64 symbols fold by
/// equal population. The `q2`/`q3` fields (3/2 bits, mid/rich tiers) are
/// equal-population folds of the marginal histogram, so each coarse bucket carries
/// comparable mass (the fqzcomp `qtab` heuristic).
///
/// Returns `(quantizer, table)` where `table` is `2 * k` bytes: for each present
/// symbol (in ascending value order) `[t1, (t2 << 3) | t3]`.
fn build_quant_ctx(binned: &[u8], syms: &[u8], qmin: u8) -> (binmix::QCtx, Vec<u8>) {
    let k = syms.len();
    let mut cnt = [0u64; 256];
    for &b in binned {
        cnt[b as usize] += 1;
    }
    let total: u64 = cnt.iter().sum::<u64>().max(1);
    let mut g1 = [0u8; 256];
    let mut g2 = [0u8; 256];
    let mut g3 = [0u8; 256];
    let mut table = Vec::with_capacity(2 * k);
    let mut cum = 0u64;
    for (i, &b) in syms.iter().enumerate() {
        // `cum` is the count of all lower-value symbols (equal-population boundary).
        let t1 = if k <= 64 {
            i as u8 // full resolution: distinct value -> distinct 6-bit bucket
        } else {
            ((cum * 64 / total).min(63)) as u8
        };
        let t2 = ((cum * 8 / total).min(7)) as u8;
        let t3 = ((cum * 4 / total).min(3)) as u8;
        let cv = (b - qmin) as usize;
        g1[cv] = t1;
        g2[cv] = t2;
        g3[cv] = t3;
        table.push(t1);
        table.push((t2 << 3) | t3);
        cum += cnt[b as usize];
    }
    (binmix::QCtx::from_tables(g1, g2, g3), table)
}

/// Leading whole-read bases to code under both quantizers when deciding whether
/// the quantizer is worth a full second encode. The quantizer's structural edge on
/// skewed HiFi/Revio quality shows up in the first few MiB; 8 MiB is a small slice
/// of a 256 MiB block (so the probe costs a few percent) yet enough reads to be
/// representative. Blocks at or below this size are trialled whole (the "probe"
/// would be the entire block, so `keep_smaller` just decides directly).
const QUANT_PROBE_BASES: usize = 8 << 20;

/// Decide whether the per-block context quantizer is worth a full second encode,
/// by coding a bounded leading prefix under both the flat and quantized contexts
/// and checking whether the quantizer codes it strictly smaller.
///
/// This is a **cost gate only**: a `false` keeps the baseline stream (never-worse
/// is preserved regardless), so a mis-prediction can only forgo a win, never
/// enlarge a block. It exists so Nanopore — where the quantizer never wins — does
/// not pay for a second full quality encode.
#[allow(clippy::too_many_arguments)] // mirrors the coder's input list
fn quant_wins_probe(
    lens: &[u32],
    binned: &[u8],
    seq: &[u8],
    dense: &[u8; 256],
    qmin: u8,
    k: usize,
    flat: &binmix::QCtx,
    qc: &binmix::QCtx,
) -> bool {
    let total = binned.len();
    if total <= QUANT_PROBE_BASES {
        // Small block: the probe would be the whole block, so skip it and let
        // `keep_smaller` decide on the full encodes — cheaper than coding twice.
        return true;
    }
    // Take leading whole reads until they cover at least the probe budget.
    let mut nreads = 0usize;
    let mut acc = 0usize;
    for &l in lens {
        acc += l as usize;
        nreads += 1;
        if acc >= QUANT_PROBE_BASES {
            break;
        }
    }
    let plens = &lens[..nreads];
    let pq = &binned[..acc];
    let ps = &seq[..acc];
    let flat_probe = binmix::encode(plens, pq, ps, dense, qmin, k, flat);
    let quant_probe = binmix::encode(plens, pq, ps, dense, qmin, k, qc);
    quant_probe.len() < flat_probe.len()
}

// --- FQXV_DIAG_QCHUNK: chunked-quality cost sweep (measurement only) ---------
//
// Encoder-side diagnostics for the parallel-first quality design
// (docs/design/parallel-decode.md): when `FQXV_DIAG_QCHUNK` is set, every
// quality block additionally re-codes its (binned) qualities under a grid of
// chunk layouts — {reset, warm-clone, warm-frozen} x K in {4, 8, 16, 32} x
// warmup in {total/K, 8 MiB} — and prints one machine-parsable stderr line per
// cell. The probe is REPORTING ONLY: it reads the same inputs the winning coder
// consumed and never touches the returned stream, so the emitted archive is
// byte-identical with the flag set or unset (bench/scripts/qchunk_diag.sh
// enforces that with `cmp` on every sweep dataset; precedent: the container's
// FQXV_DIAG_SEQ). Off by default, zero cost when unset.

/// The chunk counts the sweep measures. For the warm variants `K` counts the
/// warmup segment: warmup + `K - 1` parallel chunks, so the decode-speedup
/// model — `T_qual/K` for reset, `w*T_qual + (1-w)*T_qual/(K-1)` for warm —
/// reads straight off the cells.
const DIAG_QCHUNK_KS: [usize; 4] = [4, 8, 16, 32];
/// The fixed-size warmup policy: 8 MiB of leading whole reads (the same budget
/// the quantizer probe uses — a few percent of a 256 MiB block).
const DIAG_QCHUNK_WARM_FIXED: usize = 8 << 20;

/// Whether the chunked-quality sweep is enabled (`FQXV_DIAG_QCHUNK` set).
fn diag_qchunk_enabled() -> bool {
    std::env::var_os("FQXV_DIAG_QCHUNK").is_some()
}

/// Arrival-order id for diag lines, so the ~20 cells of one block can be
/// grouped. NOT a stable block index: blocks are coded in parallel, so the
/// numbering depends on thread scheduling — it identifies lines from one call,
/// nothing more (aggregation sums over all lines and never needs block
/// identity).
fn next_diag_blk() -> u64 {
    use std::sync::atomic::{AtomicU64, Ordering};
    static DIAG_QCHUNK_BLK: AtomicU64 = AtomicU64::new(0);
    DIAG_QCHUNK_BLK.fetch_add(1, Ordering::Relaxed)
}

/// LEB128 length of `v` — the bytes a varint of it would occupy on the wire.
fn varint_len(mut v: u64) -> usize {
    let mut n = 1;
    while v >= 0x80 {
        v >>= 7;
        n += 1;
    }
    n
}

/// Contiguous whole-read chunk boundaries with ~equal cumulative bases: returns
/// `n_chunks + 1` read-index boundaries `[0, ..., lens.len()]`, boundary `c`
/// being the first read at or past `c/n_chunks` of the total bases. A pure
/// function of `(lens, n_chunks)` — never thread- or content-derived — so both
/// sides of a chunked format would recompute identical boundaries from the
/// length array they already share (the determinism template
/// `GlobalReference::encode_blocked` uses in fqxv-lroverlap).
fn chunk_boundaries(lens: &[u32], n_chunks: usize) -> Vec<usize> {
    let total: u64 = lens.iter().map(|&l| u64::from(l)).sum();
    let n = n_chunks.max(1);
    let mut out = Vec::with_capacity(n + 1);
    out.push(0usize);
    let mut cum = 0u64;
    let mut i = 0usize;
    for c in 1..n {
        let target = total * c as u64 / n as u64;
        while i < lens.len() && cum < target {
            cum += u64::from(lens[i]);
            i += 1;
        }
        out.push(i);
    }
    out.push(lens.len());
    out
}

/// Number of leading whole reads covering at least `min_bases` (all reads when
/// the block is smaller): the warmup prefix of the warm chunk variants.
fn warmup_reads(lens: &[u32], min_bases: usize) -> usize {
    let mut acc = 0usize;
    for (i, &l) in lens.iter().enumerate() {
        acc += l as usize;
        if acc >= min_bases {
            return i + 1;
        }
    }
    lens.len()
}

/// Emit one sweep cell. `chunk_sizes` are the parallel chunks' coded bytes
/// (excluding the warmup segment). The charged header (`hdr`) is what a chunked
/// mode would transmit where the qtable sits: a variant byte, the segment
/// count, the warmup boundary (bases) and warmup segment length, and each
/// parallel chunk's byte length as varints — chunk boundaries themselves are
/// recomputed from the length array both sides already share, never
/// transmitted.
#[allow(clippy::too_many_arguments)] // a report row, one field per arg
fn qchunk_emit(
    blk: u64,
    mode: u8,
    n_reads: usize,
    bases: usize,
    k: usize,
    baseline: usize,
    variant: &str,
    n_segments: usize,
    warmup: &str,
    warm_bases: usize,
    warm_bytes: usize,
    chunk_sizes: &[usize],
) {
    let hdr = 1
        + varint_len(n_segments as u64)
        + varint_len(warm_bases as u64)
        + if warm_bytes > 0 {
            varint_len(warm_bytes as u64)
        } else {
            0
        }
        + chunk_sizes
            .iter()
            .map(|&c| varint_len(c as u64))
            .sum::<usize>();
    let chunk_total: usize = chunk_sizes.iter().sum();
    let total = warm_bytes + chunk_total + hdr;
    let delta_pct = (total as f64 - baseline as f64) * 100.0 / baseline.max(1) as f64;
    eprintln!(
        "[diag qchunk] blk={blk} mode={mode} reads={n_reads} bases={bases} k={k} \
         baseline={baseline} variant={variant} K={n_segments} warmup={warmup} \
         warm_bases={warm_bases} warm_bytes={warm_bytes} chunk_bytes={chunk_total} \
         hdr={hdr} total={total} delta_pct={delta_pct:.4}"
    );
}

/// The long-read (binmix) sweep, run for the WINNING mode/quantizer of a block
/// so `baseline` is the payload that would actually ship. Costs ~20 extra
/// full-block quality encodes; snapshots are warmed once per (policy, K) and
/// shared by the warm-clone/warm-frozen cells.
#[allow(clippy::too_many_arguments)] // mirrors the coder's input list
fn qchunk_probe_binmix(
    lens: &[u32],
    binned: &[u8],
    seq: &[u8],
    dense: &[u8; 256],
    qmin: u8,
    k: usize,
    qc: &binmix::QCtx,
    mode: u8,
    baseline: usize,
) {
    let total = binned.len();
    if lens.is_empty() || total == 0 {
        return;
    }
    let blk = next_diag_blk();
    let mut offs = Vec::with_capacity(lens.len() + 1);
    offs.push(0usize);
    for &l in lens {
        offs.push(offs.last().copied().unwrap_or(0) + l as usize);
    }

    // reset: K independent cold chunks (a decoder would run all K in parallel,
    // each owning a mutable model — the memory-hungry variant).
    for &kc in &DIAG_QCHUNK_KS {
        let bounds = chunk_boundaries(lens, kc);
        let sizes: Vec<usize> = bounds
            .windows(2)
            .map(|w| {
                binmix::chunk_bytes::<true>(
                    &mut binmix::BinMixer::new(k),
                    &lens[w[0]..w[1]],
                    &binned[offs[w[0]]..offs[w[1]]],
                    &seq[offs[w[0]]..offs[w[1]]],
                    dense,
                    qmin,
                    qc,
                )
            })
            .collect();
        qchunk_emit(
            blk,
            mode,
            lens.len(),
            total,
            k,
            baseline,
            "reset",
            kc,
            "none",
            0,
            0,
            &sizes,
        );
    }

    // Warm variants: segment 0 codes the warmup prefix from scratch (decoded
    // serially, training the snapshot); the remaining reads split into K - 1
    // chunks, each coded from the snapshot — cloned-and-adapting (warm-clone)
    // or shared read-only (warm-frozen, the O(1)-memory variant).
    let warm_cells = |policy: &str, warm_target: usize, kc: usize| {
        let wreads = warmup_reads(lens, warm_target.min(total));
        let wbases = offs[wreads];
        let (warm, warm_seg) = binmix::warm_mixer(
            &lens[..wreads],
            &binned[..wbases],
            &seq[..wbases],
            dense,
            qmin,
            k,
            qc,
        );
        let warm_bytes = warm_seg.len();
        let bounds = chunk_boundaries(&lens[wreads..], kc - 1);
        for (variant, adapt) in [("warmclone", true), ("warmfrozen", false)] {
            let mut frozen = warm.clone(); // one read-only state for every frozen chunk
            let sizes: Vec<usize> = bounds
                .windows(2)
                .map(|w| {
                    let (a, b) = (wreads + w[0], wreads + w[1]);
                    let (l, q, s) = (
                        &lens[a..b],
                        &binned[offs[a]..offs[b]],
                        &seq[offs[a]..offs[b]],
                    );
                    if adapt {
                        binmix::chunk_bytes::<true>(&mut warm.clone(), l, q, s, dense, qmin, qc)
                    } else {
                        binmix::chunk_bytes::<false>(&mut frozen, l, q, s, dense, qmin, qc)
                    }
                })
                .collect();
            qchunk_emit(
                blk,
                mode,
                lens.len(),
                total,
                k,
                baseline,
                variant,
                kc,
                policy,
                wbases,
                warm_bytes,
                &sizes,
            );
        }
    };
    for &kc in &DIAG_QCHUNK_KS {
        warm_cells("perk", total.div_ceil(kc), kc);
        warm_cells("8mib", DIAG_QCHUNK_WARM_FIXED, kc);
    }
}

/// One chunk of the `MODE_POS` sweep coded from `models` with a fresh range
/// coder (per-chunk flush included) — the `SimpleModel`-table twin of
/// `binmix::chunk_bytes`.
fn pos_chunk_bytes<const NM: usize, const ADAPT: bool>(
    models: &mut [SimpleModel<NM>],
    lens: &[u32],
    binned: &[u8],
    dense: &[u8; 256],
    qmin: u8,
) -> usize {
    let mut enc = Encoder::new();
    encode_payload_into::<NM, ADAPT>(models, &mut enc, lens, binned, None, dense, qmin);
    enc.finish().len()
}

/// The short-read (`MODE_POS`) sweep — the NovaSeq control. Same grid as
/// [`qchunk_probe_binmix`], with the per-context `SimpleModel` table (8–48 MiB)
/// as the snapshotted state.
fn qchunk_probe_pos<const NM: usize>(
    lens: &[u32],
    binned: &[u8],
    dense: &[u8; 256],
    qmin: u8,
    k: usize,
    baseline: usize,
) {
    let total = binned.len();
    if lens.is_empty() || total == 0 {
        return;
    }
    let blk = next_diag_blk();
    let mut offs = Vec::with_capacity(lens.len() + 1);
    offs.push(0usize);
    for &l in lens {
        offs.push(offs.last().copied().unwrap_or(0) + l as usize);
    }

    for &kc in &DIAG_QCHUNK_KS {
        let bounds = chunk_boundaries(lens, kc);
        let sizes: Vec<usize> = bounds
            .windows(2)
            .map(|w| {
                let mut models = vec![SimpleModel::<NM>::with_active(k); N_CTX];
                pos_chunk_bytes::<NM, true>(
                    &mut models,
                    &lens[w[0]..w[1]],
                    &binned[offs[w[0]]..offs[w[1]]],
                    dense,
                    qmin,
                )
            })
            .collect();
        qchunk_emit(
            blk,
            MODE_POS,
            lens.len(),
            total,
            k,
            baseline,
            "reset",
            kc,
            "none",
            0,
            0,
            &sizes,
        );
    }

    let warm_cells = |policy: &str, warm_target: usize, kc: usize| {
        let wreads = warmup_reads(lens, warm_target.min(total));
        let wbases = offs[wreads];
        let mut warm = vec![SimpleModel::<NM>::with_active(k); N_CTX];
        let warm_bytes = {
            let mut enc = Encoder::new();
            encode_payload_into::<NM, true>(
                &mut warm,
                &mut enc,
                &lens[..wreads],
                &binned[..wbases],
                None,
                dense,
                qmin,
            );
            enc.finish().len()
        };
        let bounds = chunk_boundaries(&lens[wreads..], kc - 1);
        for (variant, adapt) in [("warmclone", true), ("warmfrozen", false)] {
            let mut frozen = warm.clone();
            let sizes: Vec<usize> = bounds
                .windows(2)
                .map(|w| {
                    let (a, b) = (wreads + w[0], wreads + w[1]);
                    let (l, q) = (&lens[a..b], &binned[offs[a]..offs[b]]);
                    if adapt {
                        pos_chunk_bytes::<NM, true>(&mut warm.clone(), l, q, dense, qmin)
                    } else {
                        pos_chunk_bytes::<NM, false>(&mut frozen, l, q, dense, qmin)
                    }
                })
                .collect();
            qchunk_emit(
                blk,
                MODE_POS,
                lens.len(),
                total,
                k,
                baseline,
                variant,
                kc,
                policy,
                wbases,
                warm_bytes,
                &sizes,
            );
        }
    };
    for &kc in &DIAG_QCHUNK_KS {
        warm_cells("perk", total.div_ceil(kc), kc);
        warm_cells("8mib", DIAG_QCHUNK_WARM_FIXED, kc);
    }
}

/// [`qchunk_probe_pos`] behind the size-class dispatch (the probe must cost
/// chunks with the same `NM` the shipped payload used, or the trained state —
/// and so the measured deltas — would differ from reality).
fn dispatch_qchunk_probe_pos(
    lens: &[u32],
    binned: &[u8],
    dense: &[u8; 256],
    qmin: u8,
    k: usize,
    baseline: usize,
) {
    by_size_class!(k, qchunk_probe_pos, lens, binned, dense, qmin, k, baseline)
}

/// Rebuild the [`MODE_SEQ_BINMIX_Q`] context quantizer from its transmitted table
/// (`2 * k` bytes, see [`build_quant_ctx`]) and the present symbols.
fn read_quant_ctx(table: &[u8], syms: &[u8], qmin: u8) -> Result<binmix::QCtx> {
    let k = syms.len();
    if table.len() != 2 * k {
        return Err(Error::Malformed("quality quantizer table length mismatch"));
    }
    let mut g1 = [0u8; 256];
    let mut g2 = [0u8; 256];
    let mut g3 = [0u8; 256];
    for (i, &b) in syms.iter().enumerate() {
        let t1 = table[2 * i];
        let packed = table[2 * i + 1];
        let cv = (b - qmin) as usize;
        // Masks are defensive; `keys` masks each field again, so an out-of-range
        // byte can only mis-predict, never index out of bounds.
        g1[cv] = t1 & 0x3F;
        g2[cv] = (packed >> 3) & 0x7;
        g3[cv] = packed & 0x3;
    }
    Ok(binmix::QCtx::from_tables(g1, g2, g3))
}

/// Whether a stream was coded in a sequence-conditioned mode and so needs the
/// block's decoded sequence at [`decode_seq`] time. Lets the container peek the
/// header cheaply and decide whether to serialize seq → qual or keep decoding
/// them in parallel, without paying that serialization on the short-read common
/// case. A truncated or foreign stream reads as `false`;
/// [`decode`]/[`decode_seq`] then reject it.
pub fn needs_sequence(src: &[u8]) -> bool {
    // Header layout: version(0), binning tag(1), mode(2), ...
    //
    // Fail-safe on the mode byte: only `MODE_POS` is known sequence-free, so
    // any OTHER mode under the recognized version — including one this reader
    // has never heard of (a future coder, e.g. a chunked mode 5) — reports
    // `true`. The container then hands [`decode_seq`] the decoded bases and it
    // surfaces the clean [`Error::UnsupportedMode`] ("needs a newer fqxv")
    // instead of a spurious quality-only decode failure. The previous
    // enumerated match was an implicit extension point: an unknown mode read as
    // `false`, so a future stream would have been decoded WITHOUT its sequence
    // and died with a generic corruption error.
    src.first() == Some(&FORMAT_VERSION) && matches!(src.get(2), Some(&m) if m != MODE_POS)
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
    // Every sequence-conditioned mode needs the decoded bases; `MODE_POS` does not.
    let seq_mode = match mode {
        MODE_POS => false,
        MODE_SEQ
        | MODE_SEQ_BINMIX
        | MODE_SEQ_BINMIX_Q
        | MODE_SEQ_BINMIX_CHUNKED
        | MODE_SEQ_BINMIX_CHUNKED_Q => true,
        // Not corruption: the mode byte is the stream's extension point, so an
        // unknown value is a future coder this reader must refuse by name.
        m => return Err(Error::UnsupportedMode(m)),
    };
    let k = r.u8()? as usize;
    if k == 0 || k > QMAX {
        return Err(Error::AlphabetTooLarge(k));
    }
    let syms = r.take(k)?.to_vec();
    let qmin = syms[0];
    // The quantizer modes transmit the per-block context quantizer here, after the
    // symbol alphabet and before the length array (see `assemble_quality_stream`).
    let qtable: Vec<u8> = if mode == MODE_SEQ_BINMIX_Q || mode == MODE_SEQ_BINMIX_CHUNKED_Q {
        r.take(2 * k)?.to_vec()
    } else {
        Vec::new()
    };
    // The chunked modes transmit their segment table next (after the qtable
    // slot, before the length array — see `MODE_SEQ_BINMIX_CHUNKED`).
    let chunk_table = if mode == MODE_SEQ_BINMIX_CHUNKED || mode == MODE_SEQ_BINMIX_CHUNKED_Q {
        Some(read_chunk_table(&mut r)?)
    } else {
        None
    };
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
        MODE_SEQ_BINMIX => {
            let qc = binmix::QCtx::flat();
            quals.resize(total, 0);
            binmix::decode(&lens, payload, seq, &syms, qmin, k, &qc, &mut quals)?;
        }
        MODE_SEQ_BINMIX_Q => {
            let qc = read_quant_ctx(&qtable, &syms, qmin)?;
            quals.resize(total, 0);
            binmix::decode(&lens, payload, seq, &syms, qmin, k, &qc, &mut quals)?;
        }
        MODE_SEQ_BINMIX_CHUNKED | MODE_SEQ_BINMIX_CHUNKED_Q => {
            let qc = if mode == MODE_SEQ_BINMIX_CHUNKED_Q {
                read_quant_ctx(&qtable, &syms, qmin)?
            } else {
                binmix::QCtx::flat()
            };
            // Parsed above whenever the mode is chunked, so always present here.
            let t = chunk_table
                .as_ref()
                .ok_or(Error::Malformed("missing chunk table"))?;
            quals.resize(total, 0);
            decode_chunked(t, &lens, payload, seq, &syms, qmin, k, &qc, &mut quals)?;
        }
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

/// Parsed segment table of a chunked quality stream (see
/// [`MODE_SEQ_BINMIX_CHUNKED`]): the variant tag, the total segment count `K`,
/// the warmup boundary in bases, and the `K` segment byte lengths
/// (`seg_lens[0]` is the warmup segment's).
struct ChunkTable {
    variant: u8,
    segments: usize,
    warm_bases: usize,
    seg_lens: Vec<usize>,
}

/// Read a [`ChunkTable`], refusing an unknown variant tag by name
/// ([`Error::UnsupportedChunkVariant`] — the chunked layout's extension point,
/// with tag `0` reserved for a reset variant no writer emits) and bounding the
/// declared segment count before any per-segment allocation.
fn read_chunk_table(r: &mut ByteReader<'_>) -> Result<ChunkTable> {
    let variant = r.u8()?;
    if variant != CHUNK_VARIANT_WARM_CLONE && variant != CHUNK_VARIANT_WARM_FROZEN {
        return Err(Error::UnsupportedChunkVariant(variant));
    }
    let segments = r.varint()?;
    if !(2..=MAX_CHUNK_SEGMENTS as u64).contains(&segments) {
        return Err(Error::Malformed("chunk segment count out of range"));
    }
    let segments = segments as usize;
    let warm_bases = usize::try_from(r.varint()?)
        .map_err(|_| Error::Malformed("chunk warmup boundary overflows"))?;
    let mut seg_lens = Vec::with_capacity(segments);
    for _ in 0..segments {
        let l = usize::try_from(r.varint()?)
            .map_err(|_| Error::Malformed("chunk segment length overflows"))?;
        seg_lens.push(l);
    }
    Ok(ChunkTable {
        variant,
        segments,
        warm_bases,
        seg_lens,
    })
}

/// Decode a chunked quality payload into `quals` (already sized to the total):
/// recompute the segmentation from `(lens, K)`, decode the warmup serially
/// from a fresh model, snapshot it, then fan the `K - 1` chunks out on the
/// rayon pool — each into its own disjoint output range, so the result is
/// byte-for-byte the serial order at any thread count.
///
/// Decode memory: under warm-clone each **running** chunk task owns one clone
/// of the snapshot (up to ~336 MiB at long-read alphabet depth), and rayon
/// executes at most pool-threads tasks at once — so peak model memory is
/// bounded by `threads × model + one snapshot per in-flight block`, never by
/// `K` (the design doc's "schedule chunks, never spawn all K"). Warm-frozen
/// shares the one snapshot read-only: O(1) model memory at any fan-out.
#[allow(clippy::too_many_arguments)] // mirrors the coder's input list
fn decode_chunked(
    t: &ChunkTable,
    lens: &[u32],
    payload: &[u8],
    seq: &[u8],
    syms: &[u8],
    qmin: u8,
    k: usize,
    qc: &binmix::QCtx,
    quals: &mut [u8],
) -> Result<()> {
    use rayon::prelude::*;
    let total = quals.len();
    // The transmitted segment lengths must tile the payload exactly — a table
    // whose lengths overflow, fall short, or overrun is hostile or corrupt.
    let declared = t
        .seg_lens
        .iter()
        .try_fold(0usize, |a, &l| a.checked_add(l))
        .ok_or(Error::Malformed("chunk segment lengths overflow"))?;
    if declared != payload.len() {
        return Err(Error::Malformed(
            "chunk segment lengths do not match the payload",
        ));
    }
    // The warmup boundary must land exactly on a read boundary — it is a
    // prefix sum of `lens` by construction, and `warmup_reads` re-derives the
    // same minimal read count the encoder used.
    if t.warm_bases > total {
        return Err(Error::Malformed("chunk warmup boundary past the block"));
    }
    let offs = prefix_offsets(lens);
    let wreads = warmup_reads(lens, t.warm_bases);
    if offs[wreads] != t.warm_bases {
        return Err(Error::Malformed(
            "chunk warmup boundary not on a read boundary",
        ));
    }
    let n_chunks = t.segments - 1;
    let bounds = chunk_boundaries(&lens[wreads..], n_chunks);
    // Serial warmup: train a fresh model over the leading reads, exactly as
    // the encoder did; the state at its end is the chunk snapshot.
    let (warm_out, mut chunk_out) = quals.split_at_mut(t.warm_bases);
    let mut mx = binmix::BinMixer::new(k);
    binmix::decode_into(
        &mut binmix::Adapting(&mut mx),
        &lens[..wreads],
        &payload[..t.seg_lens[0]],
        &seq[..t.warm_bases],
        syms,
        qmin,
        qc,
        warm_out,
    )?;
    // Slice each chunk's payload, reads, bases, and disjoint output range.
    let mut tasks = Vec::with_capacity(n_chunks);
    let mut poff = t.seg_lens[0];
    for i in 0..n_chunks {
        let (a, b) = (wreads + bounds[i], wreads + bounds[i + 1]);
        let (head, tail) = std::mem::take(&mut chunk_out).split_at_mut(offs[b] - offs[a]);
        chunk_out = tail;
        let pl = &payload[poff..poff + t.seg_lens[i + 1]];
        poff += t.seg_lens[i + 1];
        tasks.push((pl, &lens[a..b], &seq[offs[a]..offs[b]], head));
    }
    match t.variant {
        CHUNK_VARIANT_WARM_CLONE => tasks.into_par_iter().try_for_each(|(pl, l, s, o)| {
            let mut clone = mx.clone();
            binmix::decode_into(
                &mut binmix::Adapting(&mut clone),
                l,
                pl,
                s,
                syms,
                qmin,
                qc,
                o,
            )
        }),
        CHUNK_VARIANT_WARM_FROZEN => {
            let shared = &mx;
            tasks.into_par_iter().try_for_each(|(pl, l, s, o)| {
                binmix::decode_into(&mut binmix::Frozen(shared), l, pl, s, syms, qmin, qc, o)
            })
        }
        // `read_chunk_table` admits only the two tags above.
        v => Err(Error::UnsupportedChunkVariant(v)),
    }
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
    #[test]
    fn unknown_mode_is_an_upgrade_error_and_fails_safe() {
        // A stream whose mode byte is from the future: decode must refuse it by
        // name (UnsupportedMode, "needs a newer fqxv") rather than as
        // corruption, and needs_sequence must fail SAFE — report true, so the
        // container serializes seq -> qual and the decoder gets to say why it
        // stopped instead of failing a parallel quality-only decode.
        let future = [
            FORMAT_VERSION,
            0,  /* lossless */
            99, /* future mode */
        ];
        assert!(
            needs_sequence(&future),
            "an unknown mode must be assumed sequence-conditioned"
        );
        match decode_seq(&future, &[]) {
            Err(Error::UnsupportedMode(99)) => {}
            other => panic!("want UnsupportedMode(99), got {other:?}"),
        }
        // The known modes keep their meaning: MODE_POS is sequence-free...
        assert!(!needs_sequence(&[FORMAT_VERSION, 0, MODE_POS]));
        // ...the sequence modes still report true...
        for m in [MODE_SEQ, MODE_SEQ_BINMIX, MODE_SEQ_BINMIX_Q] {
            assert!(needs_sequence(&[FORMAT_VERSION, 0, m]));
        }
        // ...and a wrong-version or truncated stream still reads as false.
        assert!(!needs_sequence(&[FORMAT_VERSION + 1, 0, 99]));
        assert!(!needs_sequence(&[FORMAT_VERSION]));
        assert!(!needs_sequence(&[]));
    }

    #[test]
    fn chunk_boundaries_are_a_pure_function_of_lens() {
        let lens: Vec<u32> = (0..1000u32).map(|i| 50 + (i * 37) % 400).collect();
        let total: u64 = lens.iter().map(|&l| u64::from(l)).sum();
        let max_len = u64::from(*lens.iter().max().unwrap());
        for k in [1usize, 2, 4, 8, 16, 32] {
            let a = chunk_boundaries(&lens, k);
            // Deterministic (a pure function of its inputs — nothing else may
            // feed it, or encoder and decoder would disagree on the split).
            assert_eq!(a, chunk_boundaries(&lens, k));
            // Shape: k+1 monotone boundaries covering every read exactly once.
            assert_eq!(a.len(), k + 1);
            assert_eq!(a[0], 0);
            assert_eq!(*a.last().unwrap(), lens.len());
            assert!(a.windows(2).all(|w| w[0] <= w[1]), "monotone boundaries");
            // Balance: a chunk overshoots its base target by at most one read.
            for w in a.windows(2) {
                let bases: u64 = lens[w[0]..w[1]].iter().map(|&l| u64::from(l)).sum();
                assert!(
                    bases <= total / k as u64 + max_len,
                    "chunk [{}, {}) holds {bases} bases (target {})",
                    w[0],
                    w[1],
                    total / k as u64
                );
            }
        }
        // Degenerate inputs must not panic or produce out-of-range boundaries.
        assert_eq!(chunk_boundaries(&[], 8), vec![0usize; 9]);
        assert_eq!(chunk_boundaries(&[10], 4), vec![0, 1, 1, 1, 1]);
    }

    #[test]
    fn qchunk_probes_run_on_small_inputs() {
        // The sweep probes must stay panic-free across every cell even when the
        // block is far smaller than the warmup budgets (all-empty chunks, the
        // warmup swallowing every read, K exceeding the read count). Output is
        // stderr-only; the archive-identity guarantee (flag set vs unset) is
        // enforced end-to-end by bench/scripts/qchunk_diag.sh with cmp.
        let (mut lens, mut quals, mut seq) = (Vec::new(), Vec::new(), Vec::new());
        for r in 0..40u32 {
            let l = 600 + (r as usize * 13) % 100;
            lens.push(l as u32);
            for i in 0..l {
                seq.push(b"ACGT"[(r as usize + i) % 4]);
                quals.push(33 + ((i * 7 + r as usize) % 40) as u8);
            }
        }
        let (_, dense) = dense_alphabet(&quals).expect("alphabet");
        let syms: Vec<u8> = {
            let mut present = [false; 256];
            for &b in &quals {
                present[b as usize] = true;
            }
            (0..=255u8).filter(|&b| present[b as usize]).collect()
        };
        let (qmin, k) = (syms[0], syms.len());
        let flat = binmix::QCtx::flat();
        qchunk_probe_binmix(
            &lens,
            &quals,
            &seq,
            &dense,
            qmin,
            k,
            &flat,
            MODE_SEQ_BINMIX,
            quals.len(),
        );
        dispatch_qchunk_probe_pos(&lens, &quals, &dense, qmin, k, quals.len());
        // Empty input: both probes must return without emitting.
        qchunk_probe_binmix(&[], &[], &[], &dense, qmin, k, &flat, MODE_SEQ_BINMIX, 0);
        dispatch_qchunk_probe_pos(&[], &[], &dense, qmin, k, 0);
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
        let enc = encode_seq(&lens, &quals, &seq, QualityBinning::Lossless, true).expect("encode");
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
        let with_seq =
            encode_seq(&lens, &quals, &seq, QualityBinning::Lossless, true).expect("seq");
        let blind = encode(&lens, &quals, QualityBinning::Lossless).expect("blind");
        assert!(
            with_seq.len() < blind.len(),
            "sequence context ({} B) should beat blind ({} B) on base-correlated quality",
            with_seq.len(),
            blind.len()
        );
    }

    #[test]
    fn quality_trial_keeps_smaller_and_roundtrips() {
        // The per-block quantizer trial is never-worse by construction: the kept
        // stream must be no larger than the pure baseline (MODE_SEQ_BINMIX) stream,
        // and must round-trip whichever candidate wins.
        let (lens, seq, quals) = longread_fixture(40, 2000);
        let (syms, dense) = dense_alphabet(&quals).expect("alphabet");
        let (qmin, k) = (syms[0], syms.len());

        let flat = binmix::QCtx::flat();
        let base_payload = binmix::encode(&lens, &quals, &seq, &dense, qmin, k, &flat);
        let baseline = assemble_quality_stream(
            QualityBinning::Lossless,
            MODE_SEQ_BINMIX,
            k,
            &syms,
            &[],
            &[],
            &lens,
            &base_payload,
        );

        let kept = encode_seq(&lens, &quals, &seq, QualityBinning::Lossless, true).expect("encode");
        assert!(
            kept.len() <= baseline.len(),
            "kept stream ({} B) must not exceed the baseline ({} B)",
            kept.len(),
            baseline.len(),
        );
        assert!(
            matches!(
                kept.get(2),
                Some(&MODE_SEQ_BINMIX) | Some(&MODE_SEQ_BINMIX_Q)
            ),
            "long reads must pick one of the two binmix modes"
        );
        let (out_lens, out_quals) = decode_seq(&kept, &seq).expect("decode");
        assert_eq!(out_lens, lens);
        assert_eq!(out_quals, quals);
    }

    #[test]
    fn quality_trial_is_deterministic() {
        // Encoding the same block twice must be byte-identical — the trial picks a
        // fixed candidate order with a fixed tie-break, so it is order-free.
        let (lens, seq, quals) = longread_fixture(24, 1800);
        let a = encode_seq(&lens, &quals, &seq, QualityBinning::Lossless, true).expect("a");
        let b = encode_seq(&lens, &quals, &seq, QualityBinning::Lossless, true).expect("b");
        assert_eq!(a, b, "encode_seq must be deterministic");
    }

    #[test]
    fn quant_mode_roundtrips_directly() {
        // Exercise the MODE_SEQ_BINMIX_Q encode+decode path explicitly, so it is
        // covered even when the trial happens not to select it on a fixture.
        let (lens, seq, quals) = longread_fixture(30, 1500);
        let (syms, dense) = dense_alphabet(&quals).expect("alphabet");
        let (qmin, k) = (syms[0], syms.len());
        let (qc, qtable) = build_quant_ctx(&quals, &syms, qmin);
        let payload = binmix::encode(&lens, &quals, &seq, &dense, qmin, k, &qc);
        let stream = assemble_quality_stream(
            QualityBinning::Lossless,
            MODE_SEQ_BINMIX_Q,
            k,
            &syms,
            &qtable,
            &[],
            &lens,
            &payload,
        );
        assert!(needs_sequence(&stream), "quant mode needs the sequence");
        assert_eq!(stream.get(2), Some(&MODE_SEQ_BINMIX_Q));
        let (out_lens, out_quals) = decode_seq(&stream, &seq).expect("decode");
        assert_eq!(out_lens, lens);
        assert_eq!(out_quals, quals);
        // A stream truncated inside the quantizer table must be rejected cleanly,
        // not panic: header is version|binning|mode|k|syms[k], then the 2*k table.
        let cut = 4 + k + k; // stop partway through the 2*k-byte table
        assert!(decode_seq(&stream[..cut], &seq).is_err());
    }

    #[test]
    fn seq_context_decode_requires_sequence() {
        // A MODE_SEQ stream cannot be decoded blind: the plain `decode` (empty
        // sequence) must reject it cleanly rather than produce wrong output.
        let (lens, seq, quals) = longread_fixture(20, 1000);
        let enc = encode_seq(&lens, &quals, &seq, QualityBinning::Lossless, true).expect("encode");
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
        let with_seq =
            encode_seq(&lens, &quals, &seq, QualityBinning::Lossless, true).expect("seq");
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

    // --- chunked quality (MODE_SEQ_BINMIX_CHUNKED / _Q) ----------------------

    /// A chunking request usable on the small in-crate fixtures (the public
    /// gate demands `MIN_CHUNKED_BASES`, so tests drive the internal
    /// entry points directly to keep them fast).
    fn chunking(segments: usize, warmup: ChunkWarmup, variant: ChunkVariant) -> QualityChunking {
        QualityChunking {
            segments,
            warmup,
            variant,
        }
    }

    /// Chunked-encode a fixture through the internal (gate-free) path.
    fn encode_chunked_fixture(
        lens: &[u32],
        seq: &[u8],
        quals: &[u8],
        try_quantizer: bool,
        c: QualityChunking,
    ) -> Vec<u8> {
        let (syms, dense) = dense_alphabet(quals).expect("alphabet");
        let (qmin, k) = (syms[0], syms.len());
        encode_binmix_chunked(
            lens,
            quals,
            seq,
            &dense,
            qmin,
            k,
            &syms,
            QualityBinning::Lossless,
            try_quantizer,
            c,
        )
        .expect("chunked encode")
    }

    #[test]
    fn chunked_roundtrips_across_variants_ks_and_warmups() {
        let (lens, seq, quals) = longread_fixture(60, 900);
        for segments in [2usize, 3, 8, 16] {
            for warmup in [
                ChunkWarmup::TotalOverSegments,
                ChunkWarmup::Bases(10_000),
                ChunkWarmup::Bases(0), // degenerate: one-read warmup
            ] {
                for variant in [ChunkVariant::WarmClone, ChunkVariant::WarmFrozen] {
                    let stream = encode_chunked_fixture(
                        &lens,
                        &seq,
                        &quals,
                        false,
                        chunking(segments, warmup, variant),
                    );
                    assert_eq!(
                        stream.get(2),
                        Some(&MODE_SEQ_BINMIX_CHUNKED),
                        "flat chunked stream must carry mode 5"
                    );
                    assert!(needs_sequence(&stream), "chunked modes need the sequence");
                    let (out_lens, out_quals) = decode_seq(&stream, &seq).unwrap_or_else(|e| {
                        panic!("decode K={segments} {warmup:?} {variant:?}: {e}")
                    });
                    assert_eq!(out_lens, lens);
                    assert_eq!(out_quals, quals, "K={segments} {warmup:?} {variant:?}");
                }
            }
        }
    }

    #[test]
    fn chunked_quantizer_trial_roundtrips_and_is_deterministic() {
        let (lens, seq, quals) = longread_fixture(50, 1200);
        let c = chunking(4, ChunkWarmup::TotalOverSegments, ChunkVariant::WarmClone);
        let a = encode_chunked_fixture(&lens, &seq, &quals, true, c);
        let b = encode_chunked_fixture(&lens, &seq, &quals, true, c);
        assert_eq!(a, b, "chunked encode must be deterministic");
        assert!(
            matches!(
                a.get(2),
                Some(&MODE_SEQ_BINMIX_CHUNKED) | Some(&MODE_SEQ_BINMIX_CHUNKED_Q)
            ),
            "the chunked trial must keep one of the two chunked modes"
        );
        let (out_lens, out_quals) = decode_seq(&a, &seq).expect("decode");
        assert_eq!(out_lens, lens);
        assert_eq!(out_quals, quals);
    }

    #[test]
    fn chunked_q_mode_roundtrips_directly() {
        // Pin the mode-6 wire path explicitly (qtable + chunk table together),
        // independent of whether the trial would select it on this fixture.
        let (lens, seq, quals) = longread_fixture(40, 1000);
        let (syms, dense) = dense_alphabet(&quals).expect("alphabet");
        let (qmin, k) = (syms[0], syms.len());
        let (qc, qtable) = build_quant_ctx(&quals, &syms, qmin);
        let c = chunking(3, ChunkWarmup::Bases(8_000), ChunkVariant::WarmClone);
        let (chunk_table, payload) =
            encode_chunked_payload(&lens, &quals, &seq, &dense, qmin, k, &qc, c);
        let stream = assemble_quality_stream(
            QualityBinning::Lossless,
            MODE_SEQ_BINMIX_CHUNKED_Q,
            k,
            &syms,
            &qtable,
            &chunk_table,
            &lens,
            &payload,
        );
        assert_eq!(stream.get(2), Some(&MODE_SEQ_BINMIX_CHUNKED_Q));
        assert!(needs_sequence(&stream));
        let (out_lens, out_quals) = decode_seq(&stream, &seq).expect("decode");
        assert_eq!(out_lens, lens);
        assert_eq!(out_quals, quals);
    }

    #[test]
    fn chunked_bytes_and_decode_are_thread_count_invariant() {
        // The hard repo invariant, at the codec level: the same input must
        // produce identical chunked bytes — and decode back identically —
        // whether rayon has 1 worker or several.
        let (lens, seq, quals) = longread_fixture(48, 1100);
        let c = chunking(8, ChunkWarmup::Bases(12_000), ChunkVariant::WarmClone);
        let one = rayon::ThreadPoolBuilder::new()
            .num_threads(1)
            .build()
            .expect("pool");
        let four = rayon::ThreadPoolBuilder::new()
            .num_threads(4)
            .build()
            .expect("pool");
        let enc1 = one.install(|| encode_chunked_fixture(&lens, &seq, &quals, false, c));
        let enc4 = four.install(|| encode_chunked_fixture(&lens, &seq, &quals, false, c));
        assert_eq!(enc1, enc4, "chunked encode must be thread-count invariant");
        let dec1 = one.install(|| decode_seq(&enc1, &seq).expect("decode@1"));
        let dec4 = four.install(|| decode_seq(&enc1, &seq).expect("decode@4"));
        assert_eq!(dec1, dec4, "chunked decode must be thread-count invariant");
        assert_eq!(dec1.1, quals);
    }

    #[test]
    fn chunked_gate_falls_back_to_serial_below_minimum() {
        // Under MIN_CHUNKED_BASES the public entry point must ignore the
        // chunking request entirely — byte-identical to the serial encode.
        let (lens, seq, quals) = longread_fixture(40, 2000);
        assert!(quals.len() < MIN_CHUNKED_BASES);
        let serial = encode_seq(&lens, &quals, &seq, QualityBinning::Lossless, true).expect("s");
        let chunked_req = encode_seq_chunked(
            &lens,
            &quals,
            &seq,
            QualityBinning::Lossless,
            true,
            Some(chunking(
                8,
                ChunkWarmup::Bases(8 << 20),
                ChunkVariant::WarmClone,
            )),
        )
        .expect("c");
        assert_eq!(
            chunked_req, serial,
            "an ineligible chunking request must be a no-op"
        );
    }

    #[test]
    fn chunking_eligibility_bounds() {
        // Only the read lengths matter to eligibility; keep `total` consistent
        // with them (as the real caller guarantees) without materializing MiB
        // of data.
        let lens = vec![600_000u32; 40]; // 24 MB of quality across 40 reads
        let total = 40 * 600_000usize;
        assert!(total >= MIN_CHUNKED_BASES);
        let ok = |c: QualityChunking| chunking_eligible(&lens, total, &c);
        // Segment count out of range.
        for segments in [0usize, 1, MAX_CHUNK_SEGMENTS + 1] {
            assert!(!ok(chunking(
                segments,
                ChunkWarmup::TotalOverSegments,
                ChunkVariant::WarmClone
            )));
        }
        // More segments than reads.
        assert!(!ok(chunking(
            lens.len() + 1,
            ChunkWarmup::TotalOverSegments,
            ChunkVariant::WarmClone
        )));
        // A warmup that swallows every read leaves nothing to fan.
        assert!(!ok(chunking(
            4,
            ChunkWarmup::Bases(usize::MAX),
            ChunkVariant::WarmClone
        )));
        // Too few total bases.
        let small = vec![100_000u32; 40]; // 4 MB — under the floor
        assert!(!chunking_eligible(
            &small,
            40 * 100_000,
            &chunking(4, ChunkWarmup::TotalOverSegments, ChunkVariant::WarmClone)
        ));
        // The shapes the encoder actually emits are eligible.
        assert!(ok(chunking(
            8,
            ChunkWarmup::TotalOverSegments,
            ChunkVariant::WarmClone
        )));
        assert!(ok(chunking(
            8,
            ChunkWarmup::Bases(8 << 20),
            ChunkVariant::WarmClone
        )));
    }

    #[test]
    fn chunked_zero_length_reads_at_the_warm_boundary_roundtrip() {
        // Zero-length reads make prefix sums repeat; the warm boundary and the
        // chunk boundaries must stay consistent between encode and decode.
        let (mut lens, seq, quals) = longread_fixture(30, 800);
        // Splice zero-length reads around the warmup landing zone and the end.
        lens.insert(1, 0);
        lens.insert(5, 0);
        lens.push(0);
        for variant in [ChunkVariant::WarmClone, ChunkVariant::WarmFrozen] {
            let stream = encode_chunked_fixture(
                &lens,
                &seq,
                &quals,
                false,
                chunking(4, ChunkWarmup::Bases(800), variant),
            );
            let (out_lens, out_quals) =
                decode_seq(&stream, &seq).unwrap_or_else(|e| panic!("{variant:?}: {e}"));
            assert_eq!(out_lens, lens);
            assert_eq!(out_quals, quals);
        }
    }

    /// Build a syntactically-valid chunked header with an arbitrary table, for
    /// the adversarial parse tests: `version|binning|mode=5|k=1|sym` + `table`
    /// + a lens array of `n` fixed-length-1 reads + `payload`.
    fn chunked_stream_with_table(table: &[u8], n_reads: u64, payload: &[u8]) -> Vec<u8> {
        let mut s = vec![FORMAT_VERSION, 0, MODE_SEQ_BINMIX_CHUNKED, 1, b'I'];
        s.extend_from_slice(table);
        push_varint(&mut s, n_reads); // n reads...
        s.push(1); // ...fixed-length...
        push_varint(&mut s, 1); // ...of 1 base each
        s.extend_from_slice(payload);
        s
    }

    #[test]
    fn chunked_rejects_hostile_tables() {
        let n = 64u64; // 64 reads of 1 base
        let seq = vec![b'A'; 64];
        let table = |variant: u8, k: u64, warm_bases: u64, seg_lens: &[u64]| {
            let mut t = vec![variant];
            push_varint(&mut t, k);
            push_varint(&mut t, warm_bases);
            for &l in seg_lens {
                push_varint(&mut t, l);
            }
            t
        };
        // Reserved / unknown variant tags refuse by name, before anything else.
        for v in [0u8, 3, 99] {
            let s = chunked_stream_with_table(&table(v, 2, 8, &[4, 4]), n, &[0; 8]);
            match decode_seq(&s, &seq) {
                Err(Error::UnsupportedChunkVariant(got)) => assert_eq!(got, v),
                other => panic!("variant {v}: want UnsupportedChunkVariant, got {other:?}"),
            }
        }
        // Segment counts out of range: 0, 1, and absurd.
        for k in [0u64, 1, MAX_CHUNK_SEGMENTS as u64 + 1, u64::MAX >> 8] {
            let s = chunked_stream_with_table(&table(1, k, 8, &[4, 4]), n, &[0; 8]);
            assert!(
                matches!(decode_seq(&s, &seq), Err(Error::Malformed(_))),
                "K={k} must be rejected"
            );
        }
        // A table truncated mid-list (declares 4 segments, carries 1 length).
        let mut t = vec![1u8];
        push_varint(&mut t, 4);
        push_varint(&mut t, 8);
        push_varint(&mut t, 4);
        let mut s = vec![FORMAT_VERSION, 0, MODE_SEQ_BINMIX_CHUNKED, 1, b'I'];
        s.extend_from_slice(&t);
        assert!(matches!(decode_seq(&s, &seq), Err(Error::Malformed(_))));
        // Segment lengths that overflow usize when summed.
        let s = chunked_stream_with_table(
            &table(1, 3, 8, &[u64::MAX / 2, u64::MAX / 2, 64]),
            n,
            &[0; 8],
        );
        assert!(matches!(decode_seq(&s, &seq), Err(Error::Malformed(_))));
        // Lengths that do not tile the payload exactly (short and long).
        for lens2 in [[4u64, 3], [4, 5]] {
            let s = chunked_stream_with_table(&table(1, 2, 8, &lens2), n, &[0; 8]);
            assert!(matches!(decode_seq(&s, &seq), Err(Error::Malformed(_))));
        }
        // A warmup boundary past the block.
        let s = chunked_stream_with_table(&table(1, 2, 65, &[4, 4]), n, &[0; 8]);
        assert!(matches!(decode_seq(&s, &seq), Err(Error::Malformed(_))));
    }

    #[test]
    fn chunked_rejects_warm_boundary_off_a_read() {
        // 64 reads of 3 bases: a warm boundary of 4 lands inside read 1.
        let seq = vec![b'A'; 192];
        let mut s = vec![FORMAT_VERSION, 0, MODE_SEQ_BINMIX_CHUNKED, 1, b'I'];
        s.push(1); // warm-clone
        push_varint(&mut s, 2); // K = 2
        push_varint(&mut s, 4); // warm boundary: not a multiple of 3
        push_varint(&mut s, 4); // warm segment length
        push_varint(&mut s, 4); // chunk length
        push_varint(&mut s, 64); // 64 reads...
        s.push(1); // ...fixed-length...
        push_varint(&mut s, 3); // ...of 3 bases
        s.extend_from_slice(&[0; 8]);
        match decode_seq(&s, &seq) {
            Err(Error::Malformed(m)) => assert!(m.contains("read boundary"), "{m}"),
            other => panic!("want Malformed(read boundary), got {other:?}"),
        }
    }

    #[test]
    fn needs_sequence_knows_the_chunked_modes() {
        for m in [MODE_SEQ_BINMIX_CHUNKED, MODE_SEQ_BINMIX_CHUNKED_Q] {
            assert!(needs_sequence(&[FORMAT_VERSION, 0, m]));
        }
    }

    proptest::proptest! {
        // Chunked round-trip over arbitrary read shapes (including empty
        // reads), K, and both variants — decode(encode(x)) == x with the
        // segmentation recomputed on the decode side.
        #[test]
        fn chunked_roundtrip_arbitrary(
            reads in proptest::collection::vec(
                proptest::collection::vec(33u8..=90, 0..60), 4..24),
            segments in 2usize..6,
            frozen in proptest::bool::ANY,
        ) {
            let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
            let quals: Vec<u8> = reads.concat();
            if quals.is_empty() {
                return Ok(());
            }
            let seq: Vec<u8> = quals.iter().map(|&q| b"ACGT"[q as usize % 4]).collect();
            let variant = if frozen { ChunkVariant::WarmFrozen } else { ChunkVariant::WarmClone };
            let stream = encode_chunked_fixture(
                &lens, &seq, &quals, false,
                chunking(segments, ChunkWarmup::TotalOverSegments, variant),
            );
            let (out_lens, out_quals) = decode_seq(&stream, &seq).expect("decode");
            proptest::prop_assert_eq!(out_lens, lens);
            proptest::prop_assert_eq!(out_quals, quals);
        }
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
