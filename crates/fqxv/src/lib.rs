//! `fqxv` — reference-free FASTQ archiver.
//!
//! Defines the block-based container format and composes the codec crates
//! ([`fqxv_tokenizer`] for names, [`fqxv_seq`] for sequences, [`fqxv_fqzcomp`]
//! for qualities) into a full compressor. Blocks are compressed and decompressed
//! in parallel with `rayon`; output is deterministic regardless of thread count.
//!
//! Losslessness: read name + description, sequence, and quality are preserved
//! exactly. The `+` separator line is normalized to a bare `+` (as SPRING and
//! fqz_comp do) — its optional repeated header is not retained.
//!
//! The binary in `fqxv-cli` is a thin CLI over [`compress`], [`decompress`],
//! and [`inspect`].

mod container;
mod crc;

pub use container::{
    BlockContents, ContentStats, Estimate, GroupLoc, Index, Info, Params, Platform, QUAL_MAX,
    Record, RecordReader, Recovery, Stats, Stream, SuffixParse, VerifyCheck, VerifyReport,
    compress, compress_auto, compress_interleaved, compress_multi, content_stats,
    decode_block_contents, decode_names, decode_quality, decode_quality_with_seq, decode_sequence,
    decompress, decompress_records, decompress_recover, decompress_split, estimate, expected_reads,
    inspect, peek, quality_needs_sequence, verify, verify_quick, verify_report, verify_roundtrip,
};
pub use fqxv_fqzcomp::QualityBinning;

use thiserror::Error;

/// The four-byte magic at the head of every `.fqxv` file.
pub const MAGIC: [u8; 4] = *b"FQXV";

/// The container format **major** version this build writes.
///
/// **The on-disk format is stable at 1.0.** An archive written by any format-1.x
/// release stays readable by later ones. A reader refuses an archive whose major
/// differs from its own (wire-incompatible) but tolerates a newer
/// [`FORMAT_MINOR`] within its major, since minor bumps only add things an older
/// reader can either skip or refuse deliberately.
///
/// This version is **independent of the crate version**: the crate is pre-1.0 and
/// its Rust API may still change, while the bytes on disk are covered by the
/// guarantee above. `fqxv --version` reports the former, `fqxv info` the latter.
///
/// A major bump is a last resort, not a planned step. The additive mechanisms —
/// 63 free [`feature`] bits, a 64 KiB extension region with 127 free non-critical
/// and 128 free critical tags, and 251 unused sequence method-byte values — are
/// meant to carry foreseeable evolution without one. Should a format 2 ever
/// happen, readers will keep a format-1 read path rather than orphan existing
/// archives. See `docs/design/container.md` for the full evolution policy: which
/// mechanism a given change belongs in, and the rule that any change to the block
/// payload's stream layout or the footer's shape must set a [`feature`] bit so an
/// older reader refuses at the header instead of mis-reading the body.
///
/// (The "v2/v3/v4" you see elsewhere are the per-block *sequence codec* versions,
/// independent of this container version.)
///
/// The v1 container:
/// - appends a footer index (`[u32 n_row_groups]`, per-group
///   `[u64 offset][u32 read_count]` followed by a `[u64 offset][u32 len][u32
///   crc32c]` triple for each of the three streams (names, sequence, quality),
///   `[u64 total_reads]`, `[u32 whole_file_crc]`, `[u32 footer_crc]`) plus an EOF
///   trailer (`[u64 footer_offset]["FQXF"]`) after a zero-length terminator block,
///   so `inspect` and random access can seek straight to the row-group index — and
///   a remote client can project a single stream (e.g. read names, ~1% of the
///   archive) with one range request, verifying it against the per-stream CRC;
/// - carries a CRC-32C on every coded payload, so on-disk corruption is detected
///   and localized rather than silently decoded;
/// - prepends three xxh3-64 digests of each block's *decoded* content — one per
///   stream (names, sequence, post-binning quality) — to the block payload,
///   verified after decode, so a codec bug that turns CRC-valid bytes into
///   wrong-but-in-bounds output is caught at runtime and localized to the stream
///   that regressed, not just to the block;
/// - tags the `FLAG_GLOBAL_REFERENCE` frame with a leading method byte so the
///   shared reference can be coded by either the clean-room order-k model or an
///   xz pass (whichever is smaller), exploiting long-range repeat structure the
///   order-k model can't see;
/// - prefixes every block frame with a `BLOCK_MAGIC` sync marker so
///   [`decompress_recover`] can resynchronize to a block boundary by scanning —
///   recovering intact blocks even when the footer index is lost (a truncated
///   tail) or a block's length prefix is corrupt;
/// - tags each block's sequence stream with a leading method byte so the
///   codec is chosen per block: the order-k context model for short reads, the
///   cross-read overlap-assembly codec (`fqxv-lroverlap`) for long reads.
///
/// See `container.rs` for the full layout.
pub const FORMAT_MAJOR: u8 = 1;

/// On-disk format minor version. A reader tolerates any minor within its
/// [`FORMAT_MAJOR`]; nothing in decode branches on it, so the value is purely
/// informational (surfaced by `inspect` and `verify`).
///
/// Bump it only when a reader can *act* on the difference — a new [`feature`] bit
/// it must refuse, or a new critical extension record. A purely skippable addition
/// (a non-critical extension tag an older reader ignores while decoding every byte
/// identically) does **not** need one: the extension region is self-describing, so
/// the bump would fire on every archive, including the ones carrying none of the
/// new field, while telling a reader nothing it can use.
pub const FORMAT_MINOR: u8 = 0;

/// Coarse capability bits an archive can require in its header `required_features`
/// word. A reader that sees a set bit outside [`KNOWN_FEATURES`] refuses the
/// archive with [`Error::UnsupportedFeature`] — a precise "upgrade fqxv" signal
/// rather than a mis-decode. Set a bit only for a capability a reader must have to
/// decode the archive *at all*, and that is knowable before the blocks are
/// written; per-block codec choices are gated by the sequence stream's method byte
/// instead ([`Error::UnsupportedMethod`]).
pub mod feature {
    /// The whole-file reorder layout carries a shared global reference frame
    /// (SPRING-style); a reader without that decode path cannot reconstruct the
    /// referenced sequence blocks.
    pub const GLOBAL_REFERENCE: u64 = 1 << 0;
}

/// The union of every [`feature`] bit this build understands. A `required_features`
/// word with any bit outside this mask is rejected by
/// `read_header`.
pub const KNOWN_FEATURES: u64 = feature::GLOBAL_REFERENCE;

/// Errors returned by the archiver.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// I/O failure.
    #[error(transparent)]
    Io(#[from] std::io::Error),
    /// The file was not a recognizable `.fqxv` container.
    #[error("not an fqxv container (bad magic)")]
    BadMagic,
    /// The container format major version differs from this build's
    /// [`FORMAT_MAJOR`] — the archive is wire-incompatible and cannot be read.
    #[error(
        "unsupported fqxv format version {major}.{minor} (this build reads {}.x)",
        FORMAT_MAJOR
    )]
    UnsupportedVersion {
        /// On-disk major version (differs from [`FORMAT_MAJOR`]).
        major: u8,
        /// On-disk minor version (informational).
        minor: u8,
    },
    /// The archive's header requires one or more capabilities ([`feature`]) this
    /// build doesn't have. The payload is the set of unknown required bits.
    #[error("fqxv archive requires unsupported features (bits {0:#x}); upgrade fqxv")]
    UnsupportedFeature(u64),
    /// The header carried a critical extension record with a tag this build
    /// doesn't recognize (see the header extension region). The payload is the tag.
    #[error("fqxv archive has an unsupported critical header extension (tag {0:#x}); upgrade fqxv")]
    UnsupportedExtension(u8),
    /// The header set one or more layout flag bits outside this build's known set.
    /// Every flag changes how the archive is read, so an unknown one is refused
    /// rather than ignored. The payload is the set of unknown bits.
    #[error("fqxv archive uses unsupported header flags (bits {0:#x}); upgrade fqxv")]
    UnsupportedFlags(u8),
    /// A stream was tagged with a codec method byte this build can't decode (e.g.
    /// a sequence codec added in a newer minor). Localized to the stream and method.
    #[error("unsupported fqxv {stream} codec method {method}; upgrade fqxv")]
    UnsupportedMethod {
        /// Which stream carried the unknown method byte (e.g. `"sequence"`).
        stream: &'static str,
        /// The unrecognized method byte.
        method: u8,
    },
    /// The stream ended in the middle of a block.
    #[error("truncated fqxv stream")]
    Truncated,
    /// A CRC-32C check failed: the data was corrupted on disk or in transit.
    /// `what` names the region (e.g. `"block 12"`, `"footer"`).
    #[error("corrupted fqxv data: {what} (crc mismatch)")]
    Corrupt {
        /// Human-readable name of the region whose checksum failed.
        what: String,
    },
    /// A block referenced a codec parameter this build doesn't support.
    #[error("malformed fqxv block: {0}")]
    Malformed(&'static str),
    /// A FASTQ record's sequence and quality lines differ in length. Valid FASTQ
    /// requires them equal; caught per-record at parse time so a compensating
    /// pair of mismatches (one read long on sequence, another long on quality)
    /// can't net out at the block level and silently misalign the quality stream.
    #[error("invalid FASTQ: sequence length {seq} but quality length {qual} (must be equal)")]
    RecordLengthMismatch {
        /// The record's sequence-line length.
        seq: usize,
        /// The record's quality-line length.
        qual: usize,
    },
    /// Read-name codec failure.
    #[error(transparent)]
    Tokenizer(#[from] fqxv_tokenizer::Error),
    /// Sequence codec failure.
    #[error(transparent)]
    Seq(#[from] fqxv_seq::Error),
    /// Quality codec failure.
    #[error(transparent)]
    Fqzcomp(#[from] fqxv_fqzcomp::Error),
    /// Read-reordering codec failure.
    #[error(transparent)]
    Reorder(#[from] fqxv_reorder::Error),
    /// Long-read overlap sequence codec failure.
    #[error(transparent)]
    Lroverlap(#[from] fqxv_lroverlap::Error),
    /// rANS coder failure (permutation stream).
    #[error(transparent)]
    Rans(#[from] fqxv_rans::Error),
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;
