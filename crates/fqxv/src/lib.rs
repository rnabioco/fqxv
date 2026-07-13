//! `fqxv` — reference-free FASTQ archiver for short-read data.
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
    compress, compress_auto, compress_interleaved, compress_multi, decompress, decompress_recover,
    decompress_split, inspect, peek, verify, Info, Params, Recovery, Stats,
};
pub use fqxv_fqzcomp::QualityBinning;

use thiserror::Error;

/// The four-byte magic at the head of every `.fqxv` file.
pub const MAGIC: [u8; 4] = *b"FQXV";

/// The container format version this build writes.
///
/// v1 appends a footer index (`[u32 n_row_groups]`, per-group `[u64 offset]
/// [u32 read_count]`, `[u64 total_reads]`, `[u32 whole_file_crc]`,
/// `[u32 footer_crc]`) plus an EOF trailer (`[u64 footer_offset]["FQXF"]`) after a
/// zero-length terminator block, so `inspect` and random access can seek straight
/// to the row-group index. Every coded payload carries a CRC-32C so corruption is
/// detected and localized rather than silently decoded; see `container.rs` for
/// the full layout. Nothing on disk is stable yet (alpha).
pub const FORMAT_VERSION: u16 = 1;

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
    /// The container format version is newer than this build understands.
    #[error("unsupported fqxv format version {0}")]
    UnsupportedVersion(u16),
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
    /// rANS coder failure (permutation stream).
    #[error(transparent)]
    Rans(#[from] fqxv_rans::Error),
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;
