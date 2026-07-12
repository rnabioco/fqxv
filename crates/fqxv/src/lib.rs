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

pub use container::{compress, decompress, inspect, Info, Params, Stats};
pub use fqxv_fqzcomp::QualityBinning;

use thiserror::Error;

/// The four-byte magic at the head of every `.fqxv` file.
pub const MAGIC: [u8; 4] = *b"FQXV";

/// The container format version this build writes.
pub const FORMAT_VERSION: u16 = 0;

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
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;
