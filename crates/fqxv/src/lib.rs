//! `fqxv` — reference-free FASTQ archiver for short-read data.
//!
//! This crate defines the block-based container format and composes the codec
//! crates (`fqxv-rans`, `fqxv-range`, `fqxv-fqzcomp`, `fqxv-tokenizer`,
//! `fqxv-seq`, `fqxv-reorder`) into a full compressor. The binary in
//! `src/main.rs` is a thin CLI over this library.
//!
//! Status: **scaffold.** The container format lands in M5, after the codecs
//! prove out in the `bench/` harness.

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
    /// A code path that is not yet implemented in this scaffold.
    #[error("not yet implemented: {0}")]
    NotImplemented(&'static str),
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;
