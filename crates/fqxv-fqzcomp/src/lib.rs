//! fqzcomp-style quality-score context model.
//!
//! The per-symbol context is assembled from quality history, position, a
//! running "how noisy so far" delta counter, and an optional selector, then fed
//! to the [`fqxv-range`](https://docs.rs/fqxv-range) coder. Optional remap
//! tables (`qtab`/`ptab`/`dtab`) match the CRAM 3.1 fqzcomp parameter space.
//!
//! Lossy quality binning (Illumina 2/4/8-bin) is opt-in via [`QualityBinning`]
//! and applied before modeling; the default is lossless.
//!
//! Status: **scaffold.** Implemented in M2.

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

/// Errors returned by the quality codec.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// The compressed stream was malformed or truncated.
    #[error("malformed fqzcomp stream: {0}")]
    Malformed(&'static str),
    /// Underlying range coder failure.
    #[error(transparent)]
    Range(#[from] fqxv_range::Error),
    /// A code path that is not yet implemented in this scaffold.
    #[error("not yet implemented: {0}")]
    NotImplemented(&'static str),
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;
