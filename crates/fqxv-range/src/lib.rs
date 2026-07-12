//! Binary range coder and adaptive frequency models.
//!
//! This is the serial entropy backend for the fqzcomp quality model
//! ([`fqxv-fqzcomp`](https://docs.rs/fqxv-fqzcomp)). The range coder follows the
//! public-domain design used by the CRAM codecs; the adaptive models are the
//! `SIMPLE_MODEL` frequency tables described there.
//!
//! Status: **scaffold.** Implemented in M2.

use thiserror::Error;

/// Errors returned by the range coder.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// The compressed stream was malformed or truncated.
    #[error("malformed range-coded stream: {0}")]
    Malformed(&'static str),
    /// A code path that is not yet implemented in this scaffold.
    #[error("not yet implemented: {0}")]
    NotImplemented(&'static str),
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;
