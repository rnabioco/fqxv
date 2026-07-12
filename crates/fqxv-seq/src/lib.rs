//! Nucleotide sequence packing and order-k context modeling.
//!
//! Provides 2-bit packing for pure-ACGT reads and a self-synchronizing
//! variable-length packing that tolerates `N` and IUPAC codes, plus an order-k
//! context model entropy-coded with [`fqxv-rans`](https://docs.rs/fqxv-rans).
//! This is the sequence path used when reads are not reordered (the reordered
//! path lives in `fqxv-reorder`).
//!
//! Status: **scaffold.** Implemented in M3.

use thiserror::Error;

/// A 2-bit nucleotide code.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[repr(u8)]
pub enum Base {
    /// Adenine.
    A = 0,
    /// Cytosine.
    C = 1,
    /// Guanine.
    G = 2,
    /// Thymine.
    T = 3,
}

/// Errors returned by the sequence codec.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// The compressed stream was malformed or truncated.
    #[error("malformed sequence stream: {0}")]
    Malformed(&'static str),
    /// Underlying rANS coder failure.
    #[error(transparent)]
    Rans(#[from] fqxv_rans::Error),
    /// A code path that is not yet implemented in this scaffold.
    #[error("not yet implemented: {0}")]
    NotImplemented(&'static str),
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;
