//! Reference-free read reordering (PgRC2/SPRING-class).
//!
//! Builds a minimizer/k-mer index over the reads, assembles an approximate
//! shortest-common-superstring "pseudogenome" from high-quality reads, then maps
//! the remaining reads (reverse-complement aware) with a bounded mismatch
//! budget. Emits a reordered residual stream for `fqxv-seq`/`fqxv-rans` plus a
//! permutation so the original order can be restored (order-preserving mode).
//!
//! In-memory and `rayon`-parallel. This is the largest and highest-risk crate;
//! the rest of the toolkit ships without it.
//!
//! Status: **scaffold.** Implemented in M4.

use thiserror::Error;

/// Whether to preserve the original read order on decompression.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Ordering {
    /// Restore the exact input order (stores a permutation).
    #[default]
    Preserve,
    /// Allow reads to be emitted in reordered form (smaller, order not kept).
    Free,
}

/// Errors returned by the reordering engine.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// The compressed stream was malformed or truncated.
    #[error("malformed reorder stream: {0}")]
    Malformed(&'static str),
    /// A code path that is not yet implemented in this scaffold.
    #[error("not yet implemented: {0}")]
    NotImplemented(&'static str),
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;
