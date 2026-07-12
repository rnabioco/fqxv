//! Positional read-name tokenizer for FASTQ headers.
//!
//! Read names are split into a positional sequence of typed tokens (alpha runs,
//! digits, zero-padded digits, single chars) and each column is modeled against
//! the previous record's token at the same position — so incrementing
//! instrument/tile/x/y fields collapse to matches and small deltas. The
//! tokenized streams are entropy-coded with
//! [`fqxv-rans`](https://docs.rs/fqxv-rans).
//!
//! Status: **scaffold.** Implemented in M3.

use thiserror::Error;

/// A single token in the positional decomposition of a read name.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum TokenType {
    /// Identical to the previous record's token at this position.
    Match,
    /// A run of alphabetic/punctuation bytes.
    Alpha,
    /// A single byte.
    Char,
    /// A number beginning 1-9.
    Digits,
    /// A zero-prefixed number (value plus width, to preserve leading zeros).
    Digits0,
    /// A small numeric delta versus the previous record.
    Delta,
    /// A small numeric delta on a zero-padded field.
    Delta0,
}

/// Errors returned by the tokenizer codec.
#[derive(Debug, Error)]
#[non_exhaustive]
pub enum Error {
    /// The compressed stream was malformed or truncated.
    #[error("malformed name stream: {0}")]
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
