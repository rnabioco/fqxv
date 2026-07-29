//! Edit-distance sequence alignment primitives.
//!
//! Two implementations of the same cost model — **unit edit costs**
//! (Levenshtein: match 0, substitution / insertion / deletion 1 each) — so
//! results are directly comparable, and both report an edit script plus a
//! distance:
//!
//! - [`align_banded`] — banded Needleman-Wunsch. Work is `O(n · band)`
//!   regardless of how similar the inputs are, with an AVX2 anti-diagonal
//!   backend and a scalar fallback. Predictable cost, and the band is the
//!   knob that bounds it.
//! - [`wfa_align`] / [`wfa_align_opt`] — wavefront alignment (Marco-Sola et
//!   al. 2021 recurrences, clean-room). Work scales with the alignment
//!   **score** rather than sequence length, so it is far cheaper on
//!   low-divergence inputs and degrades toward the DP as divergence rises.
//!   Traceback storage is `O(s²)`.
//!
//! The two deliberately choose different (equal-cost) paths, so their scripts
//! are not byte-identical; only the distances agree.
//!
//! [`wfa_align_opt`] returns `None` once the score exceeds a cap, which makes
//! it usable as a bounded "is this within N edits?" test that abandons hopeless
//! pairs early rather than aligning them to completion.
//!
//! # Why this is a separate crate
//!
//! These primitives began inside `fqxv-lroverlap`, where they align the gaps
//! between chained anchors. They carry no codec state and no dependencies, so
//! they live here as a leaf crate that consumers can take on their own.
//! `fqxv-lroverlap` re-exports the same names, so nothing downstream of it
//! needs to change.

// `deny` rather than `forbid`: this crate is unsafe-free except for the AVX2
// alignment backend in `align`, which opts back in at three narrowly annotated
// sites (`#[allow(unsafe_code)]`) and is proptested byte-identical to scalar.
#![deny(unsafe_code)]

mod align;
mod wfa;

pub use align::{Alignment, Op, align_banded, apply};
pub use wfa::{wfa_align, wfa_align_opt, wfa_cells};
