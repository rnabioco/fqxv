//! Cross-read overlap detection for long-read (ONT / PacBio) sequence coding.
//!
//! `fqxv-seq` models each read independently with an order-k context model. At
//! long-read coverage that leaves the dominant redundancy on the table: at ~300x
//! the same locus is read hundreds of times, and each copy is coded from
//! scratch. Measured on `ecoli_hifi` (1.55 Gbase, ~300x), that costs **0.653
//! bits/base**, against **0.068** for a codec that models reads against each
//! other — a 9.7x gap that *is* the entire lossless gap to the state of the art
//! (quality is already at parity). See `docs/design/longread.md`.
//!
//! This crate finds the overlaps that close it: which reads share a locus, at
//! what offset, in which orientation. The shape is the field's
//! (CoLoRd/NanoSpring/minimap2): **minimizers → chain colinear anchors → keep
//! the best-scoring targets**.
//!
//! ## Consensus vs read-vs-read: where the margin comes from
//!
//! Measured on this data, HiFi read error is **0.0025 edits/base** and a crude
//! edit coder costs **12.39 bits/edit**. That decomposes the design space:
//!
//! | approach | edits/base | total bits/base |
//! | --- | --- | --- |
//! | read-vs-read (code A against read B) | ~0.005 — *both* reads' errors | ~0.068 |
//! | consensus (code A against a voted consensus) | 0.0025 — one read's error | ~0.040 |
//!
//! Those are exactly CoLoRd's measured 0.0676 and our measured 0.0401. **CoLoRd
//! is read-vs-read, and the entire gap to it is the factor of two you pay for
//! coding against another erroneous read rather than a voted consensus.** Its
//! edit model is not weaker than ours; it is solving a harder problem than it
//! needs to.
//!
//! So: read-vs-read is the simpler build and lands at parity — the fallback, not
//! the goal. The consensus is the whole margin, and it is affordable only
//! because the assembly collapses (miniasm on this data: 1.01x the genome in 7
//! unitigs, so the reference costs ~0.006 bits/base — far less than the ~0.03
//! that halving the edit term saves).
//!
//! ## Why chaining, not a single anchor
//!
//! `fqxv-reorder` anchors each read on one global minimum k-mer and compares
//! ungapped within a ±8 offset window. That is sound for near-identical short
//! reads and useless here: one 15-mer survives ONT's ~10% error with
//! probability `0.9^15 ≈ 0.21`, and after the first indel an ungapped compare
//! diverges and mismatches everything downstream. Chaining inverts the odds — a
//! read carries thousands of minimizers, so hundreds survive, and a chain of
//! many weak anchors is decisive even when no single anchor is trustworthy. It
//! also absorbs indels, since a chain simply tolerates a gap between anchors.
//!
//! ## Determinism
//!
//! Output must not depend on thread count or hash iteration order (a workspace
//! invariant). The rules here: index maps are probed by key and never iterated;
//! every candidate set is sorted by a **total** order before use; parallel work
//! is split into a fixed number of chunks combined in chunk order; and chain DP
//! ties break on the smallest predecessor index.

#![forbid(unsafe_code)]

mod align;
mod chain;
mod consensus;
mod index;
mod layout;
mod minimizer;
mod overlap;
mod refine;
mod script;

pub use align::{align_banded, apply, Alignment, Op};
pub use chain::{chain, Anchor, Chain, ChainOpts};
pub use consensus::{consensus, Consensus, ConsensusOpts};
pub use index::{Index, Occ, Repeat};
pub use layout::{layout, Contig, LayoutOpts, Placement};
pub use minimizer::{minimizers, Minimizer};
pub use overlap::{find_overlaps, Overlap};
pub use refine::{place_against, Anchored};
pub use script::{chain_span, script_from_chain, ScriptOpts};

/// Errors returned by this crate.
#[derive(Debug, thiserror::Error)]
#[non_exhaustive]
pub enum Error {
    /// Read lengths do not sum to the sequence buffer size.
    #[error("read lengths ({lens}) do not match sequence bytes ({seq})")]
    LengthMismatch {
        /// Sum of the provided read lengths.
        lens: usize,
        /// Number of sequence bytes provided.
        seq: usize,
    },
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;

/// Minimizer sketching parameters.
///
/// Presets follow minimap2's, which are the field-tested operating points for
/// each error regime; `(w, k)` is recorded in the stream so decode never
/// re-derives it from data.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Sketch {
    /// Window: one minimizer is selected per `w` consecutive k-mers.
    pub w: usize,
    /// K-mer length (`1..=31`).
    pub k: usize,
}

impl Sketch {
    /// Oxford Nanopore (`w = 10, k = 15`) — minimap2's `map-ont`. Dense, because
    /// at ~10% error only a fifth of any given 15-mer's occurrences survive.
    #[must_use]
    pub const fn ont() -> Self {
        Self { w: 10, k: 15 }
    }

    /// PacBio HiFi (`w = 19, k = 19`) — minimap2's `map-hifi`. Sparser and
    /// longer: at <1% error nearly every k-mer survives, so a smaller sketch
    /// suffices and the index stays cheap.
    #[must_use]
    pub const fn hifi() -> Self {
        Self { w: 19, k: 19 }
    }

    /// Minimizers of one sequence under these parameters.
    #[must_use]
    pub fn minimizers(&self, seq: &[u8]) -> Vec<Minimizer> {
        minimizers(seq, self.w, self.k)
    }
}

impl Default for Sketch {
    /// ONT, the conservative choice: its denser sketch also works on HiFi (it
    /// only costs index size), whereas HiFi's sparse sketch misses ONT overlaps.
    fn default() -> Self {
        Self::ont()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn presets_match_minimap2_operating_points() {
        assert_eq!(Sketch::ont(), Sketch { w: 10, k: 15 });
        assert_eq!(Sketch::hifi(), Sketch { w: 19, k: 19 });
        assert_eq!(Sketch::default(), Sketch::ont());
    }

    #[test]
    fn hifi_sketch_is_sparser_than_ont() {
        let mut s = Vec::new();
        let mut x: u32 = 11;
        for _ in 0..50_000 {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            s.push(b"ACGT"[(x >> 16) as usize % 4]);
        }
        let ont = Sketch::ont().minimizers(&s).len();
        let hifi = Sketch::hifi().minimizers(&s).len();
        assert!(
            hifi < ont,
            "hifi sketch ({hifi}) must be sparser than ont ({ont})"
        );
    }
}
