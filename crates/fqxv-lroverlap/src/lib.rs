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

// `deny` rather than `forbid`: the crate is unsafe-free except for the AVX2
// alignment backend in `align`, which opts back in at three narrowly annotated
// sites (`#[allow(unsafe_code)]`) and is proptested byte-identical to scalar.
#![deny(unsafe_code)]

mod align;
mod anchorgap;
mod chain;
mod codec;
mod consensus;
mod index;
mod layout;
mod minimizer;
mod overlap;
mod radix;
mod refine;
mod script;
mod tile;
mod wfa;

pub use align::{Alignment, Op, align_banded, apply};
pub use chain::{Anchor, Chain, ChainOpts, Chainer};
pub use codec::{
    EncodeOpts, Reference, build_reference, decode, decode_against, encode, encode_against,
};
pub use consensus::{Consensus, ConsensusOpts, consensus};
pub use index::{Index, Occ, Repeat};
pub use layout::{Contig, LayoutOpts, Placement, layout};
pub use minimizer::{Minimizer, minimizers, syncmers};
pub use overlap::{Overlap, find_overlaps};
pub use refine::{Anchored, place_against, place_all};
pub use script::{ScriptOpts, chain_span, script_from_chain};
pub use tile::{tile_decode, tile_encode};
pub use wfa::{wfa_align, wfa_align_opt, wfa_cells};

/// Ceiling on decoded bases per coded byte, used to reject a hostile `total_bases`
/// header before it drives a `vec![0u8; total_bases]` output allocation (issue
/// #142). Shared by the two long-read decoders that size their output buffer from
/// an untrusted header count — the consensus overlap codec ([`decode`] /
/// [`decode_against`]) and the tiler ([`tile_decode`]). A real block seeds every
/// distinct base at ~2 bits/base, so even a pathological all-identical-reads block
/// stays far under this (~6 K observed worst case vs the ~6 bases/byte of real
/// ONT); it only fails a crafted length. Deliberately far below `1 << 18` so the
/// bound also caps peak decode memory (`~3 × total_bases`) on a large hostile
/// input, not just `u64::MAX`.
pub(crate) const MAX_BASES_PER_BYTE: usize = 1 << 14;

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
    /// A compressed block is malformed, truncated, or from an unknown version.
    #[error("malformed lroverlap block")]
    Corrupt,
}

/// The result type for this crate.
pub type Result<T> = std::result::Result<T, Error>;

/// How read positions are selected as anchors under a [`Sketch`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeedScheme {
    /// Window minimizers: the minimum-hash k-mer of every `w`-wide window. Cheap
    /// and robust where nearly every k-mer survives (HiFi). Its weakness on noisy
    /// reads is that a base error in *any* k-mer of a window can displace an
    /// otherwise-intact k-mer, so an error-free shared k-mer is not always
    /// co-selected by two overlapping reads.
    Minimizer,
    /// Closed syncmers with `s = k - w`: a k-mer is an anchor when its minimal
    /// `s`-mer lies at its first or last position. Selection is a function of the
    /// k-mer's own bases alone, so an intact shared k-mer is co-selected
    /// regardless of neighbour errors — ~1.3–1.6× more surviving anchors at ~10%
    /// error (ONT). Density is `2/(k - s + 1) = 2/(w + 1)`, identical to the
    /// minimizer at the same `w`, so `k`, anchor count, and specificity are
    /// unchanged; only *conservation* improves. Seeding is encode-only, so this
    /// never affects an on-disk block's decodability.
    Syncmer,
}

/// Minimizer sketching parameters.
///
/// Presets follow minimap2's, which are the field-tested operating points for
/// each error regime. Seeding is encode-only — the decoder replays stored edit
/// scripts against stored consensi and never re-sketches — so `(w, k, scheme)`
/// affects ratio and speed only, never correctness, and is not written to the
/// stream.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Sketch {
    /// Window: one anchor per `w` consecutive k-mers (density `2/(w + 1)`). For
    /// [`SeedScheme::Syncmer`] this sets the s-mer length as `s = k - w`.
    pub w: usize,
    /// K-mer length (`1..=31`).
    pub k: usize,
    /// Which positions become anchors (see [`SeedScheme`]).
    pub scheme: SeedScheme,
}

impl Sketch {
    /// Oxford Nanopore (`w = 10, k = 15`) — minimap2's `map-ont` operating point,
    /// but seeded with **closed syncmers** (`s = k - w = 5`) instead of
    /// minimizers. At ~10% error a minimizer is often deselected by an error in a
    /// window *neighbour*; a syncmer survives whenever its own k-mer does, so more
    /// genuine overlaps are found at the same density and `k`. See
    /// [`SeedScheme::Syncmer`].
    #[must_use]
    pub const fn ont() -> Self {
        Self {
            w: 10,
            k: 15,
            scheme: SeedScheme::Syncmer,
        }
    }

    /// PacBio HiFi (`w = 19, k = 19`) — minimap2's `map-hifi`. Sparser and
    /// longer: at <1% error nearly every k-mer survives, so window minimizers are
    /// already near-optimal and a smaller sketch keeps the index cheap.
    #[must_use]
    pub const fn hifi() -> Self {
        Self {
            w: 19,
            k: 19,
            scheme: SeedScheme::Minimizer,
        }
    }

    /// Anchors of one sequence under these parameters — [`minimizers`] or
    /// [`syncmers`] per [`self.scheme`](SeedScheme). Both return the same
    /// [`Minimizer`] shape, so the rest of the pipeline is scheme-agnostic.
    #[must_use]
    pub fn seeds(&self, seq: &[u8]) -> Vec<Minimizer> {
        match self.scheme {
            SeedScheme::Minimizer => minimizers(seq, self.w, self.k),
            // s = k - w makes closed-syncmer density 2/(k-s+1) equal the
            // minimizer's 2/(w+1) at the same w.
            SeedScheme::Syncmer => syncmers(seq, self.k, self.k - self.w),
        }
    }

    /// Whether this sketch's platform is low-divergence (HiFi-class, <~1% error)
    /// rather than noisy (ONT-class). A sparse, long sketch (`k >= 17`) is chosen
    /// only when nearly every k-mer survives, i.e. at low error. Callers use it to
    /// pick the score-proportional WFA aligner (a big win when reads sit close to
    /// their reference) over the banded DP (which wins once divergence is high).
    #[must_use]
    pub fn is_low_divergence(&self) -> bool {
        self.k >= 17
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
        assert_eq!(
            Sketch::ont(),
            Sketch {
                w: 10,
                k: 15,
                scheme: SeedScheme::Syncmer
            }
        );
        assert_eq!(
            Sketch::hifi(),
            Sketch {
                w: 19,
                k: 19,
                scheme: SeedScheme::Minimizer
            }
        );
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
        let ont = Sketch::ont().seeds(&s).len();
        let hifi = Sketch::hifi().seeds(&s).len();
        assert!(
            hifi < ont,
            "hifi sketch ({hifi}) must be sparser than ont ({ont})"
        );
    }
}
