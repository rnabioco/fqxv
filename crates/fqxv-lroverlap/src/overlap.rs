//! Read-to-read overlap search: index probe -> anchors -> chains -> overlaps.
//!
//! This is the crate's top-level operation. For a query read it finds the other
//! reads sharing a locus, the offset, and the orientation — everything the codec
//! needs to pick a reference read and code against it.

use crate::{
    Index,
    chain::{Anchor, Chain, ChainOpts, Chainer},
    radix::radix_sort_u128,
};

/// Pack an anchor's group-and-position key into one `u128` whose **ascending**
/// order is exactly the tuple order `(target, strand-same, tpos, qpos)` the old
/// `Vec<(u32, bool, Anchor)>` sort produced — `strand-same` false (`0`) before
/// true (`1`), matching `bool`'s `Ord`. Losslessly invertible (see [`unpack`]),
/// so the buckets need carry only the key. A read index fits 31 bits (`Occ`),
/// leaving the top bit of the 32-bit `target` field clear; `qpos`/`tpos` are full
/// `u32`.
#[inline]
fn pack(target: u32, same: bool, tpos: u32, qpos: u32) -> u128 {
    ((target as u128) << 65) | ((same as u128) << 64) | ((tpos as u128) << 32) | (qpos as u128)
}

/// Invert [`pack`]: recover `(target, same, tpos, qpos)`.
#[inline]
fn unpack(key: u128) -> (u32, bool, u32, u32) {
    let target = (key >> 65) as u32;
    let same = (key >> 64) & 1 == 1;
    let tpos = (key >> 32) as u32;
    let qpos = key as u32;
    (target, same, tpos, qpos)
}

/// A confident overlap between a query read and a target read.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Overlap {
    /// The other read.
    pub target: u32,
    /// True when the target is reverse-complemented relative to the query.
    pub strand: bool,
    /// Chaining score; higher is more confident.
    pub score: i32,
    /// Anchors supporting it.
    pub n_anchors: u32,
    /// Query span covered.
    pub q_start: u32,
    /// Query span end.
    pub q_end: u32,
    /// Target span covered, in the **query's** orientation (see [`Overlap::strand`]).
    pub t_start: u32,
    /// Target span end, in the query's orientation.
    pub t_end: u32,
}

/// Find overlaps of `query_seq` (which is read `query_read` in `idx`) against
/// every other indexed read.
///
/// Returns overlaps sorted by descending score, ties broken on `(target,
/// q_start)` so the order is total and thread-independent.
///
/// ## Reverse complement
///
/// Canonical minimizers hash a k-mer and its reverse complement to one value, so
/// a hit whose strand flag differs from the query's is a hit on the opposite
/// orientation. Its target offset is therefore in the *other* frame and must be
/// mapped into the query's before chaining, or nothing is colinear and the
/// overlap is silently lost. A k-mer at offset `p` of a read of length `len`
/// sits at `len - p - k` when that read is reverse-complemented.
#[must_use]
pub fn find_overlaps(
    idx: &Index,
    query_read: u32,
    query_seq: &[u8],
    opts: ChainOpts,
) -> Vec<Overlap> {
    let sketch = idx.sketch();
    let k = sketch.k as u32;
    let mins = sketch.seeds(query_seq);

    // The chainer's `k` MUST be the sketch's. An anchor is an exact match of
    // exactly `sketch.k` bases, and `chain` uses `k` to weight anchors, scale the
    // gap penalty, and size chain spans — so any other value mis-scores every
    // chain. `ChainOpts::default()` carries 15 (the ONT sketch), which silently
    // under-scored every HiFi chain by ~21% (k=19) and pushed marginal ones below
    // `min_score`. Overriding here rather than trusting callers: the index knows
    // the truth, and a caller cannot get it wrong if it is not theirs to pass.
    // Built once per query, not once per candidate target: it precomputes the
    // gap penalty, which pays for itself across a read's hundreds of targets and
    // would be pure overhead per target.
    let chainer = Chainer::new(ChainOpts { k, ..opts });

    // Group anchors by (target, same-orientation). Collected into a Vec and
    // sorted rather than a hash map, so iteration order is never a factor. The
    // group-and-position key packs into one `u128` whose ascending order is the
    // grouping order, so the sort is a radix sort, not a comparison sort — the
    // hot per-query sort at high coverage (see [`radix_sort_u128`]).
    let mut buckets: Vec<u128> = Vec::new();
    for m in &mins {
        // Cap the per-minimizer fan-out. A minimizer shared by every read (the
        // amplicon case: one locus, ~150k near-identical copies) otherwise
        // contributes an occurrence per read, so `buckets` grows as
        // `minimizers × reads` and OOM-kills the block (#139). A run longer than
        // the cap is subsampled with a fixed stride — deterministic, and an even
        // spread of targets rather than a biased prefix, so the reads picked
        // still represent the locus.
        let occ = idx.query(m.hash);
        let step = occ.len().div_ceil(opts.max_fanout.max(1)).max(1);
        for o in occ.iter().step_by(step) {
            if o.read() == query_read {
                continue;
            }
            let same = o.strand() == m.strand;
            let tlen = idx.read_len(o.read());
            // Map the hit into the query's frame when the strands disagree.
            let tpos = if same {
                o.pos
            } else {
                // `tlen - pos - k`; guard the arithmetic rather than trust it.
                match tlen.checked_sub(o.pos).and_then(|v| v.checked_sub(k)) {
                    Some(v) => v,
                    None => continue,
                }
            };
            buckets.push(pack(o.read(), same, tpos, m.pos));
        }
    }
    // Total order: (target, orientation, tpos, qpos) — exactly the packed key's
    // ascending order, so this radix sort matches the old comparison sort.
    radix_sort_u128(&mut buckets);

    let mut out = Vec::new();
    // Reused across every group of this query: the per-group anchor list and the
    // chainer's result buffer, cleared per group so no state leaks between them.
    let mut anchors: Vec<Anchor> = Vec::new();
    let mut chains: Vec<Chain> = Vec::new();
    let mut i = 0usize;
    while i < buckets.len() {
        // Two keys share a (target, strand) group iff their top 33 bits match.
        let group = buckets[i] >> 64;
        let (target, same, _, _) = unpack(buckets[i]);
        let mut j = i;
        anchors.clear();
        while j < buckets.len() && buckets[j] >> 64 == group {
            let (_, _, tpos, qpos) = unpack(buckets[j]);
            anchors.push(Anchor { tpos, qpos });
            j += 1;
        }
        // Anchors are already `(tpos, qpos)`-sorted (the group is a contiguous
        // run of the sorted keys), so skip the chainer's redundant re-sort.
        chainer.chain_presorted_into(&mut anchors, &mut chains);
        for c in &chains {
            out.push(Overlap {
                target,
                strand: !same,
                score: c.score,
                n_anchors: c.n_anchors,
                q_start: c.q_start,
                q_end: c.q_end,
                t_start: c.t_start,
                t_end: c.t_end,
            });
        }
        i = j;
    }

    out.sort_unstable_by(|a, b| {
        b.score
            .cmp(&a.score)
            .then(a.target.cmp(&b.target))
            .then(a.q_start.cmp(&b.q_start))
    });
    // Keep only the highest-scoring overlaps. This bounds the caller's
    // `Vec<Vec<Overlap>>` at `reads × max_overlaps`; layout uses far fewer per
    // read, so on normal blocks (where no read reaches the cap) it is a no-op.
    out.truncate(opts.max_overlaps);
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Repeat, Sketch};
    use proptest::prelude::*;

    proptest! {
        /// The packed key must invert exactly and sort exactly like the tuple it
        /// replaced — otherwise anchors regroup and the archive changes. A read
        /// index fits 31 bits, so restrict `target` accordingly.
        #[test]
        fn pack_is_invertible_and_order_preserving(
            a in (0u32..(1 << 31), any::<bool>(), any::<u32>(), any::<u32>()),
            b in (0u32..(1 << 31), any::<bool>(), any::<u32>(), any::<u32>()),
        ) {
            let (ka, kb) = (pack(a.0, a.1, a.2, a.3), pack(b.0, b.1, b.2, b.3));
            prop_assert_eq!(unpack(ka), a);
            prop_assert_eq!(unpack(kb), b);
            // The old sort key: (target, same, Anchor{tpos, qpos}).
            let ta = (a.0, a.1, Anchor { tpos: a.2, qpos: a.3 });
            let tb = (b.0, b.1, Anchor { tpos: b.2, qpos: b.3 });
            prop_assert_eq!(ka.cmp(&kb), ta.cmp(&tb));
        }
    }

    fn rand_seq(n: usize, seed: u32) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                b"ACGT"[(x >> 16) as usize % 4]
            })
            .collect()
    }

    use fqxv_dna::revcomp_acgt as revcomp;

    /// Apply substitutions and indels at roughly `rate`, like a noisy read.
    fn mutate(s: &[u8], rate: f64, seed: u32) -> Vec<u8> {
        let mut x = seed;
        let mut next = || {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            (x >> 16) as usize
        };
        let mut out = Vec::with_capacity(s.len());
        for &b in s {
            let r = (next() % 10_000) as f64 / 10_000.0;
            if r < rate {
                match next() % 3 {
                    0 => {} // deletion
                    1 => {
                        // insertion, then the base
                        out.push(b"ACGT"[next() % 4]);
                        out.push(b);
                    }
                    _ => out.push(b"ACGT"[next() % 4]), // substitution
                }
            } else {
                out.push(b);
            }
        }
        out
    }

    fn build(reads: &[Vec<u8>]) -> (Index, Vec<u32>) {
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let seq: Vec<u8> = reads.iter().flat_map(|r| r.iter().copied()).collect();
        (
            Index::build(&lens, &seq, Sketch::ont(), Repeat { drop_top_frac: 0.0 }).unwrap(),
            lens,
        )
    }

    #[test]
    fn finds_a_clean_overlap() {
        // Two reads sharing a 2 kb suffix/prefix of one genome.
        let genome = rand_seq(10_000, 5);
        let a = genome[0..4000].to_vec();
        let b = genome[2000..6000].to_vec();
        let (idx, _) = build(&[a.clone(), b.clone()]);
        let ov = find_overlaps(&idx, 0, &a, ChainOpts::default());
        assert!(!ov.is_empty(), "the overlap must be found");
        let best = ov[0];
        assert_eq!(best.target, 1);
        assert!(!best.strand, "same orientation");
        // Read a's overlap starts ~2000; read b's starts ~0.
        assert!(
            best.q_start > 1500 && best.q_start < 2500,
            "q_start {}",
            best.q_start
        );
        assert!(best.t_start < 500, "t_start {}", best.t_start);
    }

    #[test]
    fn finds_a_reverse_complement_overlap() {
        // The RC path is the one that silently loses overlaps if the coordinate
        // transform is wrong — nothing chains and the read just looks unique.
        let genome = rand_seq(10_000, 9);
        let a = genome[0..4000].to_vec();
        let b = revcomp(&genome[2000..6000]);
        let (idx, _) = build(&[a.clone(), b.clone()]);
        let ov = find_overlaps(&idx, 0, &a, ChainOpts::default());
        assert!(!ov.is_empty(), "an RC overlap must still be found");
        let best = ov[0];
        assert_eq!(best.target, 1);
        assert!(best.strand, "must be flagged reverse-complement");
        assert!(
            best.n_anchors > 10,
            "expected a solid chain, got {}",
            best.n_anchors
        );
    }

    #[test]
    fn the_hifi_sketch_is_chained_with_its_own_k() {
        // Regression: `ChainOpts::default()` carries k=15 (the ONT sketch), but
        // the HiFi sketch is k=19. Chaining HiFi anchors with k=15 under-scores
        // every chain by ~21%, pushing marginal ones below `min_score` — they
        // then look like "no overlap" rather than a scoring bug, which is
        // exactly how it hid.
        let genome = rand_seq(12_000, 41);
        let a = mutate(&genome[0..5000], 0.002, 1);
        let b = mutate(&genome[2500..7500], 0.002, 2);
        let lens = [a.len() as u32, b.len() as u32];
        let mut seq = a.clone();
        seq.extend_from_slice(&b);
        let idx = Index::build(&lens, &seq, Sketch::hifi(), Repeat { drop_top_frac: 0.0 }).unwrap();

        // Pass a deliberately WRONG k; find_overlaps must ignore it and use the
        // index's sketch, so the score is identical either way.
        let right = find_overlaps(&idx, 0, &a, ChainOpts::default());
        let wrong_k = find_overlaps(
            &idx,
            0,
            &a,
            ChainOpts {
                k: 3,
                ..ChainOpts::default()
            },
        );
        assert!(!right.is_empty(), "a HiFi overlap must be found");
        assert_eq!(
            right[0].score, wrong_k[0].score,
            "the caller's k must be ignored — the index's sketch is the truth"
        );
    }

    #[test]
    fn finds_overlap_through_hifi_error() {
        // ~0.2% error, HiFi-like. Should be trivially found.
        let genome = rand_seq(12_000, 11);
        let a = mutate(&genome[0..5000], 0.002, 1);
        let b = mutate(&genome[2500..7500], 0.002, 2);
        let (idx, _) = build(&[a.clone(), b.clone()]);
        let ov = find_overlaps(&idx, 0, &a, ChainOpts::default());
        assert!(!ov.is_empty(), "HiFi-error overlap must be found");
        assert_eq!(ov[0].target, 1);
    }

    #[test]
    fn finds_overlap_through_ont_error() {
        // ~10% error with indels — the regime that defeats a single anchor and
        // an ungapped compare. Chaining is what makes this work at all.
        let genome = rand_seq(20_000, 13);
        let a = mutate(&genome[0..8000], 0.10, 3);
        let b = mutate(&genome[4000..12_000], 0.10, 4);
        let (idx, _) = build(&[a.clone(), b.clone()]);
        let ov = find_overlaps(&idx, 0, &a, ChainOpts::default());
        assert!(!ov.is_empty(), "ONT-error overlap must be found");
        assert_eq!(ov[0].target, 1);
        assert!(
            ov[0].n_anchors >= 3,
            "expected a chain at 10% error, got {} anchors",
            ov[0].n_anchors
        );
    }

    #[test]
    fn unrelated_reads_do_not_overlap() {
        // The false-positive side: independent sequence must not chain.
        let a = rand_seq(8000, 21);
        let b = rand_seq(8000, 22);
        let (idx, _) = build(&[a.clone(), b.clone()]);
        let ov = find_overlaps(&idx, 0, &a, ChainOpts::default());
        assert!(
            ov.is_empty(),
            "unrelated reads must not overlap, got {ov:?}"
        );
    }

    #[test]
    fn a_read_never_overlaps_itself() {
        let a = rand_seq(5000, 31);
        let (idx, _) = build(std::slice::from_ref(&a));
        let ov = find_overlaps(&idx, 0, &a, ChainOpts::default());
        assert!(ov.is_empty(), "self must be excluded");
    }

    #[test]
    fn is_deterministic() {
        let genome = rand_seq(15_000, 41);
        let reads: Vec<Vec<u8>> = (0..6)
            .map(|i| mutate(&genome[i * 1500..i * 1500 + 5000], 0.05, 100 + i as u32))
            .collect();
        let (idx, _) = build(&reads);
        let a = find_overlaps(&idx, 0, &reads[0], ChainOpts::default());
        let b = find_overlaps(&idx, 0, &reads[0], ChainOpts::default());
        assert_eq!(a, b, "overlap search must be a pure function");
    }

    #[test]
    fn high_redundancy_overlaps_are_bounded() {
        // The #139 OOM shape: a large amplicon locus, many near-identical noisy
        // copies. Every read shares (almost) every minimizer, so an unbounded
        // search makes one query's anchor buffer grow as `minimizers × reads`.
        // The caps must bound both the returned overlaps and the work behind
        // them, while still finding the (abundant) real overlaps.
        let locus = rand_seq(1500, 61);
        let reads: Vec<Vec<u8>> = (0..400)
            .map(|i| mutate(&locus, 0.003, 400 + i as u32))
            .collect();
        let (idx, _) = build(&reads);

        let cap = 32;
        let opts = ChainOpts {
            max_overlaps: cap,
            ..ChainOpts::default()
        };
        let ov = find_overlaps(&idx, 0, &reads[0], opts);
        assert!(
            ov.len() <= cap,
            "the per-query overlap count must be capped, got {}",
            ov.len()
        );
        assert!(
            !ov.is_empty(),
            "the locus is shared by every read — overlaps must still be found"
        );
        // The kept overlaps are the highest-scoring: no overlap outside the cap
        // may outscore one inside it.
        let uncapped = find_overlaps(
            &idx,
            0,
            &reads[0],
            ChainOpts {
                max_overlaps: usize::MAX,
                ..ChainOpts::default()
            },
        );
        let worst_kept = ov.iter().map(|o| o.score).min().unwrap();
        let best_dropped = uncapped.iter().skip(cap).map(|o| o.score).max();
        if let Some(best_dropped) = best_dropped {
            assert!(
                worst_kept >= best_dropped,
                "truncation kept a worse overlap ({worst_kept}) than it dropped ({best_dropped})"
            );
        }
    }

    #[test]
    fn a_low_fanout_still_finds_overlaps() {
        // A fan-out far below the redundancy must still surface real overlaps:
        // the cap sheds duplicate anchors, not the signal. 40 copies of one
        // locus, clamped to 4 occurrences per minimizer — every kept anchor is a
        // genuine same-locus hit, so the query must still chain against copies.
        let locus = rand_seq(2000, 71);
        let reads: Vec<Vec<u8>> = (0..40)
            .map(|i| mutate(&locus, 0.003, 71 + i as u32))
            .collect();
        let (idx, _) = build(&reads);
        let ov = find_overlaps(
            &idx,
            0,
            &reads[0],
            ChainOpts {
                max_fanout: 4,
                ..ChainOpts::default()
            },
        );
        assert!(
            !ov.is_empty(),
            "real overlaps must survive aggressive fan-out clamping"
        );
    }

    #[test]
    fn finds_all_true_overlaps_in_a_tiling() {
        // `LEN` reads tiling a genome at `STEP`, so read i spans
        // [STEP*i, STEP*i + LEN). Read 0 = [0, 5000) therefore truly overlaps:
        //   read 1 = [2000, 7000)  -> 3000 bp
        //   read 2 = [4000, 9000)  -> 1000 bp
        // and NOT read 3 = [6000, 11000), which begins past read 0's end.
        // Recall (both, including the short 1 kb one) and precision (stops at
        // read 2) in one assertion.
        const LEN: usize = 5000;
        const STEP: usize = 2000;
        let genome = rand_seq(40_000, 51);
        let reads: Vec<Vec<u8>> = (0..8)
            .map(|i| mutate(&genome[i * STEP..i * STEP + LEN], 0.02, 200 + i as u32))
            .collect();
        let (idx, _) = build(&reads);
        let ov = find_overlaps(&idx, 0, &reads[0], ChainOpts::default());
        let mut hit: Vec<u32> = ov.iter().map(|o| o.target).collect();
        hit.sort_unstable();
        hit.dedup();
        assert_eq!(
            hit,
            vec![1, 2],
            "read 0 spans [0,{LEN}) and must overlap exactly reads 1 and 2, got {hit:?}"
        );
    }
}
