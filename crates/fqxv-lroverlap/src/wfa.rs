//! Edit-distance Wavefront Alignment (WFA) with traceback.
//!
//! Clean-room from Marco-Sola et al., *Fast gap-affine pairwise alignment using
//! the wavefront algorithm* (Bioinformatics 2021) — the recurrences only, never
//! any existing WFA source. This is the simplest variant: **unit edit costs**
//! (Levenshtein — match 0, substitution / insertion / deletion 1 each), the same
//! cost model as [`crate::align_banded`], so the two are drop-in comparable.
//!
//! ## Why it exists
//!
//! [`crate::align_banded`] is a banded DP: its work is `O(n · band)` no matter
//! how similar the two segments are. WFA's work instead scales with the
//! alignment **score** `s` (the number of edits). On low-divergence data — a HiFi
//! read against its consensus at ~0.5% error, a 13 kb read needing only ~60–70
//! edits — that is a different complexity class, and it degrades gracefully
//! toward the DP as divergence rises rather than hitting a band cap. This module
//! is a **standalone prototype** to measure that trade-off; it is not wired into
//! the codec's compress path, and it deliberately finds a different (equal-cost)
//! path than the DP, so its output is not byte-identical to `align_banded`.
//!
//! ## The recurrence
//!
//! A single wavefront per score `s` stores, per diagonal `k = query_pos −
//! refr_pos`, the farthest-reaching offset `M[s][k] = query_pos` (equivalently
//! `refr_pos = offset − k`). Each score step is:
//!
//! - **EXTEND**: advance along every diagonal while `refr`/`query` bases match.
//! - **TARGET**: if diagonal `k = m − n` reached offset `m`, the corner `(n, m)`
//!   is found and `s` is the edit distance.
//! - **EXPAND** to `s+1`: per diagonal `k`, take the max reach over the three
//!   unit-cost predecessors — substitution `M[s][k] + 1`, insertion
//!   `M[s][k−1] + 1`, deletion `M[s][k+1]` — then EXTEND again.
//!
//! Every wavefront is retained so the edit script can be traced back from the
//! corner. That makes traceback storage **`O(s²)`** — fine at HiFi divergence,
//! but the reason a memory-safe production form (BiWFA, `O(s)` space) would be
//! needed for the high-error ONT regime. See the crate's WFA benchmark example.

use crate::align::{Alignment, Op};

/// Sentinel offset for a diagonal not (yet) reached at a given score.
const NEG: i32 = i32::MIN;

/// Which unit-cost predecessor produced a wavefront cell — recorded implicitly
/// (recomputed in traceback) rather than stored, and mapped back to an [`Op`].
#[derive(Clone, Copy)]
enum Kind {
    /// Same diagonal, offset `+1`: a mismatch consuming one base of each.
    Sub,
    /// From diagonal `k−1`, offset `+1`: a query base absent from the reference.
    Ins,
    /// From diagonal `k+1`, offset unchanged: a reference base absent from query.
    Del,
}

/// Read `M[s][k]` from a wavefront stored with `k = 0` at index `s` (length
/// `2·s + 1`). Diagonals outside the stored range read as [`NEG`].
#[inline]
fn wf_get(wf: &[i32], s: i32, k: i32) -> i32 {
    let idx = k + s;
    if idx < 0 || idx as usize >= wf.len() {
        NEG
    } else {
        wf[idx as usize]
    }
}

/// Fold one predecessor candidate into `best`, keeping the largest offset. The
/// caller offers candidates substitution-first, then deletion, then insertion,
/// and this replaces only on a *strictly* greater offset, so equal-offset ties
/// resolve to that fixed priority — mirroring `align_banded`'s tie-break and
/// making the trace deterministic.
#[inline]
fn fold(best: &mut Option<(Kind, i32)>, kind: Kind, o: i32, k: i32, n: i32, m: i32) {
    // Validity: `query_pos = o` in `0..=m`, `refr_pos = o − k` in `0..=n`.
    if !(0..=m).contains(&o) {
        return;
    }
    let i = o - k;
    if !(0..=n).contains(&i) {
        return;
    }
    match best {
        Some((_, bo)) if *bo >= o => {}
        _ => *best = Some((kind, o)),
    }
}

/// The best predecessor of cell `(s, k)` given the previous wavefront `prev`
/// (score `s−1`): its edit [`Kind`] and the offset *before* this score's extend.
/// Returns `None` if the cell is unreachable. Used identically in the forward
/// sweep and in traceback, so the two agree by construction.
fn expand_cell(prev: &[i32], s: i32, k: i32, n: i32, m: i32) -> Option<(Kind, i32)> {
    let sp = s - 1;
    let mut best: Option<(Kind, i32)> = None;
    let sub = wf_get(prev, sp, k);
    if sub != NEG {
        fold(&mut best, Kind::Sub, sub + 1, k, n, m);
    }
    let del = wf_get(prev, sp, k + 1);
    if del != NEG {
        fold(&mut best, Kind::Del, del, k, n, m);
    }
    let ins = wf_get(prev, sp, k - 1);
    if ins != NEG {
        fold(&mut best, Kind::Ins, ins + 1, k, n, m);
    }
    best
}

/// Forward sweep: build wavefronts until the corner is reached or `max_score` is
/// exceeded. Returns every wavefront (retained for traceback) and `Some(sf)` with
/// the final score if the corner was reached, else `None` (capped). Callers pass
/// `n, m > 0`.
fn forward(
    refr: &[u8],
    query: &[u8],
    max_score: i32,
    n: i32,
    m: i32,
    k_target: i32,
) -> (Vec<Vec<i32>>, Option<i32>) {
    // Extend along diagonal `k` from offset `o` while bases match.
    let extend = |k: i32, o: i32| -> i32 {
        let mut o = o;
        let mut i = o - k;
        while i < n && o < m && refr[i as usize] == query[o as usize] {
            i += 1;
            o += 1;
        }
        o
    };

    let mut wavefronts: Vec<Vec<i32>> = Vec::new();

    // Score 0: only diagonal 0, extending through the leading match run.
    let mut wf0 = vec![NEG; 1];
    wf0[0] = extend(0, 0);
    let corner0 = k_target == 0 && wf0[0] == m;
    wavefronts.push(wf0);
    if corner0 {
        return (wavefronts, Some(0));
    }

    let mut s = 0i32;
    loop {
        s += 1;
        if s > max_score {
            return (wavefronts, None);
        }
        let mut wf = vec![NEG; (2 * s + 1) as usize];
        let klo = (-s).max(-n);
        let khi = s.min(m);
        for k in klo..=khi {
            let pre = expand_cell(&wavefronts[(s - 1) as usize], s, k, n, m);
            if let Some((_, o0)) = pre {
                wf[(k + s) as usize] = extend(k, o0);
            }
        }
        let reached = (klo..=khi).contains(&k_target) && wf[(k_target + s) as usize] == m;
        wavefronts.push(wf);
        if reached {
            return (wavefronts, Some(s));
        }
    }
}

/// Total wavefront cells retained by the forward sweep — the `O(s²)` traceback
/// footprint, in `i32` offsets (4 bytes each). A diagnostic for the memory
/// sweep in the WFA benchmark; the empty-segment fast paths store nothing.
#[must_use]
pub fn wfa_cells(refr: &[u8], query: &[u8], max_score: u32) -> usize {
    let (n_us, m_us) = (refr.len(), query.len());
    if n_us == 0 || m_us == 0 {
        return 0;
    }
    let max_score = i32::try_from(max_score).unwrap_or(i32::MAX);
    let (n, m) = (n_us as i32, m_us as i32);
    let (wavefronts, _) = forward(refr, query, max_score, n, m, m - n);
    wavefronts.iter().map(Vec::len).sum()
}

/// Align `refr` to `query` under unit edit costs with the wavefront algorithm,
/// returning the edit script that rewrites `refr` into `query` (so
/// `apply(refr, ops) == query`).
///
/// `max_score` caps the edit distance explored, exactly as `align_banded`'s band
/// bounds its work: if the true distance would exceed it, this returns the same
/// trivial fallback `align_banded` uses when no bounded band can span the corner
/// — `vec![Op::Del(n), Op::Ins(query)]` — which round-trips in `O(m)` memory and
/// keeps a pathological pair from exploding the `O(s²)` wavefront storage.
///
/// The path found is optimal (WFA finds the true edit distance when uncapped) but
/// generally *different* from `align_banded`'s equal-cost path, so output is not
/// byte-identical to the DP.
#[must_use]
pub fn wfa_align(refr: &[u8], query: &[u8], max_score: u32) -> Alignment {
    let (n_us, m_us) = (refr.len(), query.len());
    // Empty-segment fast paths, identical to `align_banded`.
    if n_us == 0 {
        return Alignment {
            ops: if m_us == 0 {
                Vec::new()
            } else {
                vec![Op::Ins(query.to_vec())]
            },
            dist: m_us as u32,
        };
    }
    if m_us == 0 {
        return Alignment {
            ops: vec![Op::Del(n_us as u32)],
            dist: n_us as u32,
        };
    }

    let (n, m) = (n_us as i32, m_us as i32);
    let k_target = m - n;
    let max_score = i32::try_from(max_score).unwrap_or(i32::MAX);

    let (wavefronts, sf) = forward(refr, query, max_score, n, m, k_target);
    let Some(sf) = sf else {
        // Capped: the same bounded, round-tripping rewrite the DP falls back to.
        return Alignment {
            ops: vec![Op::Del(n_us as u32), Op::Ins(query.to_vec())],
            dist: (n_us + m_us) as u32,
        };
    };

    // Traceback from the corner. At each score the cell is the edit that reached
    // it followed by an extend run; walking down to score 0 and reversing yields
    // the forward script. Emitting `Match(run)`, `Sub`, `Ins(one)`, `Del(1)` and
    // then run-length compacting reproduces `align_banded`'s output shape.
    let mut rev: Vec<Op> = Vec::new();
    let mut k = k_target;
    let mut s = sf;
    while s > 0 {
        let final_o = wavefronts[s as usize][(k + s) as usize];
        let (kind, pre) =
            expand_cell(&wavefronts[(s - 1) as usize], s, k, n, m).expect("traced cell reachable");
        let matches = final_o - pre; // extend run AFTER this edit (forward-last)
        if matches > 0 {
            rev.push(Op::Match(matches as u32));
        }
        match kind {
            Kind::Sub => rev.push(Op::Sub(query[(pre - 1) as usize])),
            Kind::Ins => {
                rev.push(Op::Ins(vec![query[(pre - 1) as usize]]));
                k -= 1;
            }
            Kind::Del => {
                rev.push(Op::Del(1));
                k += 1;
            }
        }
        s -= 1;
    }
    // Score 0: the leading match run from (0, 0).
    let lead = wavefronts[0][0];
    if lead > 0 {
        rev.push(Op::Match(lead as u32));
    }
    rev.reverse();

    // Compact runs exactly as `align_banded` does.
    let mut ops: Vec<Op> = Vec::with_capacity(rev.len());
    for op in rev {
        match (ops.last_mut(), op) {
            (Some(Op::Match(a)), Op::Match(b)) => *a += b,
            (Some(Op::Del(a)), Op::Del(b)) => *a += b,
            (Some(Op::Ins(a)), Op::Ins(b)) => a.extend(b),
            (_, op) => ops.push(op),
        }
    }
    ops.shrink_to_fit();
    Alignment {
        ops,
        dist: sf as u32,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::{align_banded, apply};
    use proptest::prelude::*;

    // A big cap: past any edit distance the small proptest inputs can reach, so
    // the wavefront is never spuriously capped when correctness is under test.
    const BIG: u32 = 1 << 20;

    #[test]
    fn identical_segments_are_one_match_run() {
        let a = b"ACGTACGTACGT";
        let al = wfa_align(a, a, BIG);
        assert_eq!(al.ops, vec![Op::Match(12)]);
        assert_eq!(al.dist, 0);
    }

    #[test]
    fn a_substitution_costs_one() {
        let al = wfa_align(b"ACGTACGT", b"ACGAACGT", BIG);
        assert_eq!(al.dist, 1);
        assert_eq!(apply(b"ACGTACGT", &al.ops), b"ACGAACGT");
    }

    #[test]
    fn an_insertion_is_recovered() {
        let al = wfa_align(b"ACGTACGT", b"ACGTTTACGT", BIG);
        assert_eq!(al.dist, 2);
        assert_eq!(apply(b"ACGTACGT", &al.ops), b"ACGTTTACGT");
    }

    #[test]
    fn a_deletion_is_recovered() {
        let al = wfa_align(b"ACGTTTACGT", b"ACGTACGT", BIG);
        assert_eq!(al.dist, 2);
        assert_eq!(apply(b"ACGTTTACGT", &al.ops), b"ACGTACGT");
    }

    #[test]
    fn empty_inputs() {
        assert_eq!(wfa_align(b"", b"", BIG).ops, vec![]);
        assert_eq!(
            wfa_align(b"", b"ACGT", BIG).ops,
            vec![Op::Ins(b"ACGT".to_vec())]
        );
        assert_eq!(wfa_align(b"ACGT", b"", BIG).ops, vec![Op::Del(4)]);
    }

    #[test]
    fn a_homopolymer_indel_compacts_to_one_del() {
        let refr = b"ACGAAAAAAGT";
        let al = wfa_align(refr, b"ACGGT", BIG);
        assert_eq!(al.dist, 6);
        assert_eq!(apply(refr, &al.ops), b"ACGGT");
        assert_eq!(
            al.ops.iter().filter(|o| matches!(o, Op::Del(_))).count(),
            1,
            "a homopolymer deletion must be ONE run-length op: {:?}",
            al.ops
        );
    }

    #[test]
    fn cap_returns_the_bounded_fallback() {
        // A maximally divergent pair (edit distance 40) with a tiny cap must fall
        // back to the delete-all/insert-all rewrite rather than explode — and the
        // fallback must still round-trip.
        let refr = vec![b'A'; 40];
        let query = vec![b'C'; 40];
        let al = wfa_align(&refr, &query, 5);
        assert_eq!(al.ops, vec![Op::Del(40), Op::Ins(query.clone())]);
        assert_eq!(al.dist, 80);
        assert_eq!(apply(&refr, &al.ops), query);
    }

    #[test]
    fn cells_grow_with_score_and_stay_bounded_by_the_cap() {
        // Identical sequences reach the corner at score 0 (one wavefront cell);
        // a divergent capped pair stores O(cap²) and no more.
        assert_eq!(wfa_cells(b"ACGTACGT", b"ACGTACGT", BIG), 1);
        let refr = vec![b'A'; 200];
        let query = vec![b'C'; 200];
        let cap = 8u32;
        let cells = wfa_cells(&refr, &query, cap);
        // Wavefronts 0..=cap+1 each length 2s+1: bounded, never O(n·m).
        let bound: usize = (0..=cap as usize + 1).map(|s| 2 * s + 1).sum();
        assert!(cells <= bound, "cells {cells} exceed cap bound {bound}");
    }

    proptest! {
        /// Non-negotiable: replaying a WFA alignment must reconstruct the query.
        #[test]
        fn wfa_round_trips(
            refr in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 0..80),
            query in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 0..80),
        ) {
            let al = wfa_align(&refr, &query, BIG);
            prop_assert_eq!(apply(&refr, &al.ops), query);
        }

        /// Optimality: uncapped WFA finds the true edit distance, which a wide
        /// exact band DP also finds — so their distances must agree.
        #[test]
        fn wfa_matches_wide_band_dp(
            refr in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 0..80),
            query in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 0..80),
        ) {
            let w = wfa_align(&refr, &query, BIG);
            let d = align_banded(&refr, &query, 128);
            prop_assert_eq!(w.dist, d.dist);
        }

        /// Optimality against an independent brute-force Levenshtein reference,
        /// so WFA is verified against a third implementation, not just the DP.
        #[test]
        fn wfa_matches_brute_levenshtein(
            refr in proptest::collection::vec(proptest::sample::select(&b"ACG"[..]), 0..40),
            query in proptest::collection::vec(proptest::sample::select(&b"ACG"[..]), 0..40),
        ) {
            let (n, m) = (refr.len(), query.len());
            let mut d = vec![vec![0u32; m + 1]; n + 1];
            for i in 0..=n { d[i][0] = i as u32; }
            for j in 0..=m { d[0][j] = j as u32; }
            for i in 1..=n {
                for j in 1..=m {
                    let c = u32::from(refr[i - 1] != query[j - 1]);
                    d[i][j] = (d[i - 1][j - 1] + c).min(d[i - 1][j] + 1).min(d[i][j - 1] + 1);
                }
            }
            prop_assert_eq!(wfa_align(&refr, &query, BIG).dist, d[n][m]);
        }

        /// Deterministic: identical inputs yield an identical alignment.
        #[test]
        fn wfa_deterministic(
            refr in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 0..60),
            query in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 0..60),
        ) {
            prop_assert_eq!(wfa_align(&refr, &query, BIG), wfa_align(&refr, &query, BIG));
        }
    }
}
