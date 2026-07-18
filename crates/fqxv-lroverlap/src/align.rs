//! Banded global alignment of two short segments.
//!
//! The codec never aligns whole reads. Chained anchors are exact k-mer matches,
//! so the only stretches needing base-level alignment are the gaps *between*
//! consecutive anchors — typically tens of bases at `w = 10`. That keeps the
//! quadratic DP over a short segment and a narrow band, which is what makes
//! aligning 14 kb reads at 18% error affordable at all.
//!
//! Clean-room: a banded Needleman-Wunsch under unit edit costs (Levenshtein),
//! with traceback compacted into a run-length edit script.

/// One step of an edit script transforming a reference segment into a query
/// segment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Op {
    /// `n` bases agree.
    Match(u32),
    /// Reference base replaced by this query base.
    Sub(u8),
    /// Query bases absent from the reference.
    Ins(Vec<u8>),
    /// `n` reference bases absent from the query.
    Del(u32),
}

/// Result of aligning a reference segment to a query segment.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Alignment {
    /// The edit script: apply to the reference to obtain the query.
    pub ops: Vec<Op>,
    /// Edit distance (substitutions + inserted + deleted bases).
    pub dist: u32,
}

/// Cell value meaning "outside the band" — never selected as a predecessor.
const INF: u32 = u32::MAX / 4;

/// Upper bound on the number of DP cells a single [`align_banded`] call may
/// allocate. Each cell costs 5 bytes (a `u32` cost and a `u8` traceback), so
/// this caps one call at ~160 MB — independent of read length — keeping a full
/// rayon pool of concurrent alignments to a few GB rather than an OOM. The band
/// is narrowed to honour it (see [`align_banded`]); narrowing only costs ratio,
/// never correctness. Reads short enough for the requested band to fit are
/// unaffected: at the default band this is every segment up to ~7.8 kb, and the
/// budget shrinks the band only on the ultra-long tail that would otherwise
/// blow up.
const MAX_DP_CELLS: usize = 32 << 20;

/// The widest band whose DP table fits [`MAX_DP_CELLS`] for a reference of
/// length `n`: the largest `b` with `(n + 1) * (2*b + 1) <= MAX_DP_CELLS`.
/// Returns at least 1 so a valid (if degenerate) band is always produced; the
/// caller checks the correctness floor separately.
fn max_band_for(n: usize) -> usize {
    let per_row = MAX_DP_CELLS / (n + 1); // widest stride (2*b + 1) that fits
    (per_row.saturating_sub(1) / 2).max(1)
}

/// Align `refr` to `query` under unit edit costs, restricted to a diagonal band
/// of half-width `band`.
///
/// Returns the script that rewrites `refr` into `query`. The band bounds the
/// work but also the answer: an indel larger than `band` cannot be represented,
/// so the alignment degrades to substitutions rather than finding the true gap.
/// Callers size `band` from the chain's observed diagonal drift.
///
/// A band is not merely an optimization here — an unbanded DP over two 14 kb
/// reads is 2×10^8 cells per pair, which at ~100 partners per read is not
/// affordable at any coverage.
#[must_use]
pub fn align_banded(refr: &[u8], query: &[u8], band: usize) -> Alignment {
    let (n, m) = (refr.len(), query.len());
    if n == 0 {
        return Alignment {
            ops: if m == 0 {
                Vec::new()
            } else {
                vec![Op::Ins(query.to_vec())]
            },
            dist: m as u32,
        };
    }
    if m == 0 {
        return Alignment {
            ops: vec![Op::Del(n as u32)],
            dist: n as u32,
        };
    }
    // The band must at least span the length difference, or no path reaches the
    // corner and the result is garbage rather than merely suboptimal.
    let need = n.abs_diff(m);
    let mut band = band.max(need) + 1;

    // MEMORY CEILING. The DP below is `(n + 1) * (2*band + 1)` cells (a u32 and a
    // u8 each), so on an ultra-long read a wide band is gigabytes in one call —
    // measured at ~6 GB for a single 286 kb ONT read at band 2048, which times a
    // rayon pool OOM-kills the whole compress. Bound the allocation to
    // [`MAX_DP_CELLS`] regardless of read length by *narrowing the band*. A
    // narrower band only ever costs ratio — any band `>= |n - m|` still reaches
    // the corner, so the alignment still round-trips exactly; it just degrades a
    // large indel it can no longer span into substitutions. Callers keep the
    // smaller of overlap-coded and order-k, so a ratio regression here can never
    // make a block worse than the within-read model.
    let max_band = max_band_for(n);
    if band > max_band {
        if need + 1 > max_band {
            // Even the correctness floor (`|n - m|`) exceeds the budget: the two
            // segments differ so much in length that no bounded band spans the
            // corner. That is a read hanging far off its reference, not a real
            // overlap. Emit a trivial rewrite — delete the reference, insert the
            // query — which round-trips exactly in O(m) memory. "Keep the
            // smaller" discards it against order-k, so ratio is unaffected.
            return Alignment {
                ops: vec![Op::Del(n as u32), Op::Ins(query.to_vec())],
                dist: (n + m) as u32,
            };
        }
        band = max_band;
    }

    // BAND-LIMITED STORAGE. Only row `i`'s window `[i-band, i+band]` is ever
    // computed, so only that slice is stored: memory is O(n * band), not
    // O(n * m). This is not a micro-optimization — callers align whole reads
    // (the consensus aligns each 14 kb read to its draft), and a full table
    // there is ~196M cells ≈ 1 GB per call, which times a rayon pool is an OOM.
    // The stride is capped at `m + 1` so a band wider than the row cannot make
    // the allocation *larger* than the full table.
    let stride = (2 * band + 1).min(m + 1);
    let lo_of = |i: usize| i.saturating_sub(band);
    // Row `i` stores column `j` at `i * stride + (j - lo_of(i))`; `j` outside the
    // row's window has no cell and reads as INF.
    let at = |i: usize, j: usize| -> Option<usize> {
        let lo = lo_of(i);
        let hi = (i + band).min(m);
        if j < lo || j > hi {
            return None;
        }
        Some(i * stride + (j - lo))
    };

    // Traceback pointers stay full — the traceback reads `from[i][j]` at every
    // step of the walk back from the corner. The score matrix does NOT: the
    // recurrence only ever reads row `i-1`, so keeping every row was ~4/5 of this
    // call's memory (a `u32` per cell) written and never re-read. Roll it to two
    // one-row windows, `prev` (row `i-1`) and `cur` (row `i`), each indexed by
    // `j - lo_of(row)` exactly as a full row's slice was. This makes the hot
    // loop's working set the `u8` `from` matrix plus two O(band) rows instead of
    // an O(n·band) `u32` table.
    let mut from = vec![0u8; (n + 1) * stride]; // 0=diag 1=up(del) 2=left(ins)
    let mut prev = vec![INF; stride];
    let mut cur = vec![INF; stride];

    // Row 0 (empty reference): dp[0][j] = j, all insertions. `lo_of(0) == 0`, so
    // the window offset of column `j` is `j`.
    if let Some(k) = at(0, 0) {
        prev[k] = 0;
    }
    for j in 1..=m.min(band) {
        if let Some(k) = at(0, j) {
            prev[k] = j as u32;
            from[k] = 2;
        }
    }
    // HOT LOOP. This is O(n * band) and, at 40x on ecoli_hifi, 63% of the entire
    // encode — so it is written against the band's structure rather than through
    // the `at`/`get` helpers, which recompute a row's bounds and build an
    // `Option` for every one of the three neighbours of every cell. Everything
    // those helpers derive is a per-ROW constant, hoisted here. The arithmetic is
    // otherwise untouched: same predecessors, same tie-break, same values.
    for i in 1..=n {
        let lo = lo_of(i);
        let hi = (i + band).min(m);
        // The previous row's window and flat base, so a neighbour is an index,
        // not a lookup.
        let lo_p = lo_of(i - 1);
        let hi_p = (i - 1 + band).min(m);
        let base = i * stride;
        let rb = refr[i - 1];

        for j in lo..=hi {
            // `j` is inside row `i`'s window by construction, so the old
            // `at(i, j)` could never fail here. `off` is its offset in `cur`.
            let off = j - lo;
            if j == 0 {
                cur[off] = i as u32;
                from[base + off] = 1;
                continue;
            }
            // A neighbour outside the previous row's window has no cell and
            // reads as INF, exactly as the full-matrix guard returned.
            let diag = if j > lo_p && j - 1 <= hi_p {
                prev[j - 1 - lo_p]
            } else {
                INF
            };
            let up = if j >= lo_p && j <= hi_p {
                prev[j - lo_p]
            } else {
                INF
            };
            // (i, j-1) is this row, one cell back — in-window iff j > lo, and it
            // was just written this row, so `cur[off - 1]` is live.
            let left = if j > lo { cur[off - 1] } else { INF };

            let cost = u32::from(rb != query[j - 1]);
            let diag = diag.saturating_add(cost);
            let up = up.saturating_add(1); // consume refr = del
            let left = left.saturating_add(1); // consume query = ins

            // Prefer diagonal, then deletion, then insertion — a fixed total
            // order, so the traceback is deterministic for equal-cost paths.
            let (best, f) = if diag <= up && diag <= left {
                (diag, 0u8)
            } else if up <= left {
                (up, 1u8)
            } else {
                (left, 2u8)
            };
            cur[off] = best;
            from[base + off] = f;
        }
        // `cur` becomes row `i`; reuse the old `prev` buffer as next row's
        // scratch. Stale cells in it are never read — every neighbour access is
        // window-guarded or a same-row cell written earlier in the sweep.
        std::mem::swap(&mut prev, &mut cur);
    }

    // Traceback from (n, m). After the final swap, `prev` holds row `n`; the
    // corner cell is at offset `m - lo_of(n)` (in-window because `band >= |n-m|`).
    let dist = at(n, m).map_or(INF, |_| prev[m - lo_of(n)]);
    let (mut i, mut j) = (n, m);
    let mut rev: Vec<Op> = Vec::new();
    while i > 0 || j > 0 {
        let f = if i == 0 {
            2
        } else if j == 0 {
            1
        } else {
            at(i, j).map_or(1, |k| from[k])
        };
        match f {
            0 => {
                if refr[i - 1] == query[j - 1] {
                    rev.push(Op::Match(1));
                } else {
                    rev.push(Op::Sub(query[j - 1]));
                }
                i -= 1;
                j -= 1;
            }
            1 => {
                rev.push(Op::Del(1));
                i -= 1;
            }
            _ => {
                rev.push(Op::Ins(vec![query[j - 1]]));
                j -= 1;
            }
        }
    }
    rev.reverse();

    // Compact runs: matches and deletions merge by count, insertions by bases.
    // A substitution is its own op — it carries a base, and runs of them are
    // rare enough not to be worth a length field.
    let mut ops: Vec<Op> = Vec::with_capacity(rev.len());
    for op in rev {
        match (ops.last_mut(), op) {
            (Some(Op::Match(a)), Op::Match(b)) => *a += b,
            (Some(Op::Del(a)), Op::Del(b)) => *a += b,
            (Some(Op::Ins(a)), Op::Ins(b)) => a.extend(b),
            (_, op) => ops.push(op),
        }
    }
    // The capacity above is the PER-BASE path length, because that is the worst
    // case and reserving it keeps the compaction realloc-free. Compaction is the
    // whole point though, so the result is ~two orders smaller — a 13 kb HiFi
    // read at 0.0025 edits/base walks ~13000 steps and compacts to ~150 ops.
    // Returning the fat buffer means every retained `Alignment` holds ~416 KB to
    // carry ~5 KB, and a caller with one per read holds 120k of them: measured at
    // 42 GB peak on ecoli_hifi at 300x, growing ~0.4 GB/s through the encode.
    // Hand back what the alignment actually is.
    ops.shrink_to_fit();
    Alignment { ops, dist }
}

/// Apply an edit script to a reference segment, producing the query segment.
///
/// The inverse of [`align_banded`], and the reason the codec can trust it: an
/// alignment is only usable if replaying it reconstructs the query exactly.
#[must_use]
pub fn apply(refr: &[u8], ops: &[Op]) -> Vec<u8> {
    let mut out = Vec::with_capacity(refr.len());
    let mut i = 0usize;
    for op in ops {
        match op {
            Op::Match(n) => {
                let n = *n as usize;
                out.extend_from_slice(&refr[i..(i + n).min(refr.len())]);
                i += n;
            }
            Op::Sub(b) => {
                out.push(*b);
                i += 1;
            }
            Op::Ins(bs) => out.extend_from_slice(bs),
            Op::Del(n) => i += *n as usize,
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    #[test]
    fn identical_segments_are_one_match_run() {
        let a = b"ACGTACGTACGT";
        let al = align_banded(a, a, 8);
        assert_eq!(al.ops, vec![Op::Match(12)]);
        assert_eq!(al.dist, 0);
    }

    #[test]
    fn a_substitution_costs_one() {
        let al = align_banded(b"ACGTACGT", b"ACGAACGT", 8);
        assert_eq!(al.dist, 1);
        assert_eq!(al.ops, vec![Op::Match(3), Op::Sub(b'A'), Op::Match(4)]);
    }

    #[test]
    fn an_insertion_is_recovered() {
        let al = align_banded(b"ACGTACGT", b"ACGTTTACGT", 8);
        assert_eq!(al.dist, 2);
        assert_eq!(apply(b"ACGTACGT", &al.ops), b"ACGTTTACGT");
    }

    #[test]
    fn a_deletion_is_recovered() {
        let al = align_banded(b"ACGTTTACGT", b"ACGTACGT", 8);
        assert_eq!(al.dist, 2);
        assert_eq!(apply(b"ACGTTTACGT", &al.ops), b"ACGTACGT");
    }

    #[test]
    fn empty_inputs() {
        assert_eq!(align_banded(b"", b"", 4).ops, vec![]);
        assert_eq!(
            align_banded(b"", b"ACGT", 4).ops,
            vec![Op::Ins(b"ACGT".to_vec())]
        );
        assert_eq!(align_banded(b"ACGT", b"", 4).ops, vec![Op::Del(4)]);
    }

    #[test]
    fn a_homopolymer_indel_is_one_del_op() {
        // ONT's signature error. It must compact into a single Del, not six —
        // the run length is what the entropy coder models cheaply.
        let al = align_banded(b"ACGAAAAAAGT", b"ACGGT", 8);
        assert_eq!(al.dist, 6);
        assert_eq!(apply(b"ACGAAAAAAGT", &al.ops), b"ACGGT");
        assert_eq!(
            al.ops.iter().filter(|o| matches!(o, Op::Del(_))).count(),
            1,
            "a homopolymer deletion must be ONE run-length op: {:?}",
            al.ops
        );
    }

    #[test]
    fn long_sequences_stay_within_band_memory() {
        // Storage is O(n * band), not O(n * m). At this size a full table would
        // be ~10^8 cells (~400 MB) per call; banded it is ~10^6. The consensus
        // aligns whole reads, so this is the real operating point, and the test
        // exists to keep the allocation honest as much as the answer.
        let mut a = Vec::new();
        let mut x: u32 = 5;
        for _ in 0..10_000 {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            a.push(b"ACGT"[(x >> 16) as usize % 4]);
        }
        // A copy with a handful of edits, well inside a narrow band.
        let mut b = a.clone();
        b[100] = b'A';
        b[5000] = b'C';
        b.remove(7000);
        let al = align_banded(&a, &b, 32);
        assert_eq!(apply(&a, &al.ops), b, "must still round-trip at length");
        assert!(al.dist <= 4, "a few edits, got {}", al.dist);
    }

    #[test]
    fn wide_band_on_a_long_read_stays_within_the_cell_budget() {
        // The OOM regression guard: a huge requested band on a long read must be
        // narrowed to honour MAX_DP_CELLS, not allocate O(n * band). The result
        // must still round-trip — narrowing costs ratio, never correctness.
        let mut a = Vec::new();
        let mut x: u32 = 9;
        for _ in 0..200_000 {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            a.push(b"ACGT"[(x >> 16) as usize % 4]);
        }
        let mut b = a.clone();
        b[1000] = b'A';
        b[150_000] = b'C';
        // A band of 4096 on a 200 kb read is ~1.6 G cells unbanded-capped; the
        // ceiling must pull it back to `max_band_for(n)`.
        let mb = max_band_for(a.len());
        assert!(mb < 4096, "budget must bind here: max_band={mb}");
        let al = align_banded(&a, &b, 4096);
        assert_eq!(apply(&a, &al.ops), b, "must round-trip at the capped band");
    }

    #[test]
    fn max_band_never_exceeds_the_cell_budget() {
        for &n in &[1usize, 1_000, 50_000, 286_296, 1_000_000] {
            let b = max_band_for(n);
            assert!(
                (n + 1) * (2 * b + 1) <= MAX_DP_CELLS,
                "n={n} band={b} exceeds budget"
            );
        }
    }

    #[test]
    fn degenerate_length_mismatch_falls_back_and_round_trips() {
        // When even the correctness floor |n-m| exceeds the budget (a long read
        // hanging off a short reference), the trivial delete-all/insert-all
        // rewrite must still reconstruct the query exactly.
        let refr = vec![b'A'; 8];
        let mut query = Vec::new();
        let mut x: u32 = 3;
        for _ in 0..40_000_000 {
            // long enough that need+1 > max_band_for(8)
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            query.push(b"ACGT"[(x >> 16) as usize % 4]);
        }
        // Guard: this really is the degenerate branch, not a normal band.
        assert!(query.len().abs_diff(refr.len()) + 1 > max_band_for(refr.len()));
        let al = align_banded(&refr, &query, 8);
        assert_eq!(apply(&refr, &al.ops), query, "fallback must round-trip");
    }

    #[test]
    fn band_is_widened_to_span_the_length_difference() {
        // A band narrower than |n-m| would leave the corner unreachable and the
        // result meaningless, so it is widened rather than trusted.
        let al = align_banded(b"ACGTACGTACGTACGT", b"ACGT", 1);
        assert_eq!(apply(b"ACGTACGTACGTACGT", &al.ops), b"ACGT");
    }

    proptest! {
        /// The property the codec depends on: replaying an alignment must
        /// reconstruct the query exactly. An alignment that does not round-trip
        /// is worse than useless — it silently corrupts the read.
        #[test]
        fn alignment_round_trips(
            refr in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 0..80),
            query in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 0..80),
        ) {
            let al = align_banded(&refr, &query, 16);
            prop_assert_eq!(apply(&refr, &al.ops), query);
        }

        /// Deterministic: equal-cost paths must resolve the same way every time.
        #[test]
        fn deterministic(
            refr in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 0..60),
            query in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 0..60),
        ) {
            prop_assert_eq!(align_banded(&refr, &query, 12), align_banded(&refr, &query, 12));
        }

        /// A wide band must find the true edit distance. Checked against a plain
        /// unbanded Levenshtein, so the DP is verified against a reference
        /// implementation rather than against itself.
        #[test]
        fn matches_unbanded_levenshtein(
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
            let al = align_banded(&refr, &query, 64);
            prop_assert_eq!(al.dist, d[n][m]);
        }
    }
}
