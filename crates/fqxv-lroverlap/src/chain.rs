//! Colinear anchor chaining (minimap2-style dynamic program).
//!
//! A shared minimizer between two reads is an *anchor*: weak evidence on its
//! own, since a 15-mer matches by chance and survives ONT error only ~21% of the
//! time. A run of anchors that advance together in both reads is strong
//! evidence, because chance matches are not colinear. Chaining turns many weak
//! anchors into one confident overlap, and absorbs indels for free: an indel
//! just shifts the diagonal, which the gap term pays for rather than being
//! derailed by (an ungapped compare, by contrast, mismatches everything after
//! the first indel).
//!
//! ## Determinism
//!
//! Anchors are sorted by a total order before the DP, and predecessor ties break
//! on the smallest index, so the chain set is a pure function of the input.

/// A shared minimizer between a query read and a target read.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Anchor {
    /// Offset in the target read.
    pub tpos: u32,
    /// Offset in the query read.
    pub qpos: u32,
}

/// A run of colinear anchors: a candidate overlap.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Chain {
    /// Chaining score (higher is better).
    pub score: i32,
    /// Number of anchors on the chain.
    pub n_anchors: u32,
    /// First and last query offsets covered.
    pub q_start: u32,
    /// Exclusive-ish end: the last anchor's query offset plus `k`.
    pub q_end: u32,
    /// First target offset covered.
    pub t_start: u32,
    /// The last anchor's target offset plus `k`.
    pub t_end: u32,
}

/// Chaining knobs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChainOpts {
    /// K-mer length: an anchor's maximum contribution to the score.
    pub k: u32,
    /// Predecessors examined per anchor. minimap2 uses ~50; the DP is O(n·h), so
    /// this bounds the cost. Chains rarely need to reach further back.
    pub lookback: usize,
    /// Largest allowed indel (diagonal shift) between consecutive anchors.
    pub max_gap: u32,
    /// Largest allowed advance between consecutive anchors on either read.
    pub max_dist: u32,
    /// Chains scoring below this are discarded.
    pub min_score: i32,
    /// Chains with fewer anchors than this are discarded.
    pub min_anchors: u32,
    /// Overlap-search memory bound, consumed by
    /// [`find_overlaps`](crate::find_overlaps) — NOT the chainer. The maximum
    /// occurrences of a single minimizer expanded into candidate anchors; a run
    /// longer than this is evenly subsampled down to it.
    ///
    /// Without it, high-redundancy input (amplicon: ~150k near-identical reads of
    /// one locus) makes every minimizer occur once per read, so a single query's
    /// anchor buffer grows as `minimizers × reads` — quadratic in the block, tens
    /// of GB, an OOM-kill (#139). The repeat filter cannot catch this: the whole
    /// frequency distribution is uniformly high, so there is no "top fraction" to
    /// drop. At the layout subsample's target coverage (~40x) a single-copy
    /// minimizer occurs far below this, so normal long-read blocks never reach the
    /// cap and are byte-identical; only the pathological case is bounded.
    pub max_fanout: usize,
    /// Overlap-search memory bound, consumed by
    /// [`find_overlaps`](crate::find_overlaps) — NOT the chainer. The maximum
    /// overlaps returned per query, keeping the highest-scoring. Bounds the
    /// materialized `Vec<Vec<Overlap>>` at `reads × this` regardless of how many
    /// reads share a locus. [`layout`](crate::layout) consumes only its top
    /// `max_candidates` (16) per read, so any value well above the subsample's
    /// per-read overlap count is lossless downstream.
    pub max_overlaps: usize,
}

impl Default for ChainOpts {
    fn default() -> Self {
        Self {
            k: 15,
            lookback: 50,
            max_gap: 5_000,
            max_dist: 5_000,
            min_score: 40,
            min_anchors: 3,
            max_fanout: 512,
            max_overlaps: 512,
        }
    }
}

/// Gap penalty for a diagonal shift of `gap` between two anchors.
///
/// Linear in the gap plus a log term, following minimap2: the linear part makes
/// a long indel progressively less attractive, and the log term keeps a single
/// large gap from being cheaper than the many small ones it would displace.
#[inline]
fn gap_cost(gap: u32, k: u32) -> i32 {
    if gap == 0 {
        return 0;
    }
    let g = gap as f64;
    let lin = 0.01 * f64::from(k) * g;
    let log = 0.5 * g.log2();
    (lin + log) as i32
}

/// A configured chainer: [`ChainOpts`] plus the gap penalty they imply,
/// precomputed.
///
/// Build one per [`find_overlaps`](crate::find_overlaps) and reuse it for every
/// candidate target — building it per target would cost far more than it saves.
///
/// ## Why the table exists
///
/// The DP evaluates [`gap_cost`] once per (anchor, predecessor) pair — about
/// 2.2e12 times on a 300x HiFi dataset, where `overlaps` is 83% of runtime — and
/// its log term was an `f64::log2`, a libm call in the innermost loop. Both of
/// that function's inputs are bounded before it is ever called: the DP skips any
/// `gap > max_gap`, and `k` is fixed for a whole `find_overlaps` (the index owns
/// it). So the penalty is not a function to evaluate, it is a table to read —
/// `max_gap + 1` entries, 20 KB at the default, small enough to stay in cache.
///
/// Every entry is exactly what the arithmetic produced, so chains — and the
/// archive — stay **byte-identical**; this is a speedup, not a heuristic change.
/// See `the_table_matches_the_arithmetic_it_replaced`.
///
/// Bundling the table WITH the opts it was built from is deliberate: a table
/// built for one `k` and used with another would silently mis-score every chain,
/// and this way that pairing is not the caller's to get wrong — the same reason
/// `find_overlaps` takes `k` from the index rather than the caller.
#[derive(Debug, Clone)]
pub struct Chainer {
    opts: ChainOpts,
    gap_cost: Vec<i32>,
}

impl Chainer {
    /// Precompute the gap penalty for `opts`.
    #[must_use]
    pub fn new(opts: ChainOpts) -> Self {
        let gap_cost = (0..=opts.max_gap).map(|g| gap_cost(g, opts.k)).collect();
        Self { opts, gap_cost }
    }

    /// The options this chainer was built from.
    #[must_use]
    pub fn opts(&self) -> ChainOpts {
        self.opts
    }

    /// The precomputed penalty for `gap`.
    #[inline]
    fn gap_penalty(&self, gap: u32) -> i32 {
        // The DP rejects `gap > max_gap` before asking, so this is always a hit.
        // Recompute rather than panic if some future caller does not.
        self.gap_cost
            .get(gap as usize)
            .copied()
            .unwrap_or_else(|| gap_cost(gap, self.opts.k))
    }

    /// Chain colinear anchors into candidate overlaps.
    ///
    /// `anchors` need not be sorted; it is sorted internally by `(tpos, qpos)`.
    /// Returns chains sorted by descending score, ties broken on `(q_start,
    /// t_start)` so the order is total and thread-independent.
    #[must_use]
    pub fn chain(&self, anchors: &mut Vec<Anchor>) -> Vec<Chain> {
        let opts = self.opts;
        anchors.sort_unstable();
        anchors.dedup();
        let n = anchors.len();
        let mut out = Vec::new();
        if n == 0 {
            return out;
        }

        // f[i]: best score of a chain ending at anchor i. p[i]: its predecessor.
        let mut f = vec![opts.k as i32; n];
        let mut p = vec![usize::MAX; n];

        for i in 0..n {
            let a = anchors[i];
            let lo = i.saturating_sub(opts.lookback);
            for j in (lo..i).rev() {
                let b = anchors[j];
                // Strictly colinear: both coordinates must advance.
                if b.tpos >= a.tpos || b.qpos >= a.qpos {
                    continue;
                }
                let dt = a.tpos - b.tpos;
                let dq = a.qpos - b.qpos;
                if dt > opts.max_dist || dq > opts.max_dist {
                    continue;
                }
                let gap = dt.abs_diff(dq);
                if gap > opts.max_gap {
                    continue;
                }
                // Overlapping anchors contribute only the new bases they cover.
                let weight = dt.min(dq).min(opts.k) as i32;
                let score = f[j] + weight - self.gap_penalty(gap);
                // Strictly-greater keeps the SMALLEST predecessor index on a tie,
                // since j descends — a total order, so the DP cannot depend on
                // iteration or thread scheduling.
                if score > f[i] {
                    f[i] = score;
                    p[i] = j;
                }
            }
        }

        // Backtrack from best-scoring anchors, each anchor used by one chain only.
        let mut order: Vec<usize> = (0..n).collect();
        // Total order: score DESC, then index ASC.
        order.sort_unstable_by(|&x, &y| f[y].cmp(&f[x]).then(x.cmp(&y)));
        let mut used = vec![false; n];

        for &end in &order {
            if used[end] || f[end] < opts.min_score {
                continue;
            }
            let mut path = Vec::new();
            let mut cur = end;
            loop {
                if used[cur] {
                    // Ran into an existing chain; stop rather than steal its anchors.
                    break;
                }
                path.push(cur);
                match p[cur] {
                    usize::MAX => break,
                    prev => cur = prev,
                }
            }
            if (path.len() as u32) < opts.min_anchors {
                continue;
            }
            for &i in &path {
                used[i] = true;
            }
            // `path` runs end -> start.
            let first = anchors[*path.last().expect("non-empty")];
            let last = anchors[path[0]];
            out.push(Chain {
                score: f[end],
                n_anchors: path.len() as u32,
                q_start: first.qpos,
                q_end: last.qpos + opts.k,
                t_start: first.tpos,
                t_end: last.tpos + opts.k,
            });
        }

        out.sort_unstable_by(|a, b| {
            b.score
                .cmp(&a.score)
                .then(a.q_start.cmp(&b.q_start))
                .then(a.t_start.cmp(&b.t_start))
        });
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Sketch;
    use proptest::prelude::*;

    fn opts() -> ChainOpts {
        ChainOpts {
            min_score: 0,
            min_anchors: 2,
            ..ChainOpts::default()
        }
    }

    /// Chain with the test options, the way every test here used to.
    fn chain(anchors: &mut Vec<Anchor>, opts: ChainOpts) -> Vec<Chain> {
        Chainer::new(opts).chain(anchors)
    }

    /// The table replaced an `f64::log2` in the innermost loop, and it is only a
    /// legitimate swap if it returns the SAME penalty for every gap the DP can
    /// reach. If it ever does not, chains re-score, overlaps move, and the
    /// archive changes — a silent ratio regression with no failing test to point
    /// at it. So check every entry against the arithmetic it replaced, for both
    /// sketch presets.
    #[test]
    fn the_table_matches_the_arithmetic_it_replaced() {
        for k in [Sketch::ont().k as u32, Sketch::hifi().k as u32] {
            let o = ChainOpts {
                k,
                ..ChainOpts::default()
            };
            let c = Chainer::new(o);
            for gap in 0..=o.max_gap {
                assert_eq!(
                    c.gap_penalty(gap),
                    gap_cost(gap, k),
                    "table disagrees with gap_cost at k={k} gap={gap}"
                );
            }
            // Past the table the DP never asks, but the fallback must still be
            // right rather than panic.
            assert_eq!(c.gap_penalty(o.max_gap + 1), gap_cost(o.max_gap + 1, k));
        }
    }

    #[test]
    fn empty_anchors_yield_no_chains() {
        assert!(chain(&mut Vec::new(), opts()).is_empty());
    }

    #[test]
    fn a_perfect_diagonal_forms_one_chain() {
        // Anchors marching in lockstep = an exact overlap at a fixed offset.
        let mut a: Vec<Anchor> = (0..20)
            .map(|i| Anchor {
                tpos: 100 + i * 20,
                qpos: i * 20,
            })
            .collect();
        let c = chain(&mut a, opts());
        assert_eq!(c.len(), 1, "one diagonal must give exactly one chain");
        assert_eq!(c[0].n_anchors, 20);
        assert_eq!(c[0].q_start, 0);
        assert_eq!(c[0].t_start, 100);
    }

    #[test]
    fn an_indel_stays_one_chain() {
        // A 5 bp shift mid-way: the whole point of chaining over ungapped
        // compare — the chain absorbs it and stays intact.
        let mut a: Vec<Anchor> = (0..10)
            .map(|i| Anchor {
                tpos: 100 + i * 20,
                qpos: i * 20,
            })
            .chain((10..20).map(|i| Anchor {
                tpos: 100 + i * 20 + 5,
                qpos: i * 20,
            }))
            .collect();
        let c = chain(&mut a, opts());
        assert_eq!(c.len(), 1, "an indel must not split the chain");
        assert_eq!(c[0].n_anchors, 20);
    }

    #[test]
    fn anticolinear_anchors_do_not_chain() {
        // Anchors going the wrong way are not an overlap.
        let mut a: Vec<Anchor> = (0..20)
            .map(|i| Anchor {
                tpos: 100 + i * 20,
                qpos: 400 - i * 20,
            })
            .collect();
        let c = chain(&mut a, opts());
        assert!(
            c.iter().all(|x| x.n_anchors < 3),
            "anti-colinear anchors must not form a long chain"
        );
    }

    #[test]
    fn random_scatter_gives_no_long_chain() {
        // Chance anchors are not colinear, so nothing substantial should chain.
        let mut x: u32 = 3;
        let mut a: Vec<Anchor> = (0..200)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                let t = (x >> 8) % 100_000;
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                let q = (x >> 8) % 100_000;
                Anchor { tpos: t, qpos: q }
            })
            .collect();
        let c = chain(
            &mut a,
            ChainOpts {
                min_anchors: 5,
                min_score: 100,
                ..ChainOpts::default()
            },
        );
        assert!(
            c.iter().all(|x| x.n_anchors < 10),
            "random anchors must not produce a long chain: {c:?}"
        );
    }

    #[test]
    fn a_big_gap_breaks_the_chain() {
        // Two diagonals separated by more than max_dist are separate overlaps.
        let mut a: Vec<Anchor> = (0..10)
            .map(|i| Anchor {
                tpos: i * 20,
                qpos: i * 20,
            })
            .chain((0..10).map(|i| Anchor {
                tpos: 50_000 + i * 20,
                qpos: 50_000 + i * 20,
            }))
            .collect();
        let c = chain(&mut a, opts());
        assert_eq!(c.len(), 2, "a gap beyond max_dist must split the chain");
    }

    #[test]
    fn chains_are_sorted_by_descending_score() {
        let mut a: Vec<Anchor> = (0..30)
            .map(|i| Anchor {
                tpos: i * 20,
                qpos: i * 20,
            })
            .chain((0..5).map(|i| Anchor {
                tpos: 90_000 + i * 20,
                qpos: 90_000 + i * 20,
            }))
            .collect();
        let c = chain(&mut a, opts());
        assert!(c.len() >= 2);
        for w in c.windows(2) {
            assert!(w[0].score >= w[1].score, "chains must descend by score");
        }
    }

    proptest! {
        /// Same anchors, same chains — regardless of the order they arrive in.
        #[test]
        fn deterministic_and_order_independent(
            mut pairs in proptest::collection::vec((0u32..5000, 0u32..5000), 0..200),
        ) {
            let mut a: Vec<Anchor> = pairs.iter().map(|&(t, q)| Anchor { tpos: t, qpos: q }).collect();
            let first = chain(&mut a.clone(), opts());
            let again = chain(&mut a, opts());
            prop_assert_eq!(&first, &again);
            // Shuffling the input must not change the result: `chain` sorts.
            pairs.reverse();
            let mut b: Vec<Anchor> = pairs.iter().map(|&(t, q)| Anchor { tpos: t, qpos: q }).collect();
            let reversed = chain(&mut b, opts());
            prop_assert_eq!(first, reversed);
        }

        /// Never panics, and every chain is internally consistent.
        #[test]
        fn chains_are_well_formed(
            pairs in proptest::collection::vec((0u32..2000, 0u32..2000), 0..150),
        ) {
            let mut a: Vec<Anchor> = pairs.iter().map(|&(t, q)| Anchor { tpos: t, qpos: q }).collect();
            for c in chain(&mut a, opts()) {
                prop_assert!(c.q_start < c.q_end);
                prop_assert!(c.t_start < c.t_end);
                prop_assert!(c.n_anchors >= 2);
            }
        }
    }
}
