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

/// Chain colinear anchors into candidate overlaps.
///
/// `anchors` need not be sorted; it is sorted internally by `(tpos, qpos)`.
/// Returns chains sorted by descending score, ties broken on `(q_start,
/// t_start)` so the order is total and thread-independent.
#[must_use]
pub fn chain(anchors: &mut Vec<Anchor>, opts: ChainOpts) -> Vec<Chain> {
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
            let score = f[j] + weight - gap_cost(gap, opts.k);
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

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn opts() -> ChainOpts {
        ChainOpts {
            min_score: 0,
            min_anchors: 2,
            ..ChainOpts::default()
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
