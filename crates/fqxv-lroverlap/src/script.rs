//! Edit scripts: rewrite one read as differences from another sequence.
//!
//! This is what the codec actually codes. A read that overlaps a reference (a
//! consensus, or another read) becomes a short script instead of thousands of
//! bases — the entire point of the crate.
//!
//! Anchors from `chain` are exact k-mer matches, so the spans
//! they cover are free: they emit as [`Op::Match`] with no alignment at all.
//! Only the gaps *between* consecutive anchors are aligned, and those are short.
//! That is what keeps a 14 kb read at 18% error affordable — the quadratic DP
//! only ever runs over tens of bases.

use crate::{
    align::{align_banded, Op},
    Anchor,
};

/// How to turn a chain into a script.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScriptOpts {
    /// K-mer length: how many bases an anchor vouches for.
    pub k: u32,
    /// Band half-width for inter-anchor gap alignment. Gaps are short, so this
    /// need only absorb the diagonal drift within one gap, not across the read.
    pub band: usize,
}

impl Default for ScriptOpts {
    fn default() -> Self {
        Self { k: 15, band: 32 }
    }
}

/// Rewrite `query` as edits against `refr`, guided by chained anchors.
///
/// `anchors` must be colinear and ascending in both coordinates — i.e. one
/// chain, as produced by `chain`. `refr` is the reference in the
/// same orientation as `query` (callers reverse-complement first if the overlap
/// says so).
///
/// # Anchors are trusted
///
/// **Every anchor must be a genuine exact `k`-mer match.** An anchor's span is
/// emitted as [`Op::Match`] *without comparing the bases* — that is the entire
/// performance argument, since it is what keeps the DP off the ~99% of a read
/// that anchors already cover. A false anchor therefore does not merely cost
/// ratio: it claims a match that is not there and the script silently
/// reconstructs the wrong read.
///
/// The precondition holds by construction for real input — anchors come from
/// shared minimizers, which are exact matches — and anchors are always derived
/// encoder-side from the reads, never read from a stream. A `debug_assert`
/// checks it in tests rather than trusting the argument.
///
/// The script covers `refr[t_from..t_to]` -> `query[q_from..q_to]`, where the
/// bounds are the chain's span. Bases outside the chain are the caller's
/// problem: a read hanging off the end of its reference codes that overhang
/// separately.
///
/// Returns the ops and the edit distance.
#[must_use]
pub fn script_from_chain(
    refr: &[u8],
    query: &[u8],
    anchors: &[Anchor],
    opts: ScriptOpts,
) -> (Vec<Op>, u32) {
    let mut ops: Vec<Op> = Vec::new();
    let mut dist = 0u32;
    if anchors.is_empty() {
        return (ops, 0);
    }

    // Walk the chain. `t`/`q` track how far each sequence is consumed.
    let (mut t, mut q) = (anchors[0].tpos as usize, anchors[0].qpos as usize);
    let mut first = true;

    for a in anchors {
        let (at, aq) = (a.tpos as usize, a.qpos as usize);
        // Anchors can overlap (minimizers are dense) or be redundant after a
        // previous anchor already consumed past them; skip those rather than
        // emit negative-length spans.
        if !first && (at < t || aq < q) {
            continue;
        }
        if !first {
            // The gap between the previous anchor's end and this one. Both sides
            // are short, so a banded DP is cheap; this is where every real edit
            // is found.
            let gap_t = &refr[t.min(refr.len())..at.min(refr.len())];
            let gap_q = &query[q.min(query.len())..aq.min(query.len())];
            if !gap_t.is_empty() || !gap_q.is_empty() {
                let al = align_banded(gap_t, gap_q, opts.band);
                dist += al.dist;
                ops.extend(al.ops);
            }
        }
        first = false;
        // The anchor itself: k bases that match by construction, so they are
        // emitted without alignment or comparison. See "Anchors are trusted".
        let klen = (opts.k as usize)
            .min(refr.len().saturating_sub(at))
            .min(query.len().saturating_sub(aq));
        debug_assert_eq!(
            &refr[at..at + klen],
            &query[aq..aq + klen],
            "anchor at (t={at}, q={aq}) is not an exact match — a false anchor is \
             emitted as Match and silently corrupts the read"
        );
        if klen > 0 {
            push_match(&mut ops, klen as u32);
        }
        t = at + klen;
        q = aq + klen;
    }

    // Merge adjacent match runs the anchor loop may have emitted separately.
    (compact(ops), dist)
}

fn push_match(ops: &mut Vec<Op>, n: u32) {
    match ops.last_mut() {
        Some(Op::Match(m)) => *m += n,
        _ => ops.push(Op::Match(n)),
    }
}

/// Merge adjacent same-kind ops. Dense minimizers mean consecutive anchors often
/// abut, and a run of `Match(15)` costs far more to code than one `Match(150)`.
fn compact(ops: Vec<Op>) -> Vec<Op> {
    let mut out: Vec<Op> = Vec::with_capacity(ops.len());
    for op in ops {
        match (out.last_mut(), op) {
            (Some(Op::Match(a)), Op::Match(b)) => *a += b,
            (Some(Op::Del(a)), Op::Del(b)) => *a += b,
            (Some(Op::Ins(a)), Op::Ins(b)) => a.extend(b),
            (_, op) => out.push(op),
        }
    }
    out
}

/// The chain's span: `(t_from, t_to, q_from, q_to)` — the region the script
/// covers. Empty chains yield zeros.
#[must_use]
pub fn chain_span(anchors: &[Anchor], k: u32) -> (u32, u32, u32, u32) {
    match (anchors.first(), anchors.last()) {
        (Some(f), Some(l)) => (f.tpos, l.tpos + k, f.qpos, l.qpos + k),
        _ => (0, 0, 0, 0),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::align::apply;
    use crate::{ChainOpts, Chainer, Index, Repeat, Sketch};
    use proptest::prelude::*;

    fn rand_seq(n: usize, seed: u32) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                b"ACGT"[(x >> 16) as usize % 4]
            })
            .collect()
    }

    fn mutate(s: &[u8], rate: f64, seed: u32) -> Vec<u8> {
        let mut x = seed;
        let mut next = || {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            (x >> 16) as usize
        };
        let mut out = Vec::with_capacity(s.len());
        for &b in s {
            if ((next() % 10_000) as f64 / 10_000.0) < rate {
                match next() % 3 {
                    0 => {}
                    1 => {
                        out.push(b"ACGT"[next() % 4]);
                        out.push(b);
                    }
                    _ => out.push(b"ACGT"[next() % 4]),
                }
            } else {
                out.push(b);
            }
        }
        out
    }

    /// Anchors on the exact diagonal, every `step` bases. Only valid when
    /// `query` really is `refr` on that diagonal — see "Anchors are trusted".
    fn diag(n: usize, step: usize) -> Vec<Anchor> {
        (0..n)
            .map(|i| Anchor {
                tpos: (i * step) as u32,
                qpos: (i * step) as u32,
            })
            .collect()
    }

    /// The real pipeline: index -> genuine minimizer anchors -> chain -> script.
    /// Returns the best chain's `(anchors, ops, dist)`.
    ///
    /// Tests that vary divergence MUST go through this. Synthetic anchors assert
    /// an alignment; real anchors *are* exact matches, which is the function's
    /// precondition. Handing `diag` anchors a mutated query violates it and the
    /// debug_assert fires.
    fn real_script(refr: &[u8], query: &[u8]) -> (Vec<Anchor>, Vec<Op>, u32) {
        let lens = [refr.len() as u32, query.len() as u32];
        let mut seq = refr.to_vec();
        seq.extend_from_slice(query);
        let idx = Index::build(&lens, &seq, Sketch::ont(), Repeat { drop_top_frac: 0.0 }).unwrap();
        let mut anchors: Vec<Anchor> = Vec::new();
        for m in Sketch::ont().minimizers(query) {
            for o in idx.query(m.hash) {
                if o.read == 0 && o.strand == m.strand {
                    anchors.push(Anchor {
                        tpos: o.pos,
                        qpos: m.pos,
                    });
                }
            }
        }
        let chains = Chainer::new(ChainOpts::default()).chain(&mut anchors);
        let Some(c) = chains.first().copied() else {
            return (Vec::new(), Vec::new(), 0);
        };
        let mut on_chain: Vec<Anchor> = anchors
            .iter()
            .copied()
            .filter(|a| {
                a.qpos >= c.q_start
                    && a.qpos + 15 <= c.q_end
                    && a.tpos >= c.t_start
                    && a.tpos + 15 <= c.t_end
            })
            .collect();
        on_chain.sort_unstable();
        on_chain.dedup();
        let (ops, dist) = script_from_chain(refr, query, &on_chain, ScriptOpts::default());
        (on_chain, ops, dist)
    }

    #[test]
    fn identical_sequences_give_one_match_run() {
        let s = rand_seq(500, 3);
        // 30 anchors, step 15, k=15 -> abutting, covering [0, 450).
        let a = diag(30, 15);
        let (ops, dist) = script_from_chain(&s, &s, &a, ScriptOpts::default());
        assert_eq!(dist, 0, "identical input has no edits");
        assert_eq!(
            ops,
            vec![Op::Match(450)],
            "abutting anchors must compact to ONE run: {ops:?}"
        );
    }

    #[test]
    fn empty_chain_yields_nothing() {
        let s = rand_seq(100, 5);
        let (ops, dist) = script_from_chain(&s, &s, &[], ScriptOpts::default());
        assert!(ops.is_empty());
        assert_eq!(dist, 0);
    }

    #[test]
    fn a_substitution_between_anchors_is_found() {
        let mut r = rand_seq(200, 7);
        let mut q = r.clone();
        q[100] = if r[100] == b'A' { b'C' } else { b'A' };
        // Anchors flanking the edit, not covering it.
        let a = vec![
            Anchor { tpos: 60, qpos: 60 },
            Anchor {
                tpos: 150,
                qpos: 150,
            },
        ];
        let (ops, dist) = script_from_chain(&r, &q, &a, ScriptOpts::default());
        assert_eq!(dist, 1, "one substitution: {ops:?}");
        r[100] = q[100]; // silence unused-mut
        let _ = r;
    }

    #[test]
    fn script_reconstructs_the_query_over_the_chain_span() {
        // The property the codec lives on: replaying the script must rebuild the
        // query EXACTLY over the span the chain covers. A script that does not
        // round-trip silently corrupts the read.
        //
        // Uses real anchors and a 5%-error copy WITH indels, so it exercises the
        // coordinate drift that the chainer's gap term exists to absorb.
        let r = rand_seq(2000, 11);
        let q = mutate(&r, 0.05, 12);
        let (anchors, ops, _) = real_script(&r, &q);
        assert!(!anchors.is_empty(), "a 5%-error copy must chain");

        let (t_from, t_to, q_from, q_to) = chain_span(&anchors, 15);
        let rebuilt = apply(&r[t_from as usize..(t_to as usize).min(r.len())], &ops);
        assert_eq!(
            rebuilt,
            q[q_from as usize..(q_to as usize).min(q.len())].to_vec(),
            "replaying the script must rebuild the query over the chain span"
        );
    }

    #[test]
    fn edits_scale_with_divergence() {
        // More error -> more edits. A script that does not track divergence is
        // not measuring anything.
        //
        // Must go through REAL anchors. Synthetic `diag` anchors on a mutated
        // query violate the precondition — they claim exact matches over bases
        // that were mutated — so every anchor span emits as Match without ever
        // reaching the aligner, and dist is 0 at every error rate. Real anchors
        // are exact matches by construction, so only the true gaps get aligned.
        let r = rand_seq(3000, 21);
        let clean = real_script(&r, &r).2;
        let noisy = real_script(&r, &mutate(&r, 0.02, 22)).2;
        let noisier = real_script(&r, &mutate(&r, 0.08, 23)).2;
        assert_eq!(clean, 0, "an identical copy has no edits");
        assert!(
            noisy > 0 && noisier > noisy,
            "edits must scale with divergence: {clean} < {noisy} < {noisier}"
        );
    }

    proptest! {
        /// Never panics on a redundant, duplicated or out-of-order anchor set —
        /// a chainer on noisy data produces all three.
        ///
        /// Anchors are `(p, p)` against a single sequence used as both sides, so
        /// every one is a genuine exact match and the precondition holds. Random
        /// `(t, q)` pairs over unrelated sequence would violate it, and the test
        /// would be asserting that garbage input is survivable rather than that
        /// weird-but-legal anchor ORDERING is.
        #[test]
        fn odd_anchor_sets_never_panic(
            s in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 20..200),
            raw in proptest::collection::vec(0u32..180, 0..20),
        ) {
            let n = s.len() as u32;
            let mut anchors: Vec<Anchor> = raw
                .iter()
                .filter(|&&p| p + 15 <= n)
                .map(|&p| Anchor { tpos: p, qpos: p })
                .collect();
            anchors.sort_unstable();
            let (ops, dist) = script_from_chain(&s, &s, &anchors, ScriptOpts::default());
            // Identical sequences on true anchors: no edits, whatever the order.
            prop_assert_eq!(dist, 0);
            prop_assert!(ops.iter().all(|o| matches!(o, Op::Match(_))));
        }

        /// Deterministic, through the real pipeline.
        #[test]
        fn deterministic(
            refr in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 60..200),
            seed in 0u32..1000,
        ) {
            let q = mutate(&refr, 0.05, seed);
            prop_assert_eq!(real_script(&refr, &q), real_script(&refr, &q));
        }
    }
}
