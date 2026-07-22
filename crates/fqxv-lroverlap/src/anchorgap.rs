//! Anchor-restricted (CoLoRd-style) edit-script coding.
//!
//! Both cross-read sequence coders — the shared-reference consensus codec
//! ([`crate::codec`]) and the multi-reference tiler ([`crate::tile`]) — align each
//! read against a reference window to produce its edit script. A whole-window
//! banded DP re-derives every exact match the read already shares with the
//! reference, which on low-error input (HiFi) or a good consensus is ~99% of the
//! read: pure wasted work whose cost scales with the read length rather than its
//! divergence.
//!
//! Anchor-restriction removes that waste. Seed minimizers on both sides, chain
//! the shared exact k-mers, emit each anchor span as a free [`Op::Match`]
//! copy-run, and run the banded aligner ONLY on the short inter-anchor gaps and
//! the two flanks (see [`script_from_chain`]). Alignment work then scales with
//! divergence, not with `refwin.len()`.
//!
//! Encoder-only: the produced [`Op`] stream is exactly what the DP would emit for
//! an equal-cost path, so the decoder is oblivious to which coder ran. On a read
//! whose reference shares too few exact k-mers to chain, [`anchor_chain`] returns
//! `None` and the caller aligns the whole window with the DP — the anchor path
//! never *loses* coverage, it only avoids re-proving matches the chain proved.

use std::collections::HashMap;

use crate::{
    Anchor, ChainOpts, Chainer, Op, ScriptOpts, Sketch, chain_span, codec::align_for_coding,
    script_from_chain,
};

/// Recover the collinear exact-match anchor chain between two same-orientation
/// slices.
///
/// Returns the chain as a strictly non-overlapping, ascending anchor set (ready
/// for [`script_from_chain`]), or `None` when the seeder recovers no chain (too
/// short, or too few shared exact k-mers), in which case the caller aligns the
/// whole window.
///
/// Every anchor is a genuine exact `k`-mer match: splitmix64 is injective over
/// the `2k`-bit k-mer, so equal hash + equal strand between two same-orientation
/// slices means the underlying k-mers are identical. That is
/// [`script_from_chain`]'s "anchors are trusted" precondition, which its
/// `debug_assert` also verifies.
pub(crate) fn anchor_chain(refwin: &[u8], query: &[u8], sketch: Sketch) -> Option<Vec<Anchor>> {
    let k = sketch.k as u32;
    if refwin.len() < sketch.k || query.len() < sketch.k {
        return None;
    }

    // Ref hash -> its (pos, strand) occurrences. Only looked up by key, never
    // iterated into output, so the HashMap does not affect determinism.
    let mut ref_by_hash: HashMap<u64, Vec<(u32, bool)>> = HashMap::new();
    for m in sketch.seeds(refwin) {
        ref_by_hash
            .entry(m.hash)
            .or_default()
            .push((m.pos, m.strand));
    }
    // Anchor every query seed that hits a same-strand reference seed of equal hash.
    let mut anchors: Vec<Anchor> = Vec::new();
    for m in sketch.seeds(query) {
        if let Some(hits) = ref_by_hash.get(&m.hash) {
            for &(tpos, tstrand) in hits {
                if tstrand == m.strand {
                    anchors.push(Anchor { tpos, qpos: m.pos });
                }
            }
        }
    }
    if anchors.is_empty() {
        return None;
    }

    // Chain, keep the best chain, and restrict to the anchors that lie on it.
    let chains = Chainer::new(ChainOpts {
        k,
        ..ChainOpts::default()
    })
    .chain(&mut anchors);
    let c = chains.first().copied()?;
    let mut span: Vec<Anchor> = anchors
        .iter()
        .copied()
        .filter(|a| {
            a.qpos >= c.q_start
                && a.qpos + k <= c.q_end
                && a.tpos >= c.t_start
                && a.tpos + k <= c.t_end
        })
        .collect();
    span.sort_unstable();
    span.dedup();

    // Keep a strictly non-overlapping subset (each anchor's `k`-window starts at or
    // past the previous kept anchor's window end, on BOTH sequences). Dense
    // minimizers make consecutive anchors overlap; `script_from_chain` silently
    // *skips* an anchor that a previous one already consumed past, so its true
    // coverage would end before `chain_span`'s last anchor — and the trailing flank,
    // computed from `chain_span`, would then leave a gap and desync the decoder.
    // Trimming to non-overlapping anchors makes `chain_span` the exact coverage
    // `script_from_chain` produces. Dropping a redundant anchor never loses bases:
    // the aligner recovers that span as an inter-anchor gap instead.
    let mut on_chain: Vec<Anchor> = Vec::with_capacity(span.len());
    let (mut pe_t, mut pe_q) = (0u32, 0u32);
    for a in span {
        if on_chain.is_empty() || (a.tpos >= pe_t && a.qpos >= pe_q) {
            pe_t = a.tpos + k;
            pe_q = a.qpos + k;
            on_chain.push(a);
        }
    }
    if on_chain.is_empty() {
        return None;
    }
    Some(on_chain)
}

/// Assemble the anchor-restricted edit script for `query` against `refwin` from a
/// recovered `on_chain` (see [`anchor_chain`]): emit each anchor span as a free
/// [`Op::Match`] copy-run and run the banded aligner ONLY on the short inter-anchor
/// gaps and the two flanks. Alignment work then scales with divergence rather than
/// with `refwin.len()`, which is the whole point.
///
/// # Correctness of the composition
///
/// [`script_from_chain`]'s ops are authored *relative to the chain's first
/// anchor* — applying them to `refwin[t_from..]` rebuilds `query[q_from..q_to]`.
/// The leading flank is aligned over `refwin[..t_from]` → `query[..q_from]`, which
/// consumes exactly `t_from` reference bases, so when the two op streams are
/// concatenated a linear [`crate::apply`] walk arrives at reference offset
/// `t_from` precisely as the chain ops expect. The trailing flank continues from
/// `t_to`/`q_to` the same way. The result therefore replays to `query` exactly,
/// with the same [`Op`] semantics the decoder already understands.
pub(crate) fn anchorgap_build(
    refwin: &[u8],
    query: &[u8],
    band: usize,
    on_chain: &[Anchor],
    k: u32,
) -> Vec<Op> {
    let (t_from, t_to, q_from, q_to) = chain_span(on_chain, k);
    let t_from = (t_from as usize).min(refwin.len());
    let t_to = (t_to as usize).min(refwin.len());
    let q_from = (q_from as usize).min(query.len());
    let q_to = (q_to as usize).min(query.len());

    let mut ops: Vec<Op> = Vec::new();
    // Leading flank: refwin[..t_from] -> query[..q_from].
    if t_from > 0 || q_from > 0 {
        extend_ops(
            &mut ops,
            align_for_coding(&refwin[..t_from], &query[..q_from], band).ops,
        );
    }
    // Chain interior: anchors free, only the inter-anchor gaps aligned.
    extend_ops(
        &mut ops,
        script_from_chain(refwin, query, on_chain, ScriptOpts { k, band }).0,
    );
    // Trailing flank: refwin[t_to..] -> query[q_to..].
    if t_to < refwin.len() || q_to < query.len() {
        extend_ops(
            &mut ops,
            align_for_coding(&refwin[t_to..], &query[q_to..], band).ops,
        );
    }
    ops
}

/// Append `src` onto `dst`, merging a leading op into `dst`'s trailing one when
/// they are the same kind — so the flank / chain / flank seams do not leave a
/// `Match(3), Match(4)` pair the op model would code as two symbols. Mirrors the
/// consensus/script codec's `compact`, keeping the produced stream tight.
fn extend_ops(dst: &mut Vec<Op>, src: Vec<Op>) {
    for op in src {
        match (dst.last_mut(), op) {
            (Some(Op::Match(a)), Op::Match(b)) => *a += b,
            (Some(Op::Del(a)), Op::Del(b)) => *a += b,
            (Some(Op::Ins(a)), Op::Ins(b)) => a.extend(b),
            (_, op) => dst.push(op),
        }
    }
}
