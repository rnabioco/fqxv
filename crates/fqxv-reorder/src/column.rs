//! Shared consensus-column primitives for the contig-assembly codecs.

use super::*;

/// Exact duplicate of the previous emitted read.
pub(crate) const OP_MATCH: u8 = 0;
/// A member of the current contig: placed at an offset on the growing consensus
/// reference, storing only the overlap's mismatches and the novel tail bases.
pub(crate) const OP_CONTIG: u8 = 1;
/// Seeds a new contig — coded as a literal via the `fqxv_seq` context model.
pub(crate) const OP_LITERAL: u8 = 2;
/// Minimum overlap with the contig to bother placing a read on it.
pub(crate) const MIN_CONTIG_OVERLAP: usize = 16;
/// Half-width of the offset search around the anchor-implied placement. The
/// shared minimizer fixes the offset exactly for substitution errors (an error
/// inside the minimizer k-mer would move the read to a different cluster), so
/// this window only has to absorb the small shifts that indels introduce. The
/// offset is stored explicitly, so widening the search is purely an
/// encoder-side choice the decoder never sees.
pub(crate) const OFF_SEARCH: i64 = 8;

/// A contig column: per-base A/C/G/T vote counts plus the current consensus
/// byte (the plurality base, or a first-seen non-ACGT byte until an ACGT wins).
#[derive(Clone)]
pub(crate) struct Column {
    pub(crate) votes: [u32; 4],
    pub(crate) base: u8,
}

/// Fold base `b` into a contig column, updating the consensus to the plurality
/// (ties resolve to the lowest A<C<G<T so encode and decode always agree).
/// Non-ACGT bytes don't vote, so a column keeps its first-seen byte until an
/// ACGT base wins — the same rule on both sides keeps the reference in sync.
#[inline]
pub(crate) fn cast_vote(col: &mut Column, b: u8) {
    let c = code_fold(b);
    if c < 4 {
        let cw = c as usize;
        col.votes[cw] += 1;
        // Only the just-incremented base can change the plurality: every other
        // count is unchanged and was already <= the old winner, so the new
        // argmax is whichever of {old winner, `cw`} has the larger key
        // `(votes, Reverse(index))` — highest count wins, ties to the lowest
        // index. The old winner's index is recovered from `col.base` (which the
        // previous vote set to the argmax); before any ACGT vote `col.base`
        // holds a non-ACGT byte and every count is 0, for which the argmax is
        // index 0 — and this first vote makes `cw` win regardless, so mapping
        // that state to 0 is exact. Byte-identical to the prior `max_by_key`.
        let ob = match col.base {
            b'A' => 0,
            b'C' => 1,
            b'G' => 2,
            b'T' => 3,
            _ => 0,
        };
        let cw_wins = col.votes[cw] > col.votes[ob] || (col.votes[cw] == col.votes[ob] && cw < ob);
        col.base = b"ACGT"[if cw_wins { cw } else { ob }];
    }
}

/// Seed a fresh contig column from the first read to cover a position.
#[inline]
pub(crate) fn seed_column(b: u8) -> Column {
    let mut col = Column {
        votes: [0; 4],
        base: b, // first-seen, until an ACGT vote takes over
    };
    cast_vote(&mut col, b);
    col
}

/// Decide whether read `cur` (with minimizer `anchor`) can be placed on the
/// current `contig` (whose seed read sat at `ref_anchor`). Returns
/// `Some((offset, overlap))` when the placement is cheaper than a literal, with
/// the mismatch positions left in `scratch` (cleared then filled), else `None`
/// (leaving `scratch` in an unspecified state). Threading `scratch` keeps the
/// per-read placement malloc-free — the caller reuses one buffer across reads.
/// Pure w.r.t. `contig`/`cur`, so both [`encode_clustered`] and [`op_stats`]
/// share one source of truth for the classification.
pub(crate) fn place_on_contig(
    contig: &[Column],
    cur: &[u8],
    anchor: u32,
    ref_anchor: u32,
    scratch: &mut Vec<usize>,
) -> Option<(usize, usize)> {
    if contig.is_empty() || cur.is_empty() {
        return None;
    }
    // The shared-minimizer anchor gives the structurally-correct offset, which is
    // exact for substitution errors (an error inside the minimizer k-mer would
    // move the read to another cluster). Try that offset first — the common path.
    // Only if it fails acceptance do we search a small window around it to rescue
    // reads an indel has shifted off the anchor. The chosen offset is stored
    // explicitly, so the search is invisible to the decoder.
    let center = ref_anchor as i64 - anchor as i64;
    // Acceptance test that only COUNTS mismatches (no allocation): the winning
    // placement's positions are materialized once, into `scratch`, at the end.
    let try_count = |off: usize| -> Option<(usize, usize)> {
        let overlap = cur.len().min(contig.len() - off);
        if overlap == 0 || overlap < MIN_CONTIG_OVERLAP.min(cur.len()) {
            return None;
        }
        // Mismatches vs the CONSENSUS of all reads placed so far.
        let mism = (0..overlap)
            .filter(|&j| cur[j] != contig[off + j].base)
            .count();
        let novel_n = cur.len() - overlap;
        // Cheap enough to be a real overlap, and smaller than a literal.
        (mism <= overlap / 4 && novel_n + mism * 2 < cur.len()).then_some((overlap, mism))
    };
    let chosen: Option<(usize, usize)> = (center >= 0 && center as usize <= contig.len())
        .then(|| try_count(center as usize).map(|(ov, _)| (center as usize, ov)))
        .flatten()
        .or_else(|| {
            // Anchor offset was rejected: scan the window for the placement with
            // the fewest mismatches (ties nearest the anchor).
            let lo = (center - OFF_SEARCH).max(0);
            let hi = (center + OFF_SEARCH).min(contig.len() as i64);
            let mut best: Option<(usize, usize)> = None;
            let mut best_key = (usize::MAX, i64::MAX);
            for off in lo..=hi {
                if off == center {
                    continue; // already tried
                }
                if let Some((ov, mism)) = try_count(off as usize) {
                    let key = (mism, (off - center).abs());
                    if key < best_key {
                        best_key = key;
                        best = Some((off as usize, ov));
                    }
                }
            }
            best
        });
    let (off, overlap) = chosen?;
    scratch.clear();
    scratch.extend((0..overlap).filter(|&j| cur[j] != contig[off + j].base));
    Some((off, overlap))
}
