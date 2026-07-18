//! Read-reordering plan: canonical-minimizer clustering into an emit order.

use super::*;
use rayon::prelude::*;
use std::borrow::Cow;

/// A reordering: the order to emit reads, and which reads are stored
/// reverse-complemented.
#[derive(Debug, Clone, Default)]
pub struct Plan {
    /// Original read indices in the order they should be emitted.
    pub order: Vec<u32>,
    /// Per original-read flag: store the read reverse-complemented.
    pub flip: Vec<bool>,
    /// Per original-read: start position of the clustering minimizer in the
    /// oriented read. Adjacent clustered reads share the minimizer, so the
    /// difference of anchors is their alignment shift (for overlap coding).
    pub anchor: Vec<u32>,
}

/// Default minimizer k-mer length.
pub const DEFAULT_K: usize = 15;

/// Minimum canonical k-mer of `read` and whether that minimizer sits on the
/// reverse strand (i.e. the read should be flipped to canonicalize it).
/// Returns `(min_canonical_kmer, flip, anchor)` where `anchor` is the start
/// position of that minimizer k-mer in the *oriented* read (the read as stored:
/// reverse-complemented iff `flip`). The anchor lets clustered reads be aligned
/// by their shared minimizer for shifted-overlap coding.
pub(crate) fn min_canonical(read: &[u8], k: usize) -> (u64, bool, u32) {
    if read.len() < k || k == 0 || k > 32 {
        return (u64::MAX, false, 0);
    }
    let mask: u64 = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let shift = 2 * (k as u64 - 1);
    let (mut fwd, mut rc, mut valid) = (0u64, 0u64, 0usize);
    let (mut best, mut best_flip, mut best_end) = (u64::MAX, false, 0usize);
    for (idx, &b) in read.iter().enumerate() {
        let c = code_fold(b);
        if c == 255 {
            fwd = 0;
            rc = 0;
            valid = 0;
            continue;
        }
        let c = u64::from(c);
        fwd = ((fwd << 2) | c) & mask;
        rc = ((rc >> 2) | ((3 - c) << shift)) & mask;
        valid += 1;
        if valid >= k {
            let (canon, is_rc) = if fwd <= rc { (fwd, false) } else { (rc, true) };
            if canon < best {
                best = canon;
                best_flip = is_rc;
                best_end = idx; // last base of the minimizing k-mer
            }
        }
    }
    // No valid (N-free) k-mer found — no minimizer, no anchor.
    if best == u64::MAX {
        return (u64::MAX, false, 0);
    }
    // Anchor = start of the minimizer k-mer in the oriented read.
    let len = read.len();
    let anchor = if best_flip {
        (len - 1 - best_end) as u32
    } else {
        (best_end + 1 - k) as u32
    };
    (best, best_flip, anchor)
}

/// Build a clustering [`Plan`] for the reads in `seq` (lengths in `lens`).
///
/// `k` is the minimizer length (see [`DEFAULT_K`]); it is clamped to `1..=32`.
#[must_use]
pub fn plan(lens: &[u32], seq: &[u8], k: usize) -> Plan {
    let k = k.clamp(1, 32);
    let n = lens.len();

    // Byte offset of each read (so the key build can run in parallel).
    let mut offs = Vec::with_capacity(n + 1);
    let mut acc = 0usize;
    for &l in lens {
        offs.push(acc);
        acc += l as usize;
    }
    offs.push(acc);

    // (canonical minimizer, oriented sequence, original index, flip, anchor).
    // Building each key is independent, so it runs across cores. The oriented
    // sequence is only a sort tiebreak, so borrow the input for the common
    // non-flipped case and allocate (via `revcomp`) only when a read flips.
    type Key<'a> = (u64, Cow<'a, [u8]>, u32, bool, u32);
    let mut keys: Vec<Key<'_>> = (0..n)
        .into_par_iter()
        .map(|i| {
            let read = &seq[offs[i]..offs[i + 1]];
            let (canon, flip, anchor) = min_canonical(read, k);
            let oriented: Cow<'_, [u8]> = if flip {
                Cow::Owned(revcomp(read))
            } else {
                Cow::Borrowed(read)
            };
            (canon, oriented, i as u32, flip, anchor)
        })
        .collect();

    // Parallel sort: cluster by minimizer, then within a cluster order by
    // anchor DESCENDING. Higher anchor = the shared minimizer sits later in the
    // read = the read starts earlier on the shared coordinate, so reads emerge
    // left-to-right and the contig assembler grows a reference rightward without
    // ever extending left. The oriented-sequence + original-index tie-breaks
    // make the comparator a TOTAL order, so the output is byte-identical
    // regardless of thread count (the determinism invariant).
    keys.par_sort_unstable_by(|a, b| {
        a.0.cmp(&b.0)
            .then_with(|| b.4.cmp(&a.4))
            .then_with(|| a.1.cmp(&b.1))
            .then_with(|| a.2.cmp(&b.2))
    });

    let order = keys.iter().map(|key| key.2).collect();
    let mut flip = vec![false; n];
    let mut anchor = vec![0u32; n];
    for key in &keys {
        flip[key.2 as usize] = key.3;
        anchor[key.2 as usize] = key.4;
    }
    Plan {
        order,
        flip,
        anchor,
    }
}

// --- clustered contig-assembly codec ----------------------------------------
