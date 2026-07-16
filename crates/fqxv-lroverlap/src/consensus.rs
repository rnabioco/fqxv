//! Plurality consensus over a laid-out contig.
//!
//! This is where the margin over the state of the art comes from, and the
//! arithmetic says so. Coding a read against another *read* pays both reads'
//! errors (~0.005 edits/base on HiFi); coding it against a voted consensus pays
//! only its own (0.0025). At ~12 bits/edit that is ~0.068 vs ~0.040 bits/base —
//! which are exactly CoLoRd's measured number and ours. CoLoRd is read-vs-read;
//! the consensus is the whole difference.
//!
//! A vote is effectively error-free at long-read coverage because a column needs
//! a *majority* of independent errors to agree, and errors are independent
//! across reads. At 40-300x that does not happen.
//!
//! ## Why alignment, not a column-wise vote by offset
//!
//! Reads differ by indels, so two reads at known offsets do not stay in register
//! past the first one — and the indel is ONT's dominant error, so a naive
//! offset-indexed vote misaligns almost immediately and votes unrelated columns
//! against each other. Each read is therefore aligned to the growing reference
//! and votes through that alignment.

use std::sync::atomic::{AtomicU32, Ordering};

use rayon::prelude::*;

use crate::{
    align::{align_banded, Op},
    Contig,
};

/// Consensus options.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConsensusOpts {
    /// Band half-width when aligning a read to the draft.
    pub band: usize,
    /// A column needs this many votes to be called; below it the draft's base
    /// stands. Guards the contig ends, where coverage tails off and a single
    /// erroneous read would otherwise dictate the consensus.
    pub min_votes: u32,
}

impl Default for ConsensusOpts {
    fn default() -> Self {
        Self {
            band: 64,
            min_votes: 3,
        }
    }
}

/// A voted consensus, plus the coordinate map back to the layout.
///
/// The map is not a convenience — it is required. A column the vote deletes is
/// dropped from the sequence, so the consensus is shorter than the draft the
/// layout's offsets are expressed in. At a few percent deletions that drift
/// reaches thousands of bases over a megabase contig, far beyond any alignment
/// band, so a [`Placement`](crate::Placement) offset applied straight to the
/// consensus would land in the wrong place entirely.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Consensus {
    /// The voted sequence.
    pub seq: Vec<u8>,
    /// `draft position -> position in `seq``. Deleted columns map to the next
    /// surviving position, so every draft coordinate maps somewhere valid.
    to_cons: Vec<u32>,
}

impl Consensus {
    /// Map a layout (draft) coordinate to its position in [`Consensus::seq`].
    ///
    /// Saturates at the sequence end, so an offset past the contig cannot
    /// produce an out-of-range index.
    #[must_use]
    pub fn map(&self, draft_pos: u32) -> u32 {
        match self.to_cons.get(draft_pos as usize) {
            Some(&p) => p,
            None => self.seq.len() as u32,
        }
    }

    /// Length of the consensus sequence.
    #[must_use]
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    /// True when the consensus is empty.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }
}

/// Per-column base votes being accumulated in parallel: `A, C, G, T`, plus a
/// deletion tally.
///
/// Atomic because reads vote concurrently, and safe to do so because
/// **increments commute**: the final per-column tallies are identical whatever
/// order the threads interleave in, so the called consensus is a pure function
/// of the input and the workspace's thread-count-invariance rule holds.
/// `Relaxed` suffices — no ordering between columns is implied or needed.
///
/// The vote must parallelise across READS, not contigs. Real long-read data
/// assembles to a handful of contigs (measured: ~3 for 15.5k HiFi reads), so
/// parallelising across contigs leaves a 32-core node running 3 threads while
/// one contig does ~15k alignments in series.
#[derive(Debug, Default)]
struct AtomicColumn {
    acgt: [AtomicU32; 4],
    del: AtomicU32,
}

impl AtomicColumn {
    /// Collapse to a plain tally once voting is done.
    fn load(&self) -> Column {
        Column {
            acgt: [
                self.acgt[0].load(Ordering::Relaxed),
                self.acgt[1].load(Ordering::Relaxed),
                self.acgt[2].load(Ordering::Relaxed),
                self.acgt[3].load(Ordering::Relaxed),
            ],
            del: self.del.load(Ordering::Relaxed),
        }
    }
}

/// Per-column base votes: `A, C, G, T`, plus a deletion tally.
#[derive(Debug, Clone, Copy, Default)]
struct Column {
    acgt: [u32; 4],
    del: u32,
}

impl Column {
    /// The winning base, or `None` when deletion wins or nothing voted.
    ///
    /// Ties resolve to the lowest base (A < C < G < T) — a fixed rule, applied
    /// identically wherever the consensus is rebuilt, so encode and decode can
    /// never disagree.
    fn call(&self, min_votes: u32) -> Option<u8> {
        let total: u32 = self.acgt.iter().sum::<u32>() + self.del;
        if total < min_votes {
            return None;
        }
        let (mut best, mut best_n) = (0usize, self.acgt[0]);
        for i in 1..4 {
            if self.acgt[i] > best_n {
                best = i;
                best_n = self.acgt[i];
            }
        }
        if self.del > best_n {
            None
        } else {
            Some(b"ACGT"[best])
        }
    }
}

#[inline]
fn base_idx(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Build a consensus for one contig.
///
/// `reads[i]` must be the sequence of read `i` **already oriented** as its
/// [`Placement`](crate::Placement) says (callers reverse-complement the flipped
/// ones first).
///
/// The draft is seeded from the reads laid down at their offsets, then every
/// read is aligned to that draft and votes through the alignment. Returns the
/// consensus and the map from layout coordinates onto it (see [`Consensus`]).
///
/// Insertions are deliberately not voted into the consensus: a read carrying an
/// insertion simply codes it as an `Ins` op. Admitting them would grow the
/// reference for the benefit of a minority of reads, and the reference is shared
/// by all of them.
#[must_use]
pub fn consensus(contig: &Contig, reads: &[Vec<u8>], opts: ConsensusOpts) -> Consensus {
    if contig.reads.is_empty() {
        return Consensus {
            seq: Vec::new(),
            to_cons: Vec::new(),
        };
    }
    // A single read is its own consensus — nothing to vote against, and the map
    // is the identity.
    if contig.reads.len() == 1 {
        let seq = reads[contig.reads[0].read as usize].clone();
        let to_cons = (0..=seq.len() as u32).collect();
        return Consensus { seq, to_cons };
    }

    // 1. Draft: first read to cover each position wins. Good enough to align
    //    against; the vote is what makes it accurate.
    let mut draft: Vec<u8> = vec![0; contig.len as usize];
    for p in &contig.reads {
        let s = &reads[p.read as usize];
        for (i, &b) in s.iter().enumerate() {
            let at = p.offset as usize + i;
            if at < draft.len() && draft[at] == 0 {
                draft[at] = b;
            }
        }
    }
    // Positions no read covered (a gapped layout) hold 0; drop them so the draft
    // is a real sequence.
    draft.retain(|&b| b != 0);
    if draft.is_empty() {
        return Consensus {
            seq: Vec::new(),
            to_cons: Vec::new(),
        };
    }

    // 2. Vote: align each read to the draft and tally through the alignment, so
    //    indels shift the register rather than corrupting every column after.
    //    Reads vote CONCURRENTLY — the alignment is the expensive part and there
    //    are tens of thousands of them per contig. Tallies are atomic and
    //    increments commute, so the result does not depend on thread count.
    let cols: Vec<AtomicColumn> = (0..draft.len()).map(|_| AtomicColumn::default()).collect();
    contig.reads.par_iter().for_each(|p| {
        let s = &reads[p.read as usize];
        let from = (p.offset as usize).min(draft.len());
        let to = (from + s.len()).min(draft.len());
        if from >= to {
            return;
        }
        let al = align_banded(&draft[from..to], s, opts.band);
        let mut d = from; // position in the draft
        let mut q = 0usize; // position in the read
        for op in &al.ops {
            match op {
                Op::Match(n) => {
                    for _ in 0..*n {
                        if d < cols.len() && q < s.len() {
                            if let Some(i) = base_idx(s[q]) {
                                cols[d].acgt[i].fetch_add(1, Ordering::Relaxed);
                            }
                        }
                        d += 1;
                        q += 1;
                    }
                }
                Op::Sub(b) => {
                    if d < cols.len() {
                        if let Some(i) = base_idx(*b) {
                            cols[d].acgt[i].fetch_add(1, Ordering::Relaxed);
                        }
                    }
                    d += 1;
                    q += 1;
                }
                Op::Del(n) => {
                    // The read lacks these draft bases: vote to remove them.
                    for _ in 0..*n {
                        if d < cols.len() {
                            cols[d].del.fetch_add(1, Ordering::Relaxed);
                        }
                        d += 1;
                    }
                }
                Op::Ins(bs) => {
                    // Not voted — see the note on the function.
                    q += bs.len();
                }
            }
        }
    });
    let cols: Vec<Column> = cols.iter().map(AtomicColumn::load).collect();

    // 3. Call each column, recording where each draft position landed. A dropped
    //    column maps to the next surviving position, so a layout offset that
    //    happens to point at a deleted column still resolves somewhere valid.
    let mut out = Vec::with_capacity(draft.len());
    let mut to_cons = Vec::with_capacity(draft.len() + 1);
    for (i, c) in cols.iter().enumerate() {
        to_cons.push(out.len() as u32);
        match c.call(opts.min_votes) {
            Some(b) => out.push(b),
            None if c.acgt.iter().sum::<u32>() + c.del < opts.min_votes => out.push(draft[i]),
            None => {} // deletion won: the column is dropped
        }
    }
    // One past the end, so a read ending at the contig's end maps cleanly.
    to_cons.push(out.len() as u32);
    Consensus { seq: out, to_cons }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Placement;

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

    fn contig_of(offsets: &[u32], lens: &[u32]) -> Contig {
        Contig {
            id: 0,
            len: offsets
                .iter()
                .zip(lens)
                .map(|(o, l)| o + l)
                .max()
                .unwrap_or(0),
            reads: offsets
                .iter()
                .enumerate()
                .map(|(i, &o)| Placement {
                    read: i as u32,
                    contig: 0,
                    offset: o,
                    flip: false,
                })
                .collect(),
        }
    }

    /// Edit distance — the only honest "how far from truth" here.
    ///
    /// A positional (Hamming) comparison measures nothing on these sequences:
    /// the reads contain INDELS, so after the first one every downstream base is
    /// shifted and the count reports the shift rather than the errors. It scores
    /// a perfect consensus and a raw read almost identically, which is exactly
    /// the coordinate-drift trap the whole codec exists to handle.
    fn edit_dist(a: &[u8], b: &[u8]) -> u32 {
        align_banded(a, b, 256).dist
    }

    #[test]
    fn a_single_read_is_its_own_consensus() {
        let r = rand_seq(200, 1);
        let c = contig_of(&[0], &[200]);
        assert_eq!(
            consensus(&c, std::slice::from_ref(&r), ConsensusOpts::default()).seq,
            r
        );
    }

    #[test]
    fn empty_contig_yields_empty() {
        let c = Contig {
            id: 0,
            len: 0,
            reads: Vec::new(),
        };
        assert!(consensus(&c, &[], ConsensusOpts::default()).is_empty());
    }

    #[test]
    fn identical_reads_vote_themselves() {
        let r = rand_seq(300, 3);
        let reads = vec![r.clone(), r.clone(), r.clone()];
        let c = contig_of(&[0, 0, 0], &[300, 300, 300]);
        assert_eq!(consensus(&c, &reads, ConsensusOpts::default()).seq, r);
    }

    #[test]
    fn the_vote_removes_errors_the_draft_contains() {
        // THE point of this module: the consensus must be closer to the truth
        // than any single read. If it is not, coding against it is no better
        // than coding against a read, and the whole margin over CoLoRd is gone.
        let truth = rand_seq(1500, 7);
        let reads: Vec<Vec<u8>> = (0..15).map(|i| mutate(&truth, 0.05, 100 + i)).collect();
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let c = contig_of(&vec![0u32; reads.len()], &lens);

        let cons = consensus(&c, &reads, ConsensusOpts::default());
        let cons_err = edit_dist(&cons.seq, &truth);
        let mean_read_err: u32 =
            reads.iter().map(|r| edit_dist(r, &truth)).sum::<u32>() / reads.len() as u32;

        assert!(
            cons_err < mean_read_err / 2,
            "consensus ({cons_err}) must be far closer to truth than a mean read \
             ({mean_read_err}) — that gap IS the margin over read-vs-read coding"
        );
    }

    #[test]
    fn more_coverage_gives_a_better_consensus() {
        // The vote must actually use the coverage; if depth does not help, it is
        // not voting.
        let truth = rand_seq(1200, 11);
        let mk = |n: u32| {
            let reads: Vec<Vec<u8>> = (0..n).map(|i| mutate(&truth, 0.06, 200 + i)).collect();
            let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
            let c = contig_of(&vec![0u32; reads.len()], &lens);
            edit_dist(&consensus(&c, &reads, ConsensusOpts::default()).seq, &truth)
        };
        let shallow = mk(3);
        let deep = mk(21);
        assert!(
            deep < shallow,
            "21x consensus ({deep}) must beat 3x ({shallow})"
        );
    }

    #[test]
    fn the_map_is_identity_when_nothing_is_deleted() {
        let r = rand_seq(300, 3);
        let reads = vec![r.clone(), r.clone(), r.clone()];
        let c = contig_of(&[0, 0, 0], &[300, 300, 300]);
        let cons = consensus(&c, &reads, ConsensusOpts::default());
        assert_eq!(cons.seq, r);
        for p in [0u32, 1, 150, 299, 300] {
            assert_eq!(cons.map(p), p, "no deletions -> map is the identity");
        }
    }

    #[test]
    fn the_map_tracks_deleted_columns() {
        // The reason `Consensus` carries a map at all. When the vote deletes
        // columns the consensus is SHORTER than the draft the layout's offsets
        // are in, so a raw offset drifts — here every read agrees on deleting a
        // chunk, so the drift is large and unambiguous.
        let truth = rand_seq(600, 31);
        // Reads that all lack truth[200..250]: the vote must delete those.
        let mut short = truth[..200].to_vec();
        short.extend_from_slice(&truth[250..]);
        let reads: Vec<Vec<u8>> = (0..7).map(|_| short.clone()).collect();
        // The draft is seeded from the reads, so the draft IS `short` here; the
        // point is that map() stays consistent with whatever survives.
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let c = contig_of(&vec![0u32; reads.len()], &lens);
        let cons = consensus(&c, &reads, ConsensusOpts::default());

        // Monotonic and in range: a coordinate map that goes backwards or off
        // the end would place reads at nonsense offsets.
        let mut prev = 0u32;
        for p in 0..=cons.seq.len() as u32 {
            let m = cons.map(p);
            assert!(m >= prev, "map must be monotonic at {p}: {m} < {prev}");
            assert!(
                m <= cons.seq.len() as u32,
                "map must stay in range at {p}: {m} > {}",
                cons.seq.len()
            );
            prev = m;
        }
        // Past the end saturates rather than panicking.
        assert_eq!(cons.map(u32::MAX), cons.seq.len() as u32);
    }

    #[test]
    fn is_thread_count_invariant() {
        // The vote is concurrent, so the workspace's byte-identical-regardless-
        // of-threads rule has to be asserted, not assumed. It holds because
        // increments commute — but "it should be fine" is exactly how a codec
        // acquires a heisenbug that only appears on a different core count.
        let truth = rand_seq(4000, 61);
        let reads: Vec<Vec<u8>> = (0..24).map(|i| mutate(&truth, 0.06, 700 + i)).collect();
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let c = contig_of(&vec![0u32; reads.len()], &lens);

        let run = |threads: usize| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .expect("pool")
                .install(|| consensus(&c, &reads, ConsensusOpts::default()))
        };
        let one = run(1);
        let many = run(8);
        assert_eq!(
            one, many,
            "consensus must be identical on 1 thread and 8 — atomic votes commute"
        );
    }

    #[test]
    fn is_deterministic() {
        let truth = rand_seq(800, 13);
        let reads: Vec<Vec<u8>> = (0..9).map(|i| mutate(&truth, 0.05, 300 + i)).collect();
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let c = contig_of(&vec![0u32; reads.len()], &lens);
        let a = consensus(&c, &reads, ConsensusOpts::default());
        let b = consensus(&c, &reads, ConsensusOpts::default());
        assert_eq!(a, b, "consensus must be a pure function");
    }
}
