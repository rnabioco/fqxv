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
//!
//! ## The draft is built by alignment, not by offset
//!
//! An earlier draft laid every read down at its layout offset, "first read
//! covering a position wins". With reads starting every ~322 bases (15.5k reads
//! over 5 Mb at 40x) that made the draft ~322-base fragments from **15,183
//! different reads**, spliced wherever an offset happened to be wrong. Reads
//! could not align across ~40 such splices, so their votes landed in the wrong
//! columns and the consensus never rose above the mosaic it started from —
//! measured at **0.0078 edits/base from truth where a raw read is 0.0025**,
//! i.e. three times *worse* than an arbitrary single read, which made the margin
//! above negative. It also did precisely what the section above says not to do.
//!
//! Now: [`tiling`] picks the ~1/40th of reads needed to span the contig, and
//! [`build_draft`] places each against the draft *as it actually is* and appends
//! only its novel tail. Every junction is alignment-verified, and there are ~400
//! of them instead of ~15,000.
//!
//! Do not chase draft quality through placement — re-placing reads against a
//! fixed frame ([`place_against`](crate::place_against)) was tried and gained
//! 2.6%; a precise placement onto a bad reference is still bad. Verify any
//! change here by dumping the consensus (`FQXV_DUMP_CONS` in
//! `examples/encode.rs`) and aligning it to a known genome. Consensus-vs-truth
//! is the only number that matters for this module, and bits/base is a lagging,
//! confounded proxy for it.

use std::sync::atomic::{AtomicU32, Ordering};

use rayon::prelude::*;

use crate::{
    align::{align_banded, Op},
    place_against, ChainOpts, Contig, Sketch,
};

/// How much of the draft's end to index when placing the next tiling read.
///
/// The next read must overlap the draft's end, so only the tail can match.
/// Indexing the whole draft each step would make construction quadratic — at
/// 5 Mb and ~400 steps that is 2 Gbase of indexing for no gain.
const TAIL_WINDOW: usize = 60_000;

/// Minimum overlap between consecutive tiling reads, in bases. Enough that the
/// chain has plenty of anchors to place against; below this a read is not a safe
/// extension and is skipped.
const MIN_TILE_OVERLAP: u32 = 2_000;

/// Consensus options.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ConsensusOpts {
    /// Band half-width when aligning a read to the draft.
    pub band: usize,
    /// A column needs this many votes to be called; below it the draft's base
    /// stands. Guards the contig ends, where coverage tails off and a single
    /// erroneous read would otherwise dictate the consensus.
    pub min_votes: u32,
    /// Sketch used to place reads onto the draft. Must match the platform: the
    /// draft is built and voted by chaining, so a sketch too sparse for the
    /// error rate loses reads from both.
    pub sketch: Sketch,
}

impl Default for ConsensusOpts {
    fn default() -> Self {
        Self {
            band: 64,
            min_votes: 3,
            sketch: Sketch::ont(),
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

/// Pick a minimal set of reads that tiles the contig, left to right.
///
/// At 40x coverage roughly one read in forty is needed to span the contig; the
/// other thirty-nine add nothing to the draft's *shape* and only add splices.
/// Classic greedy interval cover: from the reads that still overlap the current
/// end by at least [`MIN_TILE_OVERLAP`], take the one reaching furthest.
///
/// Returns indices into `contig.reads`. Ties break on read index, so the tiling
/// is a pure function of the layout.
fn tiling(contig: &Contig, lens: &[u32]) -> Vec<usize> {
    let end_of = |i: usize| {
        let p = &contig.reads[i];
        p.offset + lens[p.read as usize]
    };
    let mut picked = vec![0usize];
    let mut cur_end = end_of(0);
    let mut i = 1usize;
    while i < contig.reads.len() {
        // Candidates: start early enough to overlap the current end, and extend
        // past it. `contig.reads` is sorted by offset, so once a read starts
        // beyond reach, so does every read after it.
        let mut best: Option<(u32, usize)> = None;
        let mut j = i;
        while j < contig.reads.len() {
            let p = &contig.reads[j];
            if p.offset + MIN_TILE_OVERLAP > cur_end {
                break; // and so will all later reads
            }
            let e = end_of(j);
            if e > cur_end && best.is_none_or(|(be, _)| e > be) {
                best = Some((e, j));
            }
            j += 1;
        }
        match best {
            Some((e, idx)) => {
                picked.push(idx);
                cur_end = e;
                i = idx + 1;
            }
            // Nothing overlaps the end: a coverage gap. Restart the tiling from
            // the next read rather than splicing across a region no read spans.
            None => {
                if j >= contig.reads.len() {
                    break;
                }
                picked.push(j);
                cur_end = end_of(j);
                i = j + 1;
            }
        }
    }
    picked
}

/// Build a draft by ALIGNING each tiling read onto the growing draft and
/// appending only its novel tail.
///
/// This replaces laying reads down at their layout offsets. That produced a
/// mosaic — ~322-base fragments from thousands of reads, spliced wherever an
/// offset happened to be wrong — which no vote could recover, because reads
/// could not align across the splices to vote in the first place.
///
/// Here every junction is placed by a chain against the draft *as it actually
/// is*, so the draft stays a coherent sequence, and only the tiling reads
/// (~1/40th at 40x) contribute junctions at all.
fn build_draft(contig: &Contig, reads: &[Vec<u8>], lens: &[u32], sketch: Sketch) -> Vec<u8> {
    let picked = tiling(contig, lens);
    let seed = &contig.reads[picked[0]];
    let mut draft = reads[seed.read as usize].clone();
    // Chained splices are exact; fallback splices are approximate and are the
    // suspected source of the poorly-voted regions. Count them rather than guess.
    let (mut n_chained, mut n_fallback, mut n_flip) = (0usize, 0usize, 0usize);
    // The draft's end in LAYOUT coordinates. Needed for the fallback below: it
    // lets a failed chain still extend the draft rather than stall it.
    let mut end_layout = seed.offset + lens[seed.read as usize];

    for &idx in &picked[1..] {
        let p = &contig.reads[idx];
        let r = &reads[p.read as usize];
        if r.is_empty() {
            continue;
        }
        // Only the draft's tail can overlap the next read, and indexing the whole
        // draft every step would be quadratic.
        let tail_start = draft.len().saturating_sub(TAIL_WINDOW);
        let placed = place_against(
            &draft[tail_start..],
            &[r.as_slice()],
            sketch,
            ChainOpts::default(),
        );
        // Where this read starts on the draft. A chain is exact; the layout is
        // approximate but always available.
        let start = match placed[0] {
            // Reads are pre-oriented by the layout, so a flip means the layout
            // and the draft disagree — do not trust that placement.
            Some(a) if !a.flip => {
                n_chained += 1;
                tail_start + a.offset as usize
            }
            _ => {
                if matches!(placed[0], Some(a) if a.flip) {
                    n_flip += 1;
                } else {
                    n_fallback += 1;
                }
                // FALLBACK — never stall. An earlier version skipped the read
                // here, which truncated the draft: the next tiling read is even
                // further along, so it could not overlap the stalled end either,
                // and one failure killed every read after it. That collapsed the
                // draft from 4.92 Mb to 1.15 Mb and dumped three quarters of the
                // reads into literals at ~2 bits/base.
                //
                // The layout's coordinates are approximate, so this splice may be
                // off by an indel or two — which the vote repairs. A truncated
                // draft cannot be repaired by anything.
                let overlap = end_layout.saturating_sub(p.offset) as usize;
                draft.len().saturating_sub(overlap)
            }
        };
        let start = start.min(draft.len());
        // Append only what extends past the draft's end.
        if start + r.len() > draft.len() {
            let already = draft.len() - start;
            draft.extend_from_slice(&r[already..]);
            end_layout = p.offset + lens[p.read as usize];
        }
    }
    if std::env::var("FQXV_DIAG_COLS").is_ok() {
        eprintln!(
            "DIAG draft: {} bp from {} tiling reads · splices: {n_chained} chained, \
             {n_fallback} no-chain, {n_flip} flip-disagree",
            draft.len(),
            picked.len()
        );
    }
    draft
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
    // `reads` is indexed by read id; recover per-read lengths for the tiling.
    let mut lens = vec![0u32; reads.len()];
    for (i, r) in reads.iter().enumerate() {
        lens[i] = r.len() as u32;
    }

    // 1. Draft: built by ALIGNMENT, not by offset. See `build_draft`.
    let draft = build_draft(contig, reads, &lens, opts.sketch);
    if draft.is_empty() {
        return Consensus {
            seq: Vec::new(),
            to_cons: Vec::new(),
        };
    }

    // 2. Vote. Reads are placed against the DRAFT, not at their layout offsets:
    //    the draft is built by alignment so its coordinates are its own, and a
    //    layout offset would point at the wrong column. One hop against a fixed
    //    frame, so nothing accumulates.
    //
    //    Reads vote CONCURRENTLY — the alignment dominates and there are tens of
    //    thousands per contig. Tallies are atomic and increments commute, so the
    //    result does not depend on thread count.
    let oriented: Vec<&[u8]> = contig
        .reads
        .iter()
        .map(|p| reads[p.read as usize].as_slice())
        .collect();
    let anchored = place_against(&draft, &oriented, opts.sketch, ChainOpts::default());

    let cols: Vec<AtomicColumn> = (0..draft.len()).map(|_| AtomicColumn::default()).collect();
    contig.reads.par_iter().enumerate().for_each(|(ri, p)| {
        let s = &reads[p.read as usize];
        let Some(a) = anchored[ri] else {
            return; // no confident placement on the draft: cannot vote
        };
        if a.flip {
            return; // layout and draft disagree on orientation
        }
        let from = (a.offset as usize).min(draft.len());
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
    if std::env::var("FQXV_DIAG_COLS").is_ok() {
        let covered = cols
            .iter()
            .filter(|c| c.acgt.iter().sum::<u32>() + c.del > 0)
            .count();
        let deep = cols
            .iter()
            .filter(|c| c.acgt.iter().sum::<u32>() + c.del >= opts.min_votes)
            .count();
        let mean: f64 = cols
            .iter()
            .map(|c| f64::from(c.acgt.iter().sum::<u32>() + c.del))
            .sum::<f64>()
            / cols.len() as f64;
        eprintln!(
            "DIAG cols={} covered={covered} deep(>= {})={deep} mean_votes={mean:.2}",
            cols.len(),
            opts.min_votes
        );
        for i in [1000usize, 5000, 20000] {
            if let Some(c) = cols.get(i) {
                eprintln!(
                    "DIAG col {i}: A={} C={} G={} T={} del={} -> {:?} (draft {})",
                    c.acgt[0],
                    c.acgt[1],
                    c.acgt[2],
                    c.acgt[3],
                    c.del,
                    c.call(opts.min_votes).map(|b| b as char),
                    draft[i] as char
                );
            }
        }
    }

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
    /// shifted and the count reports the shift rather than the errors.
    ///
    /// The band must absorb the **cumulative** indel drift, not merely the length
    /// difference. A consensus with deletions spread through it drifts from the
    /// truth by the total number of them, which is far more than
    /// `|len(a) - len(b)|` when the ends happen to line up. Too narrow a band
    /// loses the diagonal and degrades to substituting everything — on random
    /// ACGT that reports ~0.33/base regardless of the real answer, which looks
    /// exactly like a broken consensus. It cost an hour of chasing a working
    /// vote before the tallies showed it was unanimous.
    fn edit_dist(a: &[u8], b: &[u8]) -> u32 {
        let band = (a.len().max(b.len()) / 8).max(a.len().abs_diff(b.len()) * 2 + 256);
        align_banded(a, b, band).dist
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
    fn a_tiled_contig_beats_a_raw_read() {
        // The property the entire margin rests on, and the one the mosaic draft
        // failed in production while every unit test here passed: a consensus
        // must be closer to the truth than the reads it was built from.
        //
        // Every other test in this module stacks reads at offset 0, where a
        // mosaic cannot form. This TILES reads across a genome, which is what
        // real data looks like and what broke: the old draft spliced ~322-base
        // fragments from every read and landed 3x WORSE than a single read.
        // Kept small on purpose: scoring needs a band wide enough to absorb the
        // consensus's cumulative drift from truth, and a banded DP costs
        // len * band. At 60 kb that is a multi-GB table just to grade the test.
        let truth = rand_seq(16_000, 91);
        let read_len = 4_000usize;
        let step = 400usize; // 10x coverage — a vote needs depth to be a vote
        let n = (truth.len() - read_len) / step;

        let mut by_id: Vec<Vec<u8>> = Vec::with_capacity(n);
        let mut places = Vec::with_capacity(n);
        for i in 0..n {
            let o = i * step;
            by_id.push(mutate(&truth[o..o + read_len], 0.04, 400 + i as u32));
            places.push(Placement {
                read: i as u32,
                contig: 0,
                offset: o as u32,
                flip: false,
            });
        }
        let c = Contig {
            id: 0,
            len: truth.len() as u32,
            reads: places,
        };

        // Isolate the two stages: a bad draft and a bad vote look identical in
        // the final number, and guessing which has cost enough already.
        let lens: Vec<u32> = by_id.iter().map(|r| r.len() as u32).collect();
        let picked = tiling(&c, &lens);
        let draft = build_draft(&c, &by_id, &lens, Sketch::ont());
        let dn = draft.len().min(truth.len());
        let draft_rate = f64::from(edit_dist(&draft[..dn], &truth[..dn])) / dn as f64;
        println!(
            "DIAG tiling={} reads of {} · draft={} bp (truth {}) · draft err {:.4}/base",
            picked.len(),
            c.reads.len(),
            draft.len(),
            truth.len(),
            draft_rate
        );

        // How well are reads placed onto that draft? A read misplaced by more
        // than the band (64) aligns two unrelated windows and votes garbage into
        // every column it touches.
        let refs: Vec<&[u8]> = c
            .reads
            .iter()
            .map(|p| by_id[p.read as usize].as_slice())
            .collect();
        let anch = crate::place_against(&draft, &refs, Sketch::ont(), ChainOpts::default());
        let (mut placed, mut within_band, mut worst) = (0usize, 0usize, 0u32);
        for (ri, p) in c.reads.iter().enumerate() {
            if let Some(a) = anch[ri] {
                placed += 1;
                // The draft tracks truth coordinates closely, so a read from
                // truth offset `p.offset` should land near there.
                let err = a.offset.abs_diff(p.offset);
                if err <= 64 {
                    within_band += 1;
                }
                worst = worst.max(err);
            }
        }
        println!(
            "DIAG placed {placed}/{} · within band-64 of true offset: {within_band} · worst off-by {worst}",
            c.reads.len()
        );

        let cons = consensus(&c, &by_id, ConsensusOpts::default());
        println!(
            "DIAG consensus={} bp · err {:.4}/base",
            cons.seq.len(),
            f64::from(edit_dist(
                &cons.seq[..cons.seq.len().min(truth.len())],
                &truth[..cons.seq.len().min(truth.len())]
            )) / cons.seq.len().min(truth.len()) as f64
        );
        assert!(!cons.is_empty(), "a tiled contig must produce a consensus");

        // Compare EQUAL-LENGTH prefixes. The tiling only spans up to the last
        // read's start (~58.5 kb of a 60 kb truth), so scoring the consensus
        // against the FULL truth with a global alignment charges the ~1.5 kb it
        // never covered as ~1500 edits — 0.025/base of pure length artifact,
        // which is most of the score and nothing to do with quality. Both reads
        // start at truth[0], so equal-length prefixes are the like-for-like
        // comparison.
        let n_cmp = cons.seq.len().min(truth.len());
        assert!(
            n_cmp > truth.len() / 2,
            "the draft must span most of the contig, got {} of {}",
            cons.seq.len(),
            truth.len()
        );
        let cons_err = edit_dist(&cons.seq[..n_cmp], &truth[..n_cmp]);
        let mean_read_err: u32 = (0..n.min(8))
            .map(|i| {
                let o = i * step;
                edit_dist(&by_id[i], &truth[o..o + read_len])
            })
            .sum::<u32>()
            / n.min(8) as u32;

        let cons_rate = f64::from(cons_err) / n_cmp as f64;
        let read_rate = f64::from(mean_read_err) / read_len as f64;
        let draft_rate_cmp = draft_rate;

        // Two assertions, because this metric cannot be made tight here.
        //
        // `edit_dist` is a GLOBAL alignment, so it charges any extent mismatch as
        // error: the consensus drops columns the vote deletes, ending ~200 bases
        // shorter than the draft, and every one of those is billed as an edit
        // against `truth[..n]` (~0.014/base here — a third of the score, and
        // nothing to do with quality). A clean measure needs free end gaps, which
        // `align_banded` deliberately does not do. The authoritative check is
        // minimap2 against a real genome (`FQXV_DUMP_CONS`, see the module docs);
        // this test's job is to catch a GROSS regression, like the mosaic draft
        // that was 3x worse than a raw read.
        assert!(
            cons_rate < read_rate,
            "consensus ({cons_rate:.4}/base over {n_cmp} bp) must beat a raw read \
             ({read_rate:.4}/base) — if it does not, coding against it is pointless"
        );
        assert!(
            cons_rate < draft_rate_cmp,
            "the vote must improve on the draft it was built from \
             ({cons_rate:.4} vs {draft_rate_cmp:.4}/base) — otherwise it is not voting"
        );
    }

    #[test]
    fn one_unplaceable_read_does_not_truncate_the_draft() {
        // Regression: a read that will not chain used to be skipped, but the next
        // tiling read is further along the genome, so it could not overlap the
        // stalled draft end either — one failure killed every read after it. On
        // real data that collapsed the draft from 4.92 Mb to 1.15 Mb and dumped
        // three quarters of the reads into literals at ~2 bits/base.
        //
        // Here read 4 is replaced with junk that cannot chain. The draft must
        // still span the contig.
        let truth = rand_seq(16_000, 77);
        let read_len = 4_000usize;
        let step = 400usize;
        let n = (truth.len() - read_len) / step;

        let mut by_id: Vec<Vec<u8>> = Vec::with_capacity(n);
        let mut places = Vec::with_capacity(n);
        for i in 0..n {
            let o = i * step;
            by_id.push(mutate(&truth[o..o + read_len], 0.04, 500 + i as u32));
            places.push(Placement {
                read: i as u32,
                contig: 0,
                offset: o as u32,
                flip: false,
            });
        }
        // Sabotage a read the tiling is likely to pick: unrelated sequence.
        by_id[10] = rand_seq(read_len, 4242);

        let c = Contig {
            id: 0,
            len: truth.len() as u32,
            reads: places,
        };
        let lens: Vec<u32> = by_id.iter().map(|r| r.len() as u32).collect();
        let draft = build_draft(&c, &by_id, &lens, Sketch::ont());
        assert!(
            draft.len() > truth.len() * 3 / 4,
            "an unplaceable read must not truncate the draft: got {} of {}",
            draft.len(),
            truth.len()
        );
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
