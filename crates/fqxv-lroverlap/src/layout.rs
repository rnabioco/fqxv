//! Greedy layout: overlaps -> reads placed at offsets along contigs.
//!
//! Given the overlap graph, decide where each read sits. The output is what the
//! consensus is voted over and what each read's edit script is coded against.
//!
//! This is deliberately *not* an assembler. It does not need to produce a
//! biologically correct genome — only a layout good enough that reads placed on
//! it agree, because the consensus is a coding reference, not a result. Anything
//! it cannot place is coded standalone, so a bad decision costs bits, never
//! correctness.
//!
//! ## Determinism
//!
//! Reads are seeded and extended in a fixed order derived from the input, every
//! candidate set is sorted by a total order, and no map is iterated — so the
//! layout is a pure function of the overlaps.

use crate::Overlap;

/// Where a read sits on a contig.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Placement {
    /// The read.
    pub read: u32,
    /// Contig it was placed on.
    pub contig: u32,
    /// Offset of the read's start along the contig. Signed during construction;
    /// contigs are normalized to start at 0, so this is non-negative.
    pub offset: u32,
    /// True when the read is stored reverse-complemented relative to the contig.
    pub flip: bool,
}

/// A laid-out contig: the reads on it, in ascending offset.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Contig {
    /// Contig id.
    pub id: u32,
    /// Total span in bases.
    pub len: u32,
    /// Reads placed on it.
    pub reads: Vec<Placement>,
}

/// Layout knobs.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct LayoutOpts {
    /// Minimum chain score for an overlap to be trusted for placement. Higher
    /// than the detection threshold: a weak overlap is fine as a *candidate* but
    /// a bad anchor for a placement, and a misplacement costs every read
    /// downstream of it.
    pub min_score: i32,
    /// Maximum overlaps considered per read when extending, best-scoring first.
    /// Bounds the work; the best few are the only plausible extensions anyway.
    pub max_candidates: usize,
}

impl Default for LayoutOpts {
    fn default() -> Self {
        Self {
            min_score: 200,
            max_candidates: 16,
        }
    }
}

/// Build a greedy layout from per-read overlaps.
///
/// `overlaps[i]` are read `i`'s overlaps (as returned by
/// [`find_overlaps`](crate::find_overlaps), best-scoring first). `lens` gives
/// each read's length.
///
/// Every read lands on exactly one contig; a read with no usable overlap becomes
/// a singleton contig, which is the graceful-degradation path — it is simply
/// coded standalone.
#[must_use]
pub fn layout(lens: &[u32], overlaps: &[Vec<Overlap>], opts: LayoutOpts) -> Vec<Contig> {
    let n = lens.len();
    // Contig assignment and offset within it, filled as reads are placed.
    let mut contig_of: Vec<Option<u32>> = vec![None; n];
    let mut offset_of: Vec<i64> = vec![0; n];
    let mut flip_of: Vec<bool> = vec![false; n];
    let mut members: Vec<Vec<u32>> = Vec::new();

    // Seed order: reads with the most high-scoring overlaps first — they sit in
    // well-covered regions and make the most reliable seeds. Ties on read index,
    // so the order is total and input-derived.
    let mut seeds: Vec<(usize, u32)> = (0..n)
        .map(|i| {
            let deg = overlaps[i]
                .iter()
                .filter(|o| o.score >= opts.min_score)
                .count();
            (deg, i as u32)
        })
        .collect();
    seeds.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));

    for &(_, seed) in &seeds {
        if contig_of[seed as usize].is_some() {
            continue;
        }
        // New contig rooted at this read.
        let cid = members.len() as u32;
        members.push(vec![seed]);
        contig_of[seed as usize] = Some(cid);
        offset_of[seed as usize] = 0;
        flip_of[seed as usize] = false;

        // Breadth-first extension from placed reads. A queue, not recursion:
        // contigs reach tens of thousands of reads at high coverage.
        let mut queue = vec![seed];
        while let Some(cur) = queue.pop() {
            let cur_off = offset_of[cur as usize];
            let cur_flip = flip_of[cur as usize];
            let lc = i64::from(lens[cur as usize]);
            for o in overlaps[cur as usize]
                .iter()
                .filter(|o| o.score >= opts.min_score)
                .take(opts.max_candidates)
            {
                if contig_of[o.target as usize].is_some() {
                    continue; // already placed; do not re-litigate
                }
                let lt = i64::from(lens[o.target as usize]);
                let qs = i64::from(o.q_start);
                let ts = i64::from(o.t_start);

                // `o.t_start` is a coordinate in the ORIENTED target — the target
                // itself when `!o.strand`, its reverse complement when
                // `o.strand` (`find_overlaps` mapped it into the query's frame).
                // The oriented target is therefore colinear with `cur`'s
                // ORIGINAL sequence, which is what was queried.
                let off = if !cur_flip {
                    // The contig holds `cur` as-is, so contig coords are cur's
                    // coords shifted: the oriented target's base 0 lands where
                    // cur's base (q_start - t_start) is.
                    cur_off + qs - ts
                } else {
                    // The contig holds RC(cur), so it runs opposite to cur's
                    // original coords and the oriented target is laid down
                    // reversed too — i.e. RC(oriented). cur's original base p
                    // sits at contig `cur_off + (lc - 1 - p)`, so aligning
                    // cur[q_start] with oriented[t_start] and reversing both
                    // puts RC(oriented)'s base 0 here:
                    cur_off + (lc - 1 - qs) - (lt - 1 - ts)
                };
                // Flip composes: passing through a flipped read inverts the
                // orientation the overlap reports.
                let flip = cur_flip ^ o.strand;

                contig_of[o.target as usize] = Some(cid);
                offset_of[o.target as usize] = off;
                flip_of[o.target as usize] = flip;
                members[cid as usize].push(o.target);
                queue.push(o.target);
            }
        }
    }

    // Any read never reached becomes a singleton — coded standalone.
    for r in 0..n as u32 {
        if contig_of[r as usize].is_none() {
            let cid = members.len() as u32;
            members.push(vec![r]);
            contig_of[r as usize] = Some(cid);
            offset_of[r as usize] = 0;
        }
    }

    // Normalize each contig to start at offset 0 and sort members by offset.
    let mut out = Vec::with_capacity(members.len());
    for (cid, mem) in members.into_iter().enumerate() {
        let min_off = mem
            .iter()
            .map(|&r| offset_of[r as usize])
            .min()
            .unwrap_or(0);
        let mut reads: Vec<Placement> = mem
            .iter()
            .map(|&r| Placement {
                read: r,
                contig: cid as u32,
                offset: (offset_of[r as usize] - min_off) as u32,
                flip: flip_of[r as usize],
            })
            .collect();
        // Total order: offset, then read index.
        reads.sort_unstable_by(|a, b| a.offset.cmp(&b.offset).then(a.read.cmp(&b.read)));
        let len = reads
            .iter()
            .map(|p| p.offset + lens[p.read as usize])
            .max()
            .unwrap_or(0);
        out.push(Contig {
            id: cid as u32,
            len,
            reads,
        });
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{find_overlaps, ChainOpts, Index, Repeat, Sketch};

    fn rand_seq(n: usize, seed: u32) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                b"ACGT"[(x >> 16) as usize % 4]
            })
            .collect()
    }

    fn revcomp(s: &[u8]) -> Vec<u8> {
        s.iter()
            .rev()
            .map(|&b| match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                x => x,
            })
            .collect()
    }

    fn ov(target: u32, score: i32, q_start: u32, t_start: u32) -> Overlap {
        ov_s(target, score, q_start, t_start, false)
    }

    fn ov_s(target: u32, score: i32, q_start: u32, t_start: u32, strand: bool) -> Overlap {
        Overlap {
            target,
            strand,
            score,
            n_anchors: 10,
            q_start,
            q_end: q_start + 1000,
            t_start,
            t_end: t_start + 1000,
        }
    }

    fn find(c: &[Contig], read: u32) -> Placement {
        *c.iter()
            .flat_map(|x| x.reads.iter())
            .find(|p| p.read == read)
            .expect("read placed")
    }

    #[test]
    fn no_overlaps_gives_all_singletons() {
        // The graceful-degradation path: every read codes standalone.
        let lens = [100u32, 100, 100];
        let ovs = vec![Vec::new(), Vec::new(), Vec::new()];
        let c = layout(&lens, &ovs, LayoutOpts::default());
        assert_eq!(c.len(), 3, "each unplaceable read must be its own contig");
        assert!(c.iter().all(|x| x.reads.len() == 1));
    }

    #[test]
    fn every_read_is_placed_exactly_once() {
        // The invariant the codec depends on: no read lost, none duplicated.
        let lens = [1000u32; 5];
        let ovs = vec![
            vec![ov(1, 500, 400, 0)],
            vec![ov(0, 500, 0, 400), ov(2, 500, 400, 0)],
            vec![ov(1, 500, 0, 400)],
            vec![],
            vec![],
        ];
        let c = layout(&lens, &ovs, LayoutOpts::default());
        let mut seen: Vec<u32> = c
            .iter()
            .flat_map(|x| x.reads.iter().map(|p| p.read))
            .collect();
        seen.sort_unstable();
        assert_eq!(seen, vec![0, 1, 2, 3, 4], "every read exactly once");
    }

    #[test]
    fn a_chain_of_overlaps_forms_one_contig() {
        // Three reads tiling: 0 -> 1 -> 2, each 400 further along.
        let lens = [1000u32; 3];
        let ovs = vec![
            vec![ov(1, 500, 400, 0)],
            vec![ov(0, 500, 0, 400), ov(2, 500, 400, 0)],
            vec![ov(2, 500, 0, 400)],
        ];
        let c = layout(&lens, &ovs, LayoutOpts::default());
        let big = c
            .iter()
            .find(|x| x.reads.len() > 1)
            .expect("a contig formed");
        assert_eq!(big.reads.len(), 3, "all three must land on one contig");
        // Offsets must be 0, 400, 800 — ascending and normalized to start at 0.
        let offs: Vec<u32> = big.reads.iter().map(|p| p.offset).collect();
        assert_eq!(
            offs,
            vec![0, 400, 800],
            "offsets must compose along the chain"
        );
        assert_eq!(big.len, 1800, "span = last offset + its read length");
    }

    #[test]
    fn weak_overlaps_are_not_used_for_placement() {
        // Below min_score: a weak overlap is a fine candidate but a bad anchor,
        // and a misplacement costs every read downstream of it.
        let lens = [1000u32; 2];
        let ovs = vec![vec![ov(1, 10, 400, 0)], vec![ov(0, 10, 0, 400)]];
        let c = layout(&lens, &ovs, LayoutOpts::default());
        assert_eq!(c.len(), 2, "weak overlaps must not merge contigs");
    }

    #[test]
    fn a_reverse_complement_overlap_places_the_target_flipped() {
        // Roughly half of all overlaps are RC — reads come off the sequencer in
        // both orientations. Dropping them (as an earlier cut of this function
        // did) halves the usable graph and shatters contigs.
        let lens = [1000u32; 2];
        let ovs = vec![
            vec![ov_s(1, 500, 400, 0, true)],
            vec![ov_s(0, 500, 0, 400, true)],
        ];
        let c = layout(&lens, &ovs, LayoutOpts::default());
        assert_eq!(c.len(), 1, "an RC overlap must still join one contig");
        assert!(!find(&c, 0).flip, "the seed is never flipped");
        assert!(find(&c, 1).flip, "an RC target must be stored flipped");
    }

    #[test]
    fn flip_composes_through_a_flipped_read() {
        // 0 -RC-> 1 -RC-> 2: read 1 is opposite both its neighbours, and read 2
        // is RC of an RC read, so it shares read 0's orientation. Flip composes
        // as XOR; getting that wrong silently stores a read's complement.
        //
        // Assert RELATIVE orientation only. Which read seeds the contig is
        // chosen by overlap degree (read 1 here, with two), and a contig's
        // overall orientation is arbitrary — so absolute flips are meaningless
        // and only the relationships between them carry information.
        let lens = [1000u32; 3];
        let ovs = vec![
            vec![ov_s(1, 500, 400, 0, true)],
            vec![ov_s(0, 500, 0, 400, true), ov_s(2, 500, 400, 0, true)],
            vec![ov_s(1, 500, 0, 400, true)],
        ];
        let c = layout(&lens, &ovs, LayoutOpts::default());
        assert_eq!(c.len(), 1, "all three must land on one contig");
        let (f0, f1, f2) = (find(&c, 0).flip, find(&c, 1).flip, find(&c, 2).flip);
        assert_eq!(f0, f2, "reads 0 and 2 must share an orientation (RC of RC)");
        assert_ne!(f1, f0, "read 1 must be opposite its neighbours");
    }

    #[test]
    fn placement_reconstructs_the_contig_sequence() {
        // The real invariant: at every placed read's offset, the read (flipped
        // as recorded) must agree with what the contig already holds there. This
        // is what makes the consensus meaningful — and it exercises the RC
        // offset arithmetic rather than just the flip flag.
        let genome = rand_seq(6000, 77);
        // Read 0 = [0,3000). Read 1 = RC of [1500,4500) — overlapping and flipped.
        let r0 = genome[0..3000].to_vec();
        let r1 = revcomp(&genome[1500..4500]);
        let lens = [r0.len() as u32, r1.len() as u32];

        // Read 1 oriented back (RC of r1) is genome[1500,4500), which aligns to
        // read 0 at q_start 1500 <-> oriented t_start 0.
        let ovs = vec![
            vec![ov_s(1, 500, 1500, 0, true)],
            vec![ov_s(0, 500, 0, 1500, true)],
        ];
        let c = layout(&lens, &ovs, LayoutOpts::default());
        assert_eq!(c.len(), 1);

        // Lay both reads down in contig coordinates and check they agree where
        // they overlap.
        let seqs = [r0, r1];
        let mut canvas: Vec<Option<u8>> = vec![None; c[0].len as usize];
        let mut checked = 0usize;
        for p in &c[0].reads {
            let s = if p.flip {
                revcomp(&seqs[p.read as usize])
            } else {
                seqs[p.read as usize].clone()
            };
            for (i, &b) in s.iter().enumerate() {
                let at = p.offset as usize + i;
                match canvas[at] {
                    None => canvas[at] = Some(b),
                    Some(prev) => {
                        assert_eq!(prev, b, "reads disagree at contig position {at}");
                        checked += 1;
                    }
                }
            }
        }
        assert!(
            checked > 1000,
            "expected a substantial overlap to verify, checked {checked} bases"
        );
    }

    #[test]
    fn is_deterministic() {
        let lens = [1000u32; 6];
        let ovs = vec![
            vec![ov(1, 500, 400, 0), ov(2, 450, 800, 0)],
            vec![ov(0, 500, 0, 400), ov(2, 500, 400, 0)],
            vec![ov(1, 500, 0, 400), ov(0, 450, 0, 800)],
            vec![ov(4, 600, 200, 0)],
            vec![ov(3, 600, 0, 200)],
            vec![],
        ];
        let a = layout(&lens, &ovs, LayoutOpts::default());
        let b = layout(&lens, &ovs, LayoutOpts::default());
        assert_eq!(a, b, "layout must be a pure function of the overlaps");
    }

    #[test]
    fn end_to_end_mixed_strand_reads_form_one_contig() {
        // The test that should have existed from the start: real reads through
        // the real overlap search, with every other read reverse-complemented as
        // they come off a sequencer. If RC edges are dropped, this shatters into
        // many contigs instead of one.
        let genome = rand_seq(30_000, 88);
        let reads: Vec<Vec<u8>> = (0..10)
            .map(|i| {
                let s = genome[i * 2000..i * 2000 + 6000].to_vec();
                if i % 2 == 1 {
                    revcomp(&s)
                } else {
                    s
                }
            })
            .collect();
        let lens: Vec<u32> = reads.iter().map(|r| r.len() as u32).collect();
        let seq: Vec<u8> = reads.iter().flat_map(|r| r.iter().copied()).collect();
        let idx = Index::build(&lens, &seq, Sketch::ont(), Repeat { drop_top_frac: 0.0 }).unwrap();

        let mut offs = vec![0usize];
        for &l in &lens {
            offs.push(offs.last().unwrap() + l as usize);
        }
        let ovs: Vec<Vec<Overlap>> = (0..reads.len())
            .map(|r| {
                find_overlaps(
                    &idx,
                    r as u32,
                    &seq[offs[r]..offs[r + 1]],
                    ChainOpts::default(),
                )
            })
            .collect();

        let c = layout(&lens, &ovs, LayoutOpts::default());
        assert_eq!(
            c.len(),
            1,
            "mixed-strand reads tiling one genome must form ONE contig, got {} \
             (dropping RC edges is what shatters this)",
            c.len()
        );
        assert_eq!(c[0].reads.len(), 10, "every read on the contig");
        // Odd reads were RC'd, so odd and even reads must end up in OPPOSITE
        // orientations. Which group is `flip=true` depends on the seed and is
        // meaningless, so compare the groups rather than absolute flags.
        let flip_of = |r: u32| find(&c, r).flip;
        let even = flip_of(0);
        for r in 0..10u32 {
            let want = if r % 2 == 0 { even } else { !even };
            assert_eq!(
                flip_of(r),
                want,
                "read {r} orientation must match its strand group"
            );
        }

        // And the placement must be geometrically right: reads laid down in
        // contig order must match the genome tiling order (forwards or
        // backwards, since the seed's orientation is arbitrary).
        let mut by_off: Vec<u32> = c[0].reads.iter().map(|p| p.read).collect();
        let fwd: Vec<u32> = (0..10).collect();
        let rev: Vec<u32> = (0..10).rev().collect();
        assert!(
            by_off == fwd || by_off == rev,
            "reads must lay out in tiling order, got {by_off:?}"
        );
        by_off.sort_unstable();
        assert_eq!(by_off, fwd, "no read lost or duplicated");
    }

    #[test]
    fn offsets_are_normalized_to_zero() {
        // Extension can run backwards (a target starting before its seed), so
        // offsets are signed during construction and normalized after.
        let lens = [1000u32; 2];
        // Read 1 sits BEFORE read 0: q_start 0 aligns to t_start 400.
        let ovs = vec![vec![ov(1, 500, 0, 400)], vec![ov(0, 500, 400, 0)]];
        let c = layout(&lens, &ovs, LayoutOpts::default());
        let big = c.iter().find(|x| x.reads.len() > 1).expect("contig");
        assert_eq!(big.reads[0].offset, 0, "the leftmost read must sit at 0");
        assert!(big.reads.iter().all(|p| p.offset < big.len));
    }
}
