//! Re-place reads against a fixed reference.
//!
//! The greedy layout composes offsets hop by hop (`off = cur_off + q_start -
//! t_start`), so every hop's indel error is added to the last and nothing
//! re-anchors them. This module removes that by construction: the reference is a
//! **fixed frame**, so each read is anchored to it directly, in one hop, and
//! error cannot compound. The layout's job is only to decide *which* reads
//! belong together and roughly where — which it does well (99%+ recall);
//! precision is this module's job.
//!
//! It is deliberately not a wider band. A band wide enough to swallow
//! accumulated drift pays for it quadratically and still only *masks* the
//! problem.
//!
//! ## What this did NOT fix
//!
//! This was built to close the 0.427-vs-0.040 bits/base gap, on the theory that
//! accumulated placement drift caused the ~14x phantom edits (0.036 subs/base
//! emitted where 0.0025 exist). **That theory was wrong.** Re-placing gained
//! 2.6% (0.4267 -> 0.4154).
//!
//! Measuring the intermediate artifact rather than theorising about the pipeline
//! found the real cause: the CONSENSUS is 0.0078 edits/base from the truth
//! versus 0.0025 for a raw read — three times worse than an arbitrary single
//! read (dump it with `FQXV_DUMP_CONS` and align to the reference). Re-placing
//! precisely against a bad reference cannot help, and no downstream change can.
//! See [`consensus`](crate::consensus) for the mosaic draft that causes it.
//!
//! This module is still correct and still needed — a fixed frame is the right
//! way to place — it simply was not the bottleneck.

use rayon::prelude::*;

use crate::{find_overlaps, ChainOpts, Index, Repeat, Sketch};

/// Where a read sits on a reference, found in a single hop.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Anchored {
    /// The read.
    pub read: u32,
    /// Offset of the read's start along the reference. Clamped at 0: a read
    /// hanging off the start codes its overhang as an insertion.
    pub offset: u32,
    /// True when the read aligns reverse-complemented to the reference.
    pub flip: bool,
    /// Chain score — how much to trust this placement.
    pub score: i32,
    /// How far the chain's diagonal moved between its first and last anchor, in
    /// bases: the net indel between read and reference across the chained span.
    ///
    /// [`offset`](Self::offset) anchors the read on the chain's FIRST diagonal,
    /// so an aligner working from it must be allowed to travel this far
    /// off-diagonal or it cannot reach the alignment that exists. That is a
    /// per-read number and it varies by orders of magnitude — most reads drift a
    /// handful of bases, a few drift hundreds — so a caller without it has to
    /// pick one band for every read and pay the worst case on all of them.
    /// Measured: a global band of 384 costs 4x the DP work of 96 to buy 16% of
    /// ratio, because 96 could not reach the tail of this distribution.
    ///
    /// It is a floor, not a bound: it measures drift only across the chained
    /// span, and says nothing about the read's unanchored ends. Add a margin.
    pub drift: u32,
}

/// Re-place `reads` against `refseq`, one hop each.
///
/// `reads[i]` is read `i`'s sequence in its ORIGINAL orientation; the returned
/// [`Anchored::flip`] says whether it aligns reverse-complemented. Reads with no
/// confident chain return `None` and are the caller's problem — they code
/// standalone, which is the graceful-degradation path.
///
/// Returns one entry per input read, in input order.
#[must_use]
pub fn place_against(
    refseq: &[u8],
    reads: &[&[u8]],
    sketch: Sketch,
    opts: ChainOpts,
) -> Vec<Option<Anchored>> {
    if refseq.is_empty() {
        return vec![None; reads.len()];
    }
    // Index the reference alone. Repeat filtering is off: a reference is one
    // copy of the locus, so its minimizers are not repetitive in the sense the
    // filter targets, and dropping any would blind reads to that stretch.
    let Ok(idx) = Index::build(
        &[refseq.len() as u32],
        refseq,
        sketch,
        Repeat { drop_top_frac: 0.0 },
    ) else {
        return vec![None; reads.len()];
    };

    reads
        .par_iter()
        .enumerate()
        .map(|(i, read)| {
            // The reference is read 0 in `idx`; pass a query id that is not in
            // the index so nothing is excluded as self.
            let ov = find_overlaps(&idx, u32::MAX, read, opts);
            let best = ov.first()?;
            // One hop from a fixed frame: read[q_start] aligns to ref[t_start],
            // so read[0] sits at t_start - q_start. No composition, no drift.
            let off = i64::from(best.t_start) - i64::from(best.q_start);
            // `find_overlaps` reports t_start in the QUERY's frame, so for a
            // reverse-complement hit that coordinate is in RC(reference), not
            // the reference. Map it back: a span at `off` in RC(ref) occupies
            // `[reflen - off - readlen, reflen - off)` in ref. Skipping this
            // silently places every RC read at the mirror of where it belongs —
            // which on a real contig is a completely different locus.
            let offset = if best.strand {
                let reflen = refseq.len() as i64;
                (reflen - off - read.len() as i64).max(0) as u32
            } else {
                off.max(0) as u32
            };
            // The chain's diagonal at its last anchor against its diagonal at
            // the first: the net indel across the chained span. Taken in the
            // query's frame, where both coordinates already live, and only its
            // magnitude is used — reflecting into the reference's frame for an
            // RC hit negates the shift but cannot change how far it is.
            let d_end = i64::from(best.t_end) - i64::from(best.q_end);
            let drift = (d_end - off).unsigned_abs() as u32;
            Some(Anchored {
                read: i as u32,
                offset,
                flip: best.strand,
                score: best.score,
                drift,
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Pseudo-random bases via splitmix64.
    ///
    /// NOT the LCG the other test modules use. An LCG mod 2^32 has period
    /// `2^(k+1)` in bit `k`, so the usual `(x >> 16) % 4` draws bits 16-17 and
    /// repeats every **2^18 = 262144** bases. Every other module's fixtures are
    /// shorter than that, so it never showed; this module needs a 400 kb
    /// reference, and the periodic "random" sequence made it a tandem repeat —
    /// a read at 380 kb genuinely matched an identical copy at 117856 and
    /// `place_against` correctly found it, failing the test by exactly 262144.
    ///
    /// A generator whose output silently repeats is not a fixture, it is a
    /// second bug hiding the first.
    fn rand_seq(n: usize, seed: u32) -> Vec<u8> {
        let mut x = u64::from(seed);
        (0..n)
            .map(|_| {
                x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
                let mut z = x;
                z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
                z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
                z ^= z >> 31;
                b"ACGT"[(z % 4) as usize]
            })
            .collect()
    }

    #[test]
    fn the_test_generator_is_not_periodic() {
        // Guards the fixture itself: the LCG this replaced repeated every 2^18
        // bases, which turned a "random" reference into a tandem repeat and made
        // a correct placement look like a 262144-base error.
        let s = rand_seq(600_000, 1);
        let period = 262_144;
        let a = &s[0..4_000];
        let b = &s[period..period + 4_000];
        let same = a.iter().zip(b).filter(|(x, y)| x == y).count();
        // Unrelated ACGT agrees ~25% of the time; a repeat agrees 100%.
        assert!(
            same < a.len() / 2,
            "generator repeats at 2^18: {same}/{} bases identical",
            a.len()
        );
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

    #[test]
    fn empty_reference_places_nothing() {
        let r = rand_seq(500, 1);
        let got = place_against(b"", &[&r], Sketch::ont(), ChainOpts::default());
        assert_eq!(got, vec![None]);
    }

    #[test]
    fn unrelated_reads_are_not_placed() {
        // The graceful-degradation path: no confident chain -> code standalone.
        let refseq = rand_seq(20_000, 3);
        let alien = rand_seq(5_000, 999);
        let got = place_against(&refseq, &[&alien], Sketch::ont(), ChainOpts::default());
        assert_eq!(got, vec![None], "unrelated sequence must not be placed");
    }

    #[test]
    fn places_reads_at_their_true_offsets() {
        // The property the whole module exists for: an offset accurate enough
        // that a narrow alignment band finds the read where it was put.
        let refseq = rand_seq(40_000, 5);
        let truth: Vec<u32> = vec![0, 3_000, 9_500, 21_000, 33_000];
        let reads: Vec<Vec<u8>> = truth
            .iter()
            .enumerate()
            .map(|(i, &o)| {
                mutate(
                    &refseq[o as usize..o as usize + 5_000],
                    0.03,
                    100 + i as u32,
                )
            })
            .collect();
        let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        let got = place_against(&refseq, &refs, Sketch::ont(), ChainOpts::default());

        for (i, want) in truth.iter().enumerate() {
            let a = got[i].expect("read must place");
            assert!(!a.flip, "read {i} is forward");
            let err = a.offset.abs_diff(*want);
            assert!(
                err <= 32,
                "read {i}: placed at {} but belongs at {want} (off by {err}) — \
                 a band-64 align cannot absorb more than this",
                a.offset
            );
        }
    }

    #[test]
    fn places_reverse_complemented_reads() {
        let refseq = rand_seq(30_000, 7);
        let fwd = mutate(&refseq[8_000..14_000], 0.03, 11);
        let rc = revcomp(&fwd);
        let got = place_against(&refseq, &[&rc], Sketch::ont(), ChainOpts::default());
        let a = got[0].expect("an RC read must place");
        assert!(a.flip, "must be flagged reverse-complement");
        assert!(
            a.offset.abs_diff(8_000) <= 64,
            "RC read placed at {} but belongs at ~8000",
            a.offset
        );
    }

    #[test]
    fn placement_error_does_not_grow_with_distance() {
        // THE point: against a fixed frame, a read at 200 kb is placed exactly as
        // accurately as one at 2 kb. The greedy layout's error compounds with
        // every hop, which is what this replaces.
        let refseq = rand_seq(400_000, 13);
        let near = mutate(&refseq[2_000..8_000], 0.03, 21);
        let far = mutate(&refseq[380_000..386_000], 0.03, 22);
        let refs: Vec<&[u8]> = vec![&near, &far];
        let got = place_against(&refseq, &refs, Sketch::ont(), ChainOpts::default());

        let near_err = got[0].expect("near places").offset.abs_diff(2_000);
        let far_err = got[1].expect("far places").offset.abs_diff(380_000);
        assert!(near_err <= 32, "near off by {near_err}");
        assert!(
            far_err <= 32,
            "far off by {far_err} — placement error must NOT grow with distance"
        );
    }

    #[test]
    fn is_deterministic() {
        let refseq = rand_seq(20_000, 17);
        let reads: Vec<Vec<u8>> = (0..6)
            .map(|i| mutate(&refseq[i * 2_000..i * 2_000 + 5_000], 0.04, 300 + i as u32))
            .collect();
        let refs: Vec<&[u8]> = reads.iter().map(|r| r.as_slice()).collect();
        let a = place_against(&refseq, &refs, Sketch::ont(), ChainOpts::default());
        let b = place_against(&refseq, &refs, Sketch::ont(), ChainOpts::default());
        assert_eq!(a, b, "re-placement must be a pure function");
    }
}
