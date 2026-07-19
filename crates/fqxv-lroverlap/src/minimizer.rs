//! Canonical `(w, k)` minimizers.
//!
//! A minimizer is the smallest hashed k-mer in every window of `w` consecutive
//! k-mers, so roughly `2/(w+1)` of positions are selected while any two
//! sequences sharing a long enough exact stretch are guaranteed to select the
//! same k-mer from it. That guarantee is what survives a noisy channel: a 14 kb
//! ONT read carries ~2500 minimizers at `w = 10`, so even at ~10% error (each
//! specific 15-mer surviving with probability `0.9^15 ≈ 0.21`) hundreds of
//! shared minimizers remain between two truly-overlapping reads — enough to
//! chain. Contrast a single global anchor per read, which survives with that
//! same 0.21 and has no redundancy at all.
//!
//! Minimizers are **canonical**: the k-mer and its reverse complement hash to
//! one value, and a strand flag records which orientation was seen. Two reads
//! from opposite strands of the same locus therefore select the same
//! minimizers.

/// One selected minimizer occurrence within a sequence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Minimizer {
    /// Hash of the canonical k-mer.
    pub hash: u64,
    /// Offset of the k-mer's first base within the sequence.
    pub pos: u32,
    /// True when the reverse-complement orientation was the canonical one.
    pub strand: bool,
}

/// 2-bit code for a base; `None` for anything outside ACGT (N, IUPAC, …), which
/// breaks the k-mer run rather than silently coding as `A`.
#[inline]
fn code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

/// Invertible integer hash (splitmix64). Hashing decorrelates minimizer
/// selection from raw k-mer order, which otherwise biases selection toward
/// poly-A: the all-A k-mer packs to `0`, the smallest possible code, so it would
/// win every window it appears in.
///
/// The golden-ratio increment is load-bearing, not decoration. Without it the
/// bare finalizer maps `0 -> 0`, so the all-A k-mer would hash to the *minimum
/// possible value* and homopolymer runs would remain minimizer magnets — the
/// exact bias this is here to remove.
#[inline]
pub(crate) fn hash64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9e37_79b9_7f4a_7c15);
    x = (x ^ (x >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    x ^ (x >> 31)
}

/// Compute the canonical `(w, k)` minimizers of `seq`.
///
/// `k` must be in `1..=31` (a k-mer is packed into a `u64`) and `w >= 1`.
/// Returns occurrences in ascending `pos` order, with consecutive same-**hash**
/// selections collapsed — a window slide usually re-selects the same k-mer, and
/// storing it once per run is what keeps the index small.
///
/// Collapsing on the hash rather than the whole `(hash, pos, strand)` triple is
/// what makes selection symmetric under reverse complement, which the canonical
/// property depends on. The window-minimum *hash* sequence of a revcomp is
/// exactly the forward sequence reversed, and run structure is invariant under
/// reversal — so deduplicating runs by hash yields the same multiset either way.
/// Deduplicating on `pos` too would not: among equal-hash candidates the deque
/// keeps the latest position, and "latest" flips with direction, so a repeated
/// k-mer (a tandem repeat) would emit a different number of times forward than
/// reversed.
///
/// Runs of non-ACGT reset the k-mer accumulator, so an `N` run simply produces
/// no minimizers rather than spurious ones.
#[must_use]
pub fn minimizers(seq: &[u8], w: usize, k: usize) -> Vec<Minimizer> {
    assert!((1..=31).contains(&k), "k must be 1..=31");
    assert!(w >= 1, "w must be >= 1");
    let mut out = Vec::new();
    if seq.len() < k {
        return out;
    }
    let mask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let shift = 2 * (k - 1);

    // Monotonic deque of candidate (hash, pos, strand), increasing by hash, so
    // the window minimum is always at the front. O(n) overall.
    let mut dq: std::collections::VecDeque<Minimizer> = std::collections::VecDeque::new();
    let (mut fwd, mut rc) = (0u64, 0u64);
    let mut have = 0usize; // consecutive valid bases so far
    let mut last_hash: Option<u64> = None;

    for (i, &b) in seq.iter().enumerate() {
        let Some(c) = code(b) else {
            // Break the run: no k-mer spans a non-ACGT base.
            have = 0;
            fwd = 0;
            rc = 0;
            dq.clear();
            continue;
        };
        fwd = ((fwd << 2) | c) & mask;
        rc = (rc >> 2) | ((3 - c) << shift);
        have += 1;
        if have < k {
            continue;
        }
        let pos = (i + 1 - k) as u32;
        // Canonical: the smaller of the two orientations. Ties (a palindromic
        // k-mer) resolve to the forward strand, on both encode and decode.
        let (canon, strand) = if fwd <= rc { (fwd, false) } else { (rc, true) };
        let cand = Minimizer {
            hash: hash64(canon),
            pos,
            strand,
        };

        // Evict candidates that can never again be the minimum.
        while dq.back().is_some_and(|m| m.hash >= cand.hash) {
            dq.pop_back();
        }
        dq.push_back(cand);
        // Drop the front once it leaves the window of the last `w` k-mers.
        let win_start = (pos + 1).saturating_sub(w as u32);
        while dq.front().is_some_and(|m| m.pos < win_start) {
            dq.pop_front();
        }
        // Emit only once the first full window exists.
        if pos + 1 >= w as u32 {
            let m = *dq.front().expect("deque non-empty after push");
            if last_hash != Some(m.hash) {
                out.push(m);
                last_hash = Some(m.hash);
            }
        }
    }
    out
}

/// Compute the canonical **closed syncmers** of `seq` for parameters `(k, s)`.
///
/// A closed syncmer is a k-mer whose minimal `s`-mer (by canonical hash, over the
/// `k - s + 1` s-mers it contains) sits at the **first or last** s-mer position.
/// Like [`minimizers`] the emitted `hash`/`strand` are for the canonical k-mer and
/// runs of non-ACGT reset the accumulator, so the two are interchangeable anchors.
///
/// **Why syncmers on noisy (ONT) reads.** A minimizer is selected only if it wins a
/// `w`-wide window, so a base error in any *neighbor* k-mer can deselect an intact
/// k-mer — the dominant conservation loss at ~10% error. A syncmer's selection is a
/// function of the k-mer's *own* content only, so an intact shared k-mer is
/// co-selected by two reads regardless of surrounding errors. At matched density and
/// the same `k` this recovers ~1.3–1.6× more shared anchors without changing
/// specificity or the hot-path anchor count.
///
/// **Strand symmetry.** Selection uses *canonical* s-mer hashes and the symmetric
/// {first, last} offset pair: under reverse-complement a k-mer's s-mer-hash sequence
/// reverses (canonical s-mer hashes are revcomp-invariant), so "min at first or last"
/// is invariant, and the run-dedup-by-hash keeps the multiset identical either
/// direction — exactly as [`minimizers`] guarantees.
///
/// **Density** is `2 / (k - s + 1)` (min of `k-s+1` positions lands on an endpoint),
/// so `s = k - w` matches [`minimizers`]' `2 / (w + 1)`. Requires `1 <= s < k`.
#[must_use]
pub fn syncmers(seq: &[u8], k: usize, s: usize) -> Vec<Minimizer> {
    assert!((1..=31).contains(&k), "k must be 1..=31");
    assert!(s >= 1 && s < k, "s must be in 1..k");
    let mut out = Vec::new();
    if seq.len() < k {
        return out;
    }
    let kmask = if k == 32 {
        u64::MAX
    } else {
        (1u64 << (2 * k)) - 1
    };
    let kshift = 2 * (k - 1);
    let smask = (1u64 << (2 * s)) - 1;
    let sshift = 2 * (s - 1);
    let win = k - s + 1; // s-mers per k-mer

    let (mut kf, mut kr) = (0u64, 0u64);
    let (mut sf, mut sr) = (0u64, 0u64);
    let mut have = 0usize;
    // Canonical s-mer hashes of the current k-mer's window (front = first, back = last).
    let mut ring: std::collections::VecDeque<u64> = std::collections::VecDeque::with_capacity(win);
    let mut last_hash: Option<u64> = None;

    for (i, &b) in seq.iter().enumerate() {
        let Some(c) = code(b) else {
            have = 0;
            kf = 0;
            kr = 0;
            sf = 0;
            sr = 0;
            ring.clear();
            continue;
        };
        kf = ((kf << 2) | c) & kmask;
        kr = (kr >> 2) | ((3 - c) << kshift);
        sf = ((sf << 2) | c) & smask;
        sr = (sr >> 2) | ((3 - c) << sshift);
        have += 1;
        if have >= s {
            let scanon = if sf <= sr { sf } else { sr };
            ring.push_back(hash64(scanon));
            if ring.len() > win {
                ring.pop_front();
            }
        }
        if have >= k {
            // `ring` holds exactly this k-mer's `win` s-mer hashes.
            let minh = *ring.iter().min().expect("window non-empty");
            let closed =
                *ring.front().expect("front") == minh || *ring.back().expect("back") == minh;
            if closed {
                let pos = (i + 1 - k) as u32;
                let (canon, strand) = if kf <= kr { (kf, false) } else { (kr, true) };
                let h = hash64(canon);
                if last_hash != Some(h) {
                    out.push(Minimizer {
                        hash: h,
                        pos,
                        strand,
                    });
                    last_hash = Some(h);
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    use fqxv_dna::revcomp_acgt as revcomp;

    #[test]
    fn empty_and_short() {
        assert!(minimizers(b"", 5, 5).is_empty());
        assert!(minimizers(b"ACGT", 5, 5).is_empty(), "seq shorter than k");
    }

    #[test]
    fn all_a_kmer_does_not_hash_to_zero() {
        // The all-A k-mer packs to 0. Without splitmix64's golden-ratio
        // increment the finalizer maps 0 -> 0 — the minimum possible hash — so
        // every homopolymer run would win every window it touched and poly-A
        // would dominate minimizer selection.
        assert_ne!(hash64(0), 0, "hash64(0) must not be 0");
        // And it must not be pathologically small either.
        assert!(
            hash64(0) > u64::MAX / 1000,
            "hash64(0) is suspiciously small"
        );
    }

    #[test]
    fn homopolymer_run_collapses_to_one_minimizer() {
        // Every k-mer in a poly-A run is identical, so the whole run is one
        // hash and must collapse to a single emission — not one per window.
        let m = minimizers(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 5, 11);
        assert_eq!(
            m.len(),
            1,
            "a homopolymer run must emit exactly one minimizer"
        );
    }

    #[test]
    fn n_runs_produce_no_minimizers() {
        assert!(
            minimizers(b"NNNNNNNNNNNNNNNNNNNN", 5, 5).is_empty(),
            "an N run must not yield minimizers"
        );
    }

    #[test]
    fn density_is_about_two_over_w_plus_one() {
        // A pseudo-random sequence should select ~2/(w+1) of positions.
        let mut s = Vec::new();
        let mut x: u32 = 7;
        for _ in 0..200_000 {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            s.push(b"ACGT"[(x >> 16) as usize % 4]);
        }
        for w in [5usize, 10, 19] {
            let m = minimizers(&s, w, 15);
            let density = m.len() as f64 / s.len() as f64;
            let expect = 2.0 / (w as f64 + 1.0);
            assert!(
                (density - expect).abs() < expect * 0.25,
                "w={w}: density {density:.4} vs expected ~{expect:.4}"
            );
        }
    }

    #[test]
    fn positions_are_ascending_and_in_range() {
        let s = b"ACGTTGCAACGTTGCAGGCCATATCGCGATCGGATCAGCTAGCTAGCATCGA";
        let m = minimizers(s, 5, 7);
        assert!(!m.is_empty());
        for pair in m.windows(2) {
            assert!(pair[0].pos < pair[1].pos, "positions must strictly ascend");
        }
        assert!(m.iter().all(|x| (x.pos as usize) + 7 <= s.len()));
    }

    #[test]
    fn syncmer_empty_and_short() {
        assert!(syncmers(b"", 5, 2).is_empty());
        assert!(syncmers(b"ACGT", 5, 2).is_empty(), "seq shorter than k");
    }

    #[test]
    fn syncmer_homopolymer_collapses_to_one() {
        // Every k-mer in a poly-A run is identical (min s-mer trivially at both
        // ends → every k-mer is a closed syncmer), so the run must still collapse
        // to a single emission by the same-hash dedup.
        let m = syncmers(b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", 11, 5);
        assert_eq!(
            m.len(),
            1,
            "a homopolymer run must emit exactly one syncmer"
        );
    }

    #[test]
    fn syncmer_n_runs_produce_none() {
        assert!(
            syncmers(b"NNNNNNNNNNNNNNNNNNNN", 11, 5).is_empty(),
            "an N run must not yield syncmers"
        );
    }

    #[test]
    fn syncmer_positions_are_ascending_and_in_range() {
        let s = b"ACGTTGCAACGTTGCAGGCCATATCGCGATCGGATCAGCTAGCTAGCATCGA";
        let m = syncmers(s, 11, 5);
        assert!(!m.is_empty());
        for pair in m.windows(2) {
            assert!(pair[0].pos < pair[1].pos, "positions must strictly ascend");
        }
        assert!(m.iter().all(|x| (x.pos as usize) + 11 <= s.len()));
    }

    #[test]
    fn syncmer_density_is_about_two_over_k_minus_s_plus_one() {
        // Closed syncmers select a k-mer when its minimal s-mer lands on either
        // endpoint: 2 of the `k - s + 1` positions, so density ~= 2/(k-s+1). With
        // s = k - w this equals the minimizer's 2/(w+1) at the same w.
        let mut s = Vec::new();
        let mut x: u32 = 7;
        for _ in 0..200_000 {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            s.push(b"ACGT"[(x >> 16) as usize % 4]);
        }
        let k = 15;
        for smer in [5usize, 8, 10] {
            let m = syncmers(&s, k, smer);
            let density = m.len() as f64 / s.len() as f64;
            let expect = 2.0 / ((k - smer) as f64 + 1.0);
            assert!(
                (density - expect).abs() < expect * 0.25,
                "s={smer}: density {density:.4} vs expected ~{expect:.4}"
            );
        }
    }

    #[test]
    fn syncmer_and_minimizer_match_density_at_s_equals_k_minus_w() {
        // The load-bearing relation the wiring rests on: `s = k - w` makes the
        // closed-syncmer density equal the `w`-minimizer density, so swapping the
        // scheme leaves the index size (and the derived stride) unchanged.
        let mut s = Vec::new();
        let mut x: u32 = 99;
        for _ in 0..200_000 {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            s.push(b"ACGT"[(x >> 16) as usize % 4]);
        }
        let (w, k) = (10usize, 15usize);
        let mins = minimizers(&s, w, k).len() as f64;
        let syncs = syncmers(&s, k, k - w).len() as f64;
        assert!(
            (mins - syncs).abs() < mins * 0.1,
            "syncmer count {syncs} should track minimizer count {mins} within 10%"
        );
    }

    proptest! {
        /// The canonical property for syncmers: a sequence and its reverse
        /// complement select the SAME multiset of hashes. Closed syncmers get this
        /// from the symmetric {first, last} endpoint pair over *canonical* s-mer
        /// hashes plus the same-hash run dedup — exactly the argument `minimizers`
        /// relies on. This is what lets two reads off opposite strands overlap.
        #[test]
        fn syncmer_canonical_across_revcomp(
            bases in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 60..300),
        ) {
            let fwd = syncmers(&bases, 11, 4);
            let rev = syncmers(&revcomp(&bases), 11, 4);
            let mut a: Vec<u64> = fwd.iter().map(|m| m.hash).collect();
            let mut b: Vec<u64> = rev.iter().map(|m| m.hash).collect();
            a.sort_unstable();
            b.sort_unstable();
            prop_assert_eq!(a, b);
        }

        /// Deterministic: same input, same output, always.
        #[test]
        fn syncmer_deterministic(
            bases in proptest::collection::vec(proptest::sample::select(&b"ACGTN"[..]), 0..400),
        ) {
            prop_assert_eq!(syncmers(&bases, 9, 4), syncmers(&bases, 9, 4));
        }

        /// Never panics on arbitrary bytes.
        #[test]
        fn syncmer_arbitrary_bytes_never_panic(bases in proptest::collection::vec(any::<u8>(), 0..400)) {
            let _ = syncmers(&bases, 11, 5);
        }

        /// Every emitted `hash` is the canonical hash of the k-mer actually at
        /// `pos` — the anchor identity `index`/`chain` depend on.
        #[test]
        fn syncmer_hash_matches_kmer_at_pos(
            bases in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 20..200),
        ) {
            for m in syncmers(&bases, 11, 5) {
                let kmer = &bases[m.pos as usize..m.pos as usize + 11];
                let (mut f, mut r) = (0u64, 0u64);
                for &b in kmer {
                    let c = code(b).unwrap();
                    f = (f << 2) | c;
                    r = (r >> 2) | ((3 - c) << (2 * (11 - 1)));
                }
                let canon = if f <= r { f } else { r };
                prop_assert_eq!(m.hash, hash64(canon));
            }
        }

        /// The canonical property: a sequence and its reverse complement select
        /// the SAME multiset of minimizer hashes. This is what lets two reads
        /// off opposite strands of one locus find each other.
        #[test]
        fn canonical_across_revcomp(
            bases in proptest::collection::vec(proptest::sample::select(&b"ACGT"[..]), 60..300),
        ) {
            let fwd = minimizers(&bases, 5, 11);
            let rev = minimizers(&revcomp(&bases), 5, 11);
            let mut a: Vec<u64> = fwd.iter().map(|m| m.hash).collect();
            let mut b: Vec<u64> = rev.iter().map(|m| m.hash).collect();
            a.sort_unstable();
            b.sort_unstable();
            prop_assert_eq!(a, b);
        }

        /// Deterministic: same input, same output, always.
        #[test]
        fn deterministic(
            bases in proptest::collection::vec(proptest::sample::select(&b"ACGTN"[..]), 0..400),
        ) {
            prop_assert_eq!(minimizers(&bases, 7, 9), minimizers(&bases, 7, 9));
        }

        /// Never panics on arbitrary bytes.
        #[test]
        fn arbitrary_bytes_never_panic(bases in proptest::collection::vec(any::<u8>(), 0..400)) {
            let _ = minimizers(&bases, 5, 11);
        }
    }
}
