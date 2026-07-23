//! A small deterministic LSD radix sort over packed `u128` keys.
//!
//! The overlap search sorts a query's whole anchor set once per read
//! ([`find_overlaps`](crate::find_overlaps)); at high coverage that set is large
//! and the comparison sort showed up as ~7% of ONT compress. Every anchor's sort
//! key — `(target, strand, tpos, qpos)` — is an integer tuple that packs
//! losslessly into a `u128` whose **ascending** order is exactly the tuple's, so
//! the sort is a radix sort, not a comparison sort.
//!
//! ## Byte-identical
//!
//! LSD radix over the full key is a total order identical to `sort_unstable` on
//! the same keys (the keys are distinct-or-equal integers; equal keys are
//! interchangeable, so stability is irrelevant to the *set* of orders that are
//! correct, and the caller's downstream `dedup` collapses equals either way).
//! The result is therefore byte-identical to the comparison sort it replaces —
//! this is a speed lever, not a heuristic change. `radix_matches_comparison_sort`
//! pins it.
//!
//! minimap2 radix-sorts throughout (`radix_sort_128x`); this is the same idea,
//! written clean-room for one packed-`u128` key.

use std::cell::RefCell;

thread_local! {
    /// Double-buffer scratch reused across calls, so a per-read sort does not
    /// allocate. Cleared and resized per call — never carries state between
    /// sorts, so the result stays a pure function of the input.
    static RADIX_SCRATCH: RefCell<Vec<u128>> = const { RefCell::new(Vec::new()) };
}

/// Sort `v` ascending with an LSD radix sort over its bytes.
///
/// Only the bytes up to the most significant set byte of the maximum element are
/// swept, so a key set that never sets the high bits costs proportionally fewer
/// passes. Small inputs fall back to `sort_unstable`, where the radix passes'
/// fixed overhead would not pay for itself; both produce the identical ascending
/// order, so the threshold is a pure performance knob.
pub(crate) fn radix_sort_u128(v: &mut [u128]) {
    let n = v.len();
    // Below this the O(passes·n) radix loses to a comparison sort's cache
    // behaviour; both give the same order, so switching on size is safe.
    if n < 128 {
        v.sort_unstable();
        return;
    }

    let max = v.iter().copied().max().unwrap_or(0);
    if max == 0 {
        return; // all equal (all zero) — already sorted.
    }
    let nbytes = ((128 - max.leading_zeros()) as usize).div_ceil(8);

    RADIX_SCRATCH.with(|cell| {
        let tmp = &mut *cell.borrow_mut();
        tmp.clear();
        tmp.resize(n, 0u128);

        {
            let mut src: &mut [u128] = v;
            let mut dst: &mut [u128] = tmp.as_mut_slice();
            let mut counts = [0usize; 256];
            for byte in 0..nbytes {
                let shift = byte * 8;
                counts.fill(0);
                for &x in src.iter() {
                    counts[((x >> shift) & 0xff) as usize] += 1;
                }
                // Exclusive prefix sums -> the start offset of each bucket.
                let mut sum = 0usize;
                for c in counts.iter_mut() {
                    let cnt = *c;
                    *c = sum;
                    sum += cnt;
                }
                for &x in src.iter() {
                    let d = ((x >> shift) & 0xff) as usize;
                    dst[counts[d]] = x;
                    counts[d] += 1;
                }
                std::mem::swap(&mut src, &mut dst);
            }
        }

        // After `nbytes` passes+swaps the freshest data sits in `v` when `nbytes`
        // is even and in `tmp` when it is odd; copy back only in the odd case.
        if nbytes % 2 == 1 {
            v.copy_from_slice(tmp);
        }
    });
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        /// The whole point: the radix sort's order is byte-for-byte the
        /// comparison sort's. If it ever were not, the anchor buckets would
        /// regroup, chains would move, and the archive would change — the exact
        /// failure this lever must not introduce. Force past the small-input
        /// fallback with a length floor so the radix path itself is exercised.
        #[test]
        fn radix_matches_comparison_sort(mut v in proptest::collection::vec(any::<u128>(), 128..1000)) {
            let mut expect = v.clone();
            expect.sort_unstable();
            radix_sort_u128(&mut v);
            prop_assert_eq!(v, expect);
        }

        /// Keys shaped like the real packed anchor key `(target, strand, tpos,
        /// qpos)` — most of the high bits clear, so the pass-count shortcut and
        /// the odd/even copy-back both fire on realistic data.
        #[test]
        fn radix_matches_on_packed_shaped_keys(
            triples in proptest::collection::vec(
                (0u32..4000, any::<bool>(), 0u32..60000, 0u32..60000),
                0..1500,
            ),
        ) {
            let mut v: Vec<u128> = triples
                .iter()
                .map(|&(t, s, tp, qp)| {
                    ((t as u128) << 65)
                        | ((s as u128) << 64)
                        | ((tp as u128) << 32)
                        | (qp as u128)
                })
                .collect();
            let mut expect = v.clone();
            expect.sort_unstable();
            radix_sort_u128(&mut v);
            prop_assert_eq!(v, expect);
        }
    }

    #[test]
    fn empty_and_singleton_are_noops() {
        let mut e: Vec<u128> = vec![];
        radix_sort_u128(&mut e);
        assert!(e.is_empty());
        let mut one = vec![42u128];
        radix_sort_u128(&mut one);
        assert_eq!(one, vec![42u128]);
    }

    #[test]
    fn all_equal_is_stable_and_sorted() {
        let mut v = vec![7u128; 500];
        radix_sort_u128(&mut v);
        assert_eq!(v, vec![7u128; 500]);
    }
}
