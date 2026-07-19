//! Minimizer index over a read set: hash -> the reads and offsets carrying it.
//!
//! Built once per block, then probed per read to collect anchors for chaining.
//!
//! ## Determinism
//!
//! Output must not depend on thread count (a workspace invariant). Minimizers
//! are computed with `par_iter().map().collect()`, which preserves read order
//! regardless of scheduling; the flattened occurrences are then sorted by a
//! **total** order `(hash, read, pos)`, so the layout is a pure function of the
//! input. Nothing is ever iterated out of a hash map — there is no hash map.

use rayon::prelude::*;

use crate::{Error, Result, Sketch};

/// One minimizer occurrence within the read set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Occ {
    /// Hash of the canonical k-mer.
    pub hash: u64,
    /// Index of the read carrying it.
    pub read: u32,
    /// Offset of the k-mer within that read.
    pub pos: u32,
    /// True when the reverse-complement orientation was canonical here.
    pub strand: bool,
}

/// How aggressively to discard repetitive minimizers.
///
/// A minimizer at a single-copy locus occurs about once per covering read, so at
/// ~300x coverage ~300 occurrences is *normal*, not repetitive — an absolute
/// threshold would therefore be a coverage-dependent guess. Instead drop the
/// most frequent `drop_top_frac` of distinct minimizers, which adapts to
/// coverage on its own. minimap2 uses the same shape (its default discards the
/// top ~0.02% of distinct minimizers).
///
/// Repetitive minimizers are worth dropping because they are quadratic: a
/// minimizer occurring `n` times contributes `n^2` candidate anchor pairs while
/// carrying almost no positional information.
///
/// Note this targets motifs recurring at **scattered** loci, not tandem repeats.
/// A tandem repeat's k-mers cycle with its period, so consecutive windows select
/// the same hash and [`minimizers`](crate::minimizers)' same-hash dedup already
/// collapses the run to a few emissions — it never reaches the top of the
/// frequency distribution in the first place.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Repeat {
    /// Fraction of the most frequent *distinct* minimizers to discard.
    pub drop_top_frac: f64,
}

impl Default for Repeat {
    fn default() -> Self {
        Self {
            drop_top_frac: 0.0002,
        }
    }
}

/// A minimizer index over a set of reads.
#[derive(Debug)]
pub struct Index {
    sketch: Sketch,
    /// All kept occurrences, sorted by `(hash, read, pos)`.
    occs: Vec<Occ>,
    /// Per-read lengths — needed to map a reverse-complement hit into the
    /// query's coordinate frame.
    lens: Vec<u32>,
    /// Occurrence count of the most frequent minimizer kept.
    max_kept: usize,
    /// Distinct minimizers discarded as repetitive.
    dropped: usize,
}

impl Index {
    /// Build an index over `n` reads whose sequences are concatenated in `seq`
    /// with lengths `lens`.
    ///
    /// # Errors
    ///
    /// [`Error::LengthMismatch`] if `lens` does not sum to `seq.len()`.
    pub fn build(lens: &[u32], seq: &[u8], sketch: Sketch, repeat: Repeat) -> Result<Self> {
        let total: usize = lens.iter().map(|&l| l as usize).sum();
        if total != seq.len() {
            return Err(Error::LengthMismatch {
                lens: total,
                seq: seq.len(),
            });
        }

        // Per-read byte offsets, so each read can be sketched independently.
        let mut offs = Vec::with_capacity(lens.len() + 1);
        let mut acc = 0usize;
        for &l in lens {
            offs.push(acc);
            acc += l as usize;
        }
        offs.push(acc);

        // `par_iter().map().collect()` preserves input order regardless of how
        // rayon schedules the work, so this is thread-count invariant.
        let per_read: Vec<Vec<Occ>> = (0..lens.len())
            .into_par_iter()
            .map(|r| {
                let s = &seq[offs[r]..offs[r + 1]];
                sketch
                    .seeds(s)
                    .into_iter()
                    .map(|m| Occ {
                        hash: m.hash,
                        read: r as u32,
                        pos: m.pos,
                        strand: m.strand,
                    })
                    .collect()
            })
            .collect();

        let mut occs: Vec<Occ> = per_read.into_iter().flatten().collect();
        // Total order on the whole struct (derived Ord: hash, read, pos, strand)
        // — so ties cannot depend on sort stability or thread count.
        occs.par_sort_unstable();

        // Repeat control. Count occurrences per distinct hash, then drop the
        // most frequent `drop_top_frac`. Counts are read off the sorted runs, so
        // this never touches parallel state.
        let mut counts: Vec<(usize, u64)> = Vec::new(); // (count, hash)
        let mut i = 0usize;
        while i < occs.len() {
            let h = occs[i].hash;
            let mut j = i;
            while j < occs.len() && occs[j].hash == h {
                j += 1;
            }
            counts.push((j - i, h));
            i = j;
        }
        // Sort by count DESC, then hash ASC — a total order, so which hashes sit
        // at the cutoff is fully determined.
        counts.sort_unstable_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
        let n_drop = ((counts.len() as f64) * repeat.drop_top_frac).floor() as usize;
        let n_drop = n_drop.min(counts.len());

        let mut drop_set: Vec<u64> = counts[..n_drop].iter().map(|&(_, h)| h).collect();
        drop_set.sort_unstable();
        if !drop_set.is_empty() {
            occs.retain(|o| drop_set.binary_search(&o.hash).is_err());
        }
        let max_kept = counts.get(n_drop).map_or(0, |&(c, _)| c);

        Ok(Self {
            sketch,
            occs,
            lens: lens.to_vec(),
            max_kept,
            dropped: n_drop,
        })
    }

    /// Length of read `r`, or 0 if out of range.
    #[must_use]
    pub fn read_len(&self, r: u32) -> u32 {
        self.lens.get(r as usize).copied().unwrap_or(0)
    }

    /// Number of reads indexed.
    #[must_use]
    pub fn n_reads(&self) -> usize {
        self.lens.len()
    }

    /// The occurrences of `hash`, ascending by `(read, pos)`; empty if absent or
    /// discarded as repetitive.
    #[must_use]
    pub fn query(&self, hash: u64) -> &[Occ] {
        let lo = self.occs.partition_point(|o| o.hash < hash);
        let hi = self.occs.partition_point(|o| o.hash <= hash);
        &self.occs[lo..hi]
    }

    /// The sketch parameters this index was built with.
    #[must_use]
    pub fn sketch(&self) -> Sketch {
        self.sketch
    }

    /// Total kept occurrences.
    #[must_use]
    pub fn len(&self) -> usize {
        self.occs.len()
    }

    /// True when no occurrence was kept.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.occs.is_empty()
    }

    /// Occurrence count of the most frequent minimizer that survived repeat
    /// filtering — the effective cap that filtering produced.
    #[must_use]
    pub fn max_kept(&self) -> usize {
        self.max_kept
    }

    /// Number of distinct minimizers discarded as repetitive.
    #[must_use]
    pub fn dropped(&self) -> usize {
        self.dropped
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::SeedScheme;

    /// Deterministic pseudo-random bases.
    fn rand_seq(n: usize, seed: u32) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
                b"ACGT"[(x >> 16) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn length_mismatch_is_rejected() {
        let e = Index::build(&[10, 10], b"ACGT", Sketch::ont(), Repeat::default());
        assert!(matches!(e, Err(Error::LengthMismatch { .. })));
    }

    #[test]
    fn empty_input_builds_empty_index() {
        let idx = Index::build(&[], b"", Sketch::ont(), Repeat::default()).unwrap();
        assert!(idx.is_empty());
        assert!(idx.query(12345).is_empty());
    }

    #[test]
    fn query_returns_only_the_asked_hash() {
        let a = rand_seq(2000, 1);
        let b = rand_seq(2000, 2);
        let mut seq = a.clone();
        seq.extend_from_slice(&b);
        let idx = Index::build(
            &[a.len() as u32, b.len() as u32],
            &seq,
            Sketch::ont(),
            Repeat { drop_top_frac: 0.0 },
        )
        .unwrap();
        assert!(!idx.is_empty());
        // Every occurrence of some present hash must carry that hash, and the
        // slice must be exactly the run.
        let h = idx.occs[idx.occs.len() / 2].hash;
        let got = idx.query(h);
        assert!(!got.is_empty());
        assert!(got.iter().all(|o| o.hash == h));
        assert_eq!(got.len(), idx.occs.iter().filter(|o| o.hash == h).count());
        // Absent hash -> empty.
        assert!(idx.query(u64::MAX).is_empty());
    }

    #[test]
    fn identical_reads_share_every_minimizer() {
        // Two copies of the same read must co-occur under every hash: this is
        // the property the whole overlap search rests on.
        let r = rand_seq(3000, 7);
        let mut seq = r.clone();
        seq.extend_from_slice(&r);
        let idx = Index::build(
            &[r.len() as u32, r.len() as u32],
            &seq,
            Sketch::ont(),
            Repeat { drop_top_frac: 0.0 },
        )
        .unwrap();
        let single = Sketch::ont().seeds(&r);
        for m in &single {
            let occ = idx.query(m.hash);
            let reads: Vec<u32> = occ.iter().map(|o| o.read).collect();
            assert!(
                reads.contains(&0) && reads.contains(&1),
                "hash {} must occur in both copies",
                m.hash
            );
        }
    }

    #[test]
    fn build_is_deterministic() {
        let r1 = rand_seq(4000, 11);
        let r2 = rand_seq(4000, 13);
        let mut seq = r1.clone();
        seq.extend_from_slice(&r2);
        let lens = [r1.len() as u32, r2.len() as u32];
        let a = Index::build(&lens, &seq, Sketch::ont(), Repeat::default()).unwrap();
        let b = Index::build(&lens, &seq, Sketch::ont(), Repeat::default()).unwrap();
        assert_eq!(
            a.occs, b.occs,
            "index layout must be a pure function of input"
        );
        assert_eq!(a.dropped(), b.dropped());
    }

    #[test]
    fn a_tandem_repeat_is_collapsed_by_dedup_not_by_the_filter() {
        // Counter-intuitive, so pinned: a TANDEM repeat never becomes a frequent
        // minimizer. Its k-mers cycle with the repeat period, so consecutive
        // windows select the same hash and `minimizers`' same-hash dedup
        // collapses the whole run to a handful of emissions. Tandem repeats are
        // handled upstream and are NOT what the repeat filter is for — see
        // `repeat_filter_drops_the_most_frequent`.
        let mut seq = rand_seq(2000, 17);
        for _ in 0..200 {
            seq.extend_from_slice(b"ACGTACGTACGTACGTACGTACGTACGTACGT");
        }
        // Same-hash-dedup collapse of a tandem repeat is a *window-minimizer*
        // property: consecutive windows re-select the run's single minimum. A
        // syncmer is chosen per-k-mer, so a repeat's phases don't collapse the
        // same way (the repeat *filter* bounds them there instead). ONT now
        // defaults to syncmers, so pin the minimizer scheme explicitly here.
        let idx = Index::build(
            &[seq.len() as u32],
            &seq,
            Sketch {
                w: 10,
                k: 15,
                scheme: SeedScheme::Minimizer,
            },
            Repeat { drop_top_frac: 0.0 },
        )
        .unwrap();
        // 6400 bp of tandem repeat yields only a few occurrences, so the busiest
        // minimizer stays far below the repeat's 200 copies.
        assert!(
            idx.max_kept() < 50,
            "tandem repeat should collapse under dedup, but max_kept = {}",
            idx.max_kept()
        );
    }

    #[test]
    fn repeat_filter_drops_the_most_frequent() {
        // A motif recurring at SCATTERED loci is what the filter exists for: the
        // occurrences are non-consecutive, so they survive the same-hash dedup
        // and genuinely pile up (unlike a tandem run — see the test above).
        //
        // An explicit fraction, not `Repeat::default()`: the 0.02% default is
        // scaled for the millions of distinct minimizers in a real block, and on
        // a fixture this size it floors to zero dropped (correct, but not a test
        // of the mechanism).
        let motif = rand_seq(300, 99);
        let mut seq = Vec::new();
        for i in 0..40u32 {
            seq.extend_from_slice(&rand_seq(400, 1000 + i));
            seq.extend_from_slice(&motif);
        }
        let lens = [seq.len() as u32];
        let none = Index::build(&lens, &seq, Sketch::ont(), Repeat { drop_top_frac: 0.0 }).unwrap();
        let some = Index::build(
            &lens,
            &seq,
            Sketch::ont(),
            Repeat {
                drop_top_frac: 0.05,
            },
        )
        .unwrap();
        assert_eq!(none.dropped(), 0);
        // Sanity: the motif really does recur ~40x, so there is a peak to remove.
        assert!(
            none.max_kept() >= 40,
            "motif should occur ~40x, got max_kept = {}",
            none.max_kept()
        );
        assert!(some.dropped() > 0, "a scattered repeat must be filtered");
        assert!(
            some.len() < none.len(),
            "filtering must remove occurrences ({} vs {})",
            some.len(),
            none.len()
        );
        // The peak is gone: what survives is far less frequent.
        assert!(
            some.max_kept() < none.max_kept(),
            "the most frequent minimizer must be gone ({} vs {})",
            some.max_kept(),
            none.max_kept()
        );
    }

    #[test]
    fn default_frac_drops_nothing_on_a_small_input() {
        // The 0.02% default is scaled for a real block's millions of distinct
        // minimizers; on a small input it floors to zero. Pin that so the
        // rounding behaviour is a decision, not an accident.
        let s = rand_seq(11_000, 23);
        let idx = Index::build(&[s.len() as u32], &s, Sketch::ont(), Repeat::default()).unwrap();
        assert_eq!(idx.dropped(), 0);
    }

    #[test]
    fn drop_frac_zero_keeps_everything() {
        let s = rand_seq(3000, 19);
        let idx = Index::build(
            &[s.len() as u32],
            &s,
            Sketch::ont(),
            Repeat { drop_top_frac: 0.0 },
        )
        .unwrap();
        assert_eq!(idx.dropped(), 0);
        assert_eq!(idx.len(), Sketch::ont().seeds(&s).len());
    }
}
