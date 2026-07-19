//! Overlap-merge of assembled contigs into a smaller global reference.

use super::*;
use rayon::prelude::*;

/// k-mer length for detecting contig overlaps (matches the assembly minimizer).
pub(crate) const MERGE_K: usize = RESCUE_K;
/// Shortest contig-contig overlap worth merging.
pub(crate) const MIN_MERGE_OVL: usize = 24;
/// Index each contig's first `MERGE_PREFIX` bases; a successor's start must land
/// within here for the overlap to be found. Bounds the index to ~one short-read
/// worth of prefix per contig.
pub(crate) const MERGE_PREFIX: usize = 64;
/// Probe each contig's last `MERGE_SUFFIX` bases for overlaps into a successor.
pub(crate) const MERGE_SUFFIX: usize = 220;
/// Cap the candidates kept per k-mer so a repetitive k-mer can't blow up cost.
pub(crate) const MERGE_FANOUT: usize = 16;

/// Union-find root with path halving.
pub(crate) fn uf_find(parent: &mut [u32], mut x: u32) -> u32 {
    while parent[x as usize] != x {
        parent[x as usize] = parent[parent[x as usize] as usize];
        x = parent[x as usize];
    }
    x
}

/// Tunable thresholds for [`merge_reference_with`]. [`MergeConfig::default`]
/// reproduces [`merge_reference`] byte-for-byte, so sweeping these is a pure
/// encoder-side experiment — the decoder never sees the reference shape.
#[derive(Debug, Clone, Copy)]
pub struct MergeConfig {
    /// Shortest contig-contig overlap worth merging.
    pub min_ovl: usize,
    /// Index each contig's first `prefix` bases as successor entry points.
    pub prefix: usize,
    /// Probe each contig's last `suffix` bases for overlaps into a successor.
    pub suffix: usize,
    /// Cap candidates kept per k-mer so a repetitive k-mer can't blow up cost.
    pub fanout: usize,
    /// Mismatch budget for an overlap is `overlap / mism_div` (larger = stricter).
    pub mism_div: usize,
}

impl Default for MergeConfig {
    fn default() -> Self {
        Self {
            min_ovl: MIN_MERGE_OVL,
            prefix: MERGE_PREFIX,
            suffix: MERGE_SUFFIX,
            fanout: MERGE_FANOUT,
            mism_div: 8,
        }
    }
}

/// Overlap-merge with default thresholds ([`MergeConfig::default`]). See
/// [`merge_reference_with`] for the full semantics.
#[must_use]
pub fn merge_reference(
    reads: &[&[u8]],
    reference: &GlobalReference,
    places: &[Place4],
) -> (GlobalReference, Vec<Place4>) {
    merge_reference_with(reads, reference, places, MergeConfig::default())
}

/// Number of contig chunks the merge k-mer index is built over in parallel.
/// The combined index is invariant to this (chunks are combined in contig order),
/// so it affects only parallelism, not the output.
pub(crate) const MERGE_INDEX_CHUNKS: usize = 64;
/// The prefix-k-mer index is sharded by k-mer so the per-shard combines run in
/// parallel (the serial single-map combine of ~26 M keys was the merge's biggest
/// cost). Lookups route through [`merge_shard`], so this is transparent to callers.
pub(crate) const MERGE_SHARD_BITS: u32 = 6;
pub(crate) const MERGE_SHARDS: usize = 1 << MERGE_SHARD_BITS;

/// A prefix-k-mer index: `MERGE_SHARDS` maps, keyed by k-mer, each holding up to
/// `fanout` `(contig, pos)` entries. Look up a k-mer with
/// `index[merge_shard(kmer)].get(&kmer)`.
pub(crate) type MergeIndex = Vec<IntMap<u64, Vec<(u32, u32)>>>;

/// Route a k-mer to its shard via the top bits of a multiplicative mix — the low
/// bits of a 2-bit-packed k-mer are just its last base, so hash for a uniform split.
#[inline]
pub(crate) fn merge_shard(kmer: u64) -> usize {
    (kmer.wrapping_mul(0x9E37_79B9_7F4A_7C15) >> (64 - MERGE_SHARD_BITS)) as usize
}

/// Build the prefix-k-mer index for [`merge_reference_with`], fully in parallel.
/// Each contig chunk builds one fan-out-capped partial map PER SHARD (parallel
/// over chunks); then each shard is combined across chunks in contig order and
/// re-capped (parallel over shards). Both passes are parallel, and the per-shard
/// combine keeps the exact same first-N entries per k-mer as a serial build — so
/// the index is byte-for-byte independent of chunk/shard/thread count.
pub(crate) fn build_merge_index(contigs: &[&[u8]], prefix: usize, fanout: usize) -> MergeIndex {
    let nc = contigs.len();
    let chunk = nc.div_ceil(MERGE_INDEX_CHUNKS.clamp(1, nc.max(1))).max(1);
    let partials: Vec<MergeIndex> = (0..nc)
        .step_by(chunk)
        .collect::<Vec<_>>()
        .par_iter()
        .map(|&start| {
            let end = (start + chunk).min(nc);
            const MASK: u64 = (1u64 << (2 * MERGE_K)) - 1;
            let mut shards: MergeIndex = (0..MERGE_SHARDS).map(|_| IntMap::default()).collect();
            for ci in start..end {
                let c = contigs[ci];
                let hi = c.len().min(prefix);
                // Roll the k-mer instead of recomputing all MERGE_K bases per
                // start (O(bases) not O(bases * k)); a non-ACGT base resets the
                // run, so keys and insert order match the per-start recompute.
                let (mut v, mut run) = (0u64, 0usize);
                for s in 0..hi {
                    let cb = code_fold(c[s]);
                    if cb >= 4 {
                        run = 0;
                        continue;
                    }
                    v = ((v << 2) | u64::from(cb)) & MASK;
                    run += 1;
                    if run >= MERGE_K {
                        let e = shards[merge_shard(v)].entry(v).or_default();
                        if e.len() < fanout {
                            e.push((ci as u32, (s + 1 - MERGE_K) as u32));
                        }
                    }
                }
            }
            shards
        })
        .collect();
    // Combine per shard IN PARALLEL: shard `sh` merges every chunk's shard-`sh`
    // partial in chunk (contig) order, re-capping to `fanout`.
    (0..MERGE_SHARDS)
        .into_par_iter()
        .map(|sh| {
            let mut m: IntMap<u64, Vec<(u32, u32)>> = IntMap::default();
            for part in &partials {
                for (&code, list) in &part[sh] {
                    let e = m.entry(code).or_default();
                    for &item in list {
                        if e.len() < fanout {
                            e.push(item);
                        } else {
                            break;
                        }
                    }
                }
            }
            m
        })
        .collect()
}

/// Overlap-merge a greedy reference (see the module note): returns a new
/// `(reference, placements)` with fewer, longer contigs, usable by
/// [`encode_global_block`] unchanged. After chaining, the merged consensus is
/// RE-VOTED from every read at its remapped position, so overlap columns reflect
/// all contributing reads (not just the earliest contig's bytes) — that keeps the
/// per-read mismatch cost down. `reads` are the clustered, oriented reads that
/// produced `places` (read `i` has placement `places[i]`). Purely additive
/// refinement — never splits a contig, so every read keeps a valid placement.
/// `cfg` tunes the overlap-search thresholds ([`MergeConfig`]).
#[must_use]
pub fn merge_reference_with(
    reads: &[&[u8]],
    reference: &GlobalReference,
    places: &[Place4],
    cfg: MergeConfig,
) -> (GlobalReference, Vec<Place4>) {
    let nc = reference.n_contigs();
    if nc < 2 {
        return (reference.clone(), places.to_vec());
    }
    let contigs: Vec<&[u8]> = (0..nc).map(|c| reference.contig(c)).collect();

    // 1. Index each contig's PREFIX k-mers -> [(contig, pos)] (capped fan-out).
    // Built over contig CHUNKS in parallel and combined in contig order, so the
    // fan-out cap keeps the same first-N entries as a serial build — the combined
    // index is independent of the chunk count, hence of the thread count.
    let index = build_merge_index(&contigs, cfg.prefix, cfg.fanout);

    // 2. For each contig A, find its best successor B: A's suffix overlaps B's
    //    prefix (B starts at offset `s` inside A, overlap = A.len - s reaches A's
    //    end and matches B[0..overlap] within a small mismatch budget). Prefer the
    //    longest overlap, then fewest mismatches, then smallest ids (determinism).
    // Each contig's best successor depends only on the immutable `contigs` and
    // `index`, so compute them in parallel — this is the merge's hottest loop
    // (per-contig suffix probing + mismatch scans). `best_key` is a total order
    // ((MAX-ovl, mism, bi, s) minimised), so the winner — and the whole result —
    // is independent of thread count. `succ[ai] = (contig B, shift s)`.
    let succ: Vec<Option<(u32, u32)>> = (0..nc)
        .into_par_iter()
        .map(|ai| {
            let a = contigs[ai];
            if a.len() < MERGE_K {
                return None;
            }
            let lo = a.len().saturating_sub(cfg.suffix);
            let mut best_key = (usize::MAX, usize::MAX, usize::MAX, usize::MAX);
            let mut best: Option<(u32, u32)> = None;
            let mut pos_a = lo;
            while pos_a + MERGE_K <= a.len() {
                if let Some(code) = kmer_at(a, pos_a, MERGE_K)
                    && let Some(list) = index[merge_shard(code)].get(&code)
                {
                    for &(bi_u, pos_b_u) in list {
                        let bi = bi_u as usize;
                        if bi == ai {
                            continue;
                        }
                        let pos_b = pos_b_u as usize;
                        if pos_a < pos_b {
                            continue;
                        }
                        let s = pos_a - pos_b;
                        if s == 0 || s >= a.len() {
                            continue;
                        }
                        let ovl = a.len() - s;
                        let b = contigs[bi];
                        if ovl < cfg.min_ovl || ovl > b.len() {
                            continue;
                        }
                        let budget = ovl / cfg.mism_div;
                        let mut mism = 0usize;
                        for t in 0..ovl {
                            if a[s + t] != b[t] {
                                mism += 1;
                                if mism > budget {
                                    break;
                                }
                            }
                        }
                        if mism > budget {
                            continue;
                        }
                        let key = (usize::MAX - ovl, mism, bi, s);
                        if key < best_key {
                            best_key = key;
                            best = Some((bi as u32, s as u32));
                        }
                    }
                }
                pos_a += 1;
            }
            best
        })
        .collect();

    // 3. Resolve successor edges into simple chains: each contig gets at most one
    //    successor and one predecessor, no cycles (union-find). Deterministic:
    //    accept edges in contig order.
    let mut parent: Vec<u32> = (0..nc as u32).collect();
    let mut pred_taken = vec![false; nc];
    let mut chosen: Vec<Option<(u32, u32)>> = vec![None; nc];
    for ai in 0..nc {
        if let Some((bi, s)) = succ[ai] {
            let b = bi as usize;
            if pred_taken[b] {
                continue;
            }
            if uf_find(&mut parent, ai as u32) == uf_find(&mut parent, bi) {
                continue; // would close a cycle
            }
            chosen[ai] = Some((bi, s));
            pred_taken[b] = true;
            let ra = uf_find(&mut parent, ai as u32);
            let rb = uf_find(&mut parent, bi);
            parent[ra as usize] = rb;
        }
    }

    // 4. Walk each chain head (no predecessor) into a super-contig, recording each
    //    original contig's (super id, offset). Overlap bytes come from the earlier
    //    contig; only a successor's non-overlapping tail is appended.
    let mut super_id = vec![u32::MAX; nc];
    let mut super_off = vec![0u32; nc];
    let mut new_seq: Vec<u8> = Vec::with_capacity(reference.total_bases());
    let mut new_offs: Vec<usize> = vec![0];
    let mut sid = 0u32;
    for head in 0..nc {
        if pred_taken[head] {
            continue;
        }
        let super_start = new_seq.len();
        new_seq.extend_from_slice(contigs[head]);
        super_id[head] = sid;
        super_off[head] = 0;
        let mut cur = head;
        let mut base = 0usize; // super-offset of `cur` relative to super_start
        while let Some((bi, s)) = chosen[cur] {
            let bi = bi as usize;
            let bbase = base + s as usize;
            super_id[bi] = sid;
            super_off[bi] = bbase as u32;
            let cur_super_len = new_seq.len() - super_start;
            let b = contigs[bi];
            if bbase + b.len() > cur_super_len {
                let new_from = cur_super_len - bbase; // first novel base of B
                new_seq.extend_from_slice(&b[new_from..]);
            }
            base = bbase;
            cur = bi;
        }
        new_offs.push(new_seq.len());
        sid += 1;
    }

    // Remap each read onto its merged super-contig.
    let new_places: Vec<Place4> = places
        .iter()
        .map(|p| {
            let oc = p.ci as usize;
            Place4 {
                ci: super_id[oc],
                off: super_off[oc] + p.off,
            }
        })
        .collect();

    // Re-vote the merged consensus: fold every read into its remapped position and
    // take the per-column plurality (ties to the lowest base, matching the greedy
    // assembler). Overlap columns now reflect all reads, so the reads mismatch the
    // reference less — recovering most of the block-byte cost the layout-only merge
    // would otherwise add. Columns with no ACGT vote keep their laid-down byte
    // (preserving non-ACGT reference content).
    // Scatter votes in PARALLEL. The fold is O(all bases) — previously the largest
    // serial loop in the merge. Increments into a flat atomic vote array
    // (`[pos*4 + base]`) commute, so the final per-column counts — hence the
    // plurality byte — are identical regardless of thread interleaving: the output
    // stays byte-for-byte independent of thread count. Low contention: increments
    // spread over millions of columns.
    use std::sync::atomic::{AtomicU32, Ordering};
    let votes: Vec<AtomicU32> = (0..new_seq.len() * 4).map(|_| AtomicU32::new(0)).collect();
    reads
        .par_iter()
        .zip(new_places.par_iter())
        .for_each(|(r, pl)| {
            let start = new_offs[pl.ci as usize] + pl.off as usize;
            for (j, &byte) in r.iter().enumerate() {
                let c = code_fold(byte);
                if c < 4 {
                    votes[(start + j) * 4 + c as usize].fetch_add(1, Ordering::Relaxed);
                }
            }
        });
    // Per-column plurality is independent per position, so resolve in parallel.
    // Deterministic: each output byte is a pure function of that column's counts
    // (ties to the lowest base via `Reverse(k)`).
    new_seq.par_iter_mut().enumerate().for_each(|(i, byte)| {
        let v = [
            votes[i * 4].load(Ordering::Relaxed),
            votes[i * 4 + 1].load(Ordering::Relaxed),
            votes[i * 4 + 2].load(Ordering::Relaxed),
            votes[i * 4 + 3].load(Ordering::Relaxed),
        ];
        if v.iter().any(|&x| x > 0) {
            let best = (0..4)
                .max_by_key(|&k| (v[k], std::cmp::Reverse(k)))
                .unwrap();
            *byte = b"ACGT"[best];
        }
    });

    let merged = GlobalReference {
        seq: new_seq,
        offs: new_offs,
    };
    (merged, new_places)
}
