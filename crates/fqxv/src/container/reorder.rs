//! Globally-clustered (SPRING-style) reorder layout: encode and decode.

use super::*;
use rayon::prelude::*;
use tracing::info;

/// Mean read length (bp) above which globally-clustered reorder is skipped. On
/// long-read data (nanopore/PacBio, ~10-14 kb reads) reorder yields no ratio
/// gain — the non-reorder deep-context path is actually *smaller* — for roughly
/// 10x the compress time and 6x the peak memory. Illumina reads are <= ~300 bp,
/// so this threshold cleanly separates the two regimes. See [`is_long_read`].
pub(crate) const REORDER_MAX_MEAN_LEN: u64 = 500;
/// Minimizer length for clustering in reorder mode.
pub(crate) const REORDER_K: usize = 15;
/// Reads per block in whole-file (global-cluster) reorder mode. Moderate, so the
/// sequence and name/quality blocks fan out across cores. Clustering is global,
/// so block size no longer trades against ratio — only parallelism and the
/// per-block model reset (cheap, since clustered duplicates collapse to MATCH).
pub(crate) const REORDER_BLOCK_READS: usize = 1 << 18;
/// Windows the global assembly is split into for parallel building. Fixed (never
/// derived from the thread count) so the reference is byte-identical regardless of
/// `--threads`. Enough to saturate a many-core node; cross-window dedup lost to
/// the split is recovered by the overlap-merge that follows.
pub(crate) const ASSEMBLY_WINDOWS: usize = 8;

/// Buffer every record of a single (possibly interleaved) FASTQ stream into one
/// [`RawBlock`], preserving input order. Used by the reorder path, which needs the
/// whole file resident before it can cluster globally.
pub(crate) fn buffer_records(buf: &[u8]) -> Result<RawBlock> {
    let mut all = RawBlock::default();
    let mut r: &[u8] = buf;
    let (mut def, mut seq, mut qual) = (Vec::new(), Vec::new(), Vec::new());
    while read_raw_record(&mut r, &mut def, &mut seq, &mut qual)? {
        all.push_raw(&def, &seq, &qual);
    }
    Ok(all)
}

/// Mean of a read-length sample (the first `sample` reads is plenty to tell
/// ~14 kb long reads from ~150 bp short reads). Zero for an empty block.
pub(crate) fn mean_read_len(lens: &[u32]) -> u64 {
    if lens.is_empty() {
        return 0;
    }
    let sample = lens.len().min(256);
    let sum: u64 = lens[..sample].iter().map(|&l| u64::from(l)).sum();
    sum / sample as u64
}

/// True when the data is long-read (mean length over [`REORDER_MAX_MEAN_LEN`]),
/// for which the globally-clustered reorder layout is skipped in favour of the
/// non-reorder deep-context path.
pub(crate) fn is_long_read(lens: &[u32]) -> bool {
    mean_read_len(lens) > REORDER_MAX_MEAN_LEN
}

/// Serialize a buffered block back to interleaved FASTQ, so a reorder-path block
/// can be handed to the non-reorder encoder ([`compress_buffered`]) unchanged.
/// Record order (hence any mate interleaving) is preserved.
pub(crate) fn serialize_block(all: &RawBlock) -> Vec<u8> {
    let mut buf = Vec::with_capacity(
        all.seq.len() + all.qual.len() + all.header_buf.len() + all.n_reads() * 4,
    );
    let mut off = 0usize;
    for i in 0..all.n_reads() {
        let l = all.lens[i] as usize;
        write_record(
            &mut buf,
            all.header(i),
            &all.seq[off..off + l],
            &all.qual[off..off + l],
        );
        off += l;
    }
    buf
}

/// Buffer a single reader and hand off to [`encode_reordered`] (single-end when
/// `group_size == 1`, or an already-interleaved stream for `group_size > 1`).
pub(crate) fn compress_reordered_whole<R: Read + Send, W: Write>(
    reader: R,
    writer: W,
    params: Params,
    group_size: u8,
) -> Result<Stats> {
    // Buffer every read; input order is preserved, so an interleaved stream stays
    // interleaved and the permutation can restore that spot order on decode.
    let mut all = RawBlock::default();
    let mut br = BufReader::new(reader);
    let (mut def, mut seq, mut qual) = (Vec::new(), Vec::new(), Vec::new());
    while read_raw_record(&mut br, &mut def, &mut seq, &mut qual)? {
        all.push_raw(&def, &seq, &qual);
    }
    encode_reordered(all, writer, params, group_size)
}

/// Globally cluster the buffered reads (SPRING-style) and write the whole-file
/// reorder archive: cluster once, then code the clustered sequence in independent
/// moderate blocks that fan out across cores (clustering is global, so block size
/// trades only against parallelism, not ratio). Two modes:
///
/// - `keep_order`: names+quality are coded in ORIGINAL order and a global
///   permutation (byte-plane rANS) restores it, so reads come back byte-exact in
///   input order.
/// - without `keep_order`: names+quality are coded in CLUSTERED order and NO
///   permutation is written, so decode emits reads in clustered order (records
///   preserved as a set) — smaller, but order is not restorable.
///
/// `group_size` is the mate interleaving (1 = single-end). When `group_size > 1`
/// the original spot order *is* the mate interleaving, so the permutation is
/// required to reconstruct it: `keep_order` is forced on regardless of
/// `params.keep_order`, making grouped reorder order-preserving.
pub(crate) fn encode_reordered<W: Write>(
    all: RawBlock,
    writer: W,
    params: Params,
    group_size: u8,
) -> Result<Stats> {
    // Long-read data: reorder buys nothing (the non-reorder deep-context path is
    // smaller on nanopore/PacBio) for ~10x the time and ~6x the memory. Fall back
    // to the non-reorder layout, keeping the requested effort level — its hashed
    // high-order sequence context is exactly what wins on long reads. This makes
    // reorder adaptive: `--order any` / `--max` still do the right thing on a
    // long-read input instead of paying a large cost for no benefit.
    if is_long_read(&all.lens) {
        info!(
            mean_len = mean_read_len(&all.lens),
            "long-read data: skipping reorder (no ratio benefit, high cost) — using non-reorder layout"
        );
        let mut p = params;
        p.reorder = false;
        return compress_buffered(&serialize_block(&all), writer, p, group_size);
    }
    let g = group_size.max(1);
    let pool = build_pool(params.threads)?;
    let n = all.n_reads();

    // Cumulative byte offsets into the concatenated seq/qual.
    let mut offs = Vec::with_capacity(n + 1);
    let mut acc = 0usize;
    for &l in &all.lens {
        offs.push(acc);
        acc += l as usize;
    }
    offs.push(acc);

    // `offs` is derived from `lens` alone, and everything below indexes `seq` and
    // `qual` with it — so if the lengths disagree with the buffers, every one of
    // those slices is out of bounds and the process ABORTS (the release profile is
    // `panic = "abort"`). A truncated FASTQ does exactly that: it is the most
    // ordinary corruption there is, and `--order any` core-dumped on it while
    // `--order preserve` rejected it cleanly, because preserve hands the whole
    // buffer to a codec that checks. Check here, once, before the first slice.
    //
    // Raised as the codecs' own errors rather than a new one, so both `--order`
    // modes fail with the identical message for the identical input — the two
    // modes disagreeing about whether a file is valid is the actual bug.
    if acc != all.seq.len() {
        return Err(fqxv_seq::Error::LengthMismatch {
            lens: acc,
            seq: all.seq.len(),
        }
        .into());
    }
    if acc != all.qual.len() {
        return Err(fqxv_fqzcomp::Error::LengthMismatch {
            lens: acc,
            quals: all.qual.len(),
        }
        .into());
    }

    // 2. One global clustering pass over every read.
    let plan = pool.install(|| fqxv_reorder::plan(&all.lens, &all.seq, REORDER_K));

    // 3. Clustered, oriented sequences (parallel copy/revcomp) + flip bitmap.
    //
    // A flat arena plus offsets, not a `Vec<Vec<u8>>`. The old shape cost one heap
    // allocation per read — 256Ki of them per block, `malloc`/`free` measuring
    // ~6.9% of the reorder encode — and it held each read behind its own pointer
    // when `RawBlock` already hands us exactly this layout for the input side.
    // Every consumer only ever wanted a `&[u8]`.
    let mut cl_offs = Vec::with_capacity(n + 1);
    let mut cl_acc = 0usize;
    for &oi in &plan.order {
        cl_offs.push(cl_acc);
        cl_acc += all.lens[oi as usize] as usize;
    }
    cl_offs.push(cl_acc);
    debug_assert_eq!(cl_acc, acc, "clustering is a permutation: same total bases");

    let mut cl_seq = vec![0u8; cl_acc];
    {
        // Hand each read its own disjoint destination up front. Splitting serially
        // is a few pointer moves and keeps the parallel fill safe — no `unsafe`,
        // no per-read allocation, and each worker writes only its own slice.
        let mut dsts: Vec<&mut [u8]> = Vec::with_capacity(n);
        let mut rest: &mut [u8] = &mut cl_seq;
        for j in 0..n {
            let (head, tail) = rest.split_at_mut(cl_offs[j + 1] - cl_offs[j]);
            dsts.push(head);
            rest = tail;
        }
        pool.install(|| {
            dsts.par_iter_mut().enumerate().for_each(|(j, dst)| {
                let oi = plan.order[j] as usize;
                let s = &all.seq[offs[oi]..offs[oi + 1]];
                if plan.flip[oi] {
                    fqxv_reorder::revcomp_into(s, dst);
                } else {
                    dst.copy_from_slice(s);
                }
            });
        });
    }
    let cl_read = |j: usize| -> &[u8] { &cl_seq[cl_offs[j]..cl_offs[j + 1]] };
    let mut flip_bits = vec![0u8; n.div_ceil(8)];
    for (j, &oi) in plan.order.iter().enumerate() {
        if plan.flip[oi as usize] {
            flip_bits[j / 8] |= 1 << (j % 8);
        }
    }

    // Minimizer anchors in clustered order (for shifted-overlap alignment).
    let cl_anchors: Vec<u32> = plan
        .order
        .iter()
        .map(|&oi| plan.anchor[oi as usize])
        .collect();

    // 4. Moderate blocks (same count for both partitions).
    let bsz = REORDER_BLOCK_READS.max(1);
    let ranges: Vec<(usize, usize)> = if n == 0 {
        vec![(0, 0)]
    } else {
        (0..n).step_by(bsz).map(|s| (s, (s + bsz).min(n))).collect()
    };

    // 5a. Sequence — clustered order, differential-coded per block, in parallel.
    //
    // Three coexisting codecs, arranged to be NEVER WORSE than the block-local
    // baseline:
    //   * v2 single-contig and v3 literal-rescue are BLOCK-LOCAL (each block is
    //     self-contained); the container keeps the smaller per block, as before.
    //   * v4 codes reads as positions on ONE frozen global reference assembled
    //     over every clustered read (SPRING-style), so the cross-block overlaps v3
    //     strands as literals collapse to a cheap (contig, offset, mismatches)
    //     back-reference — at the cost of a whole-file reference frame stored once.
    // Pass 1 builds the reference; pass 2 codes every block against it in parallel;
    // then ONE whole-file choice keeps the reference layout only when it pays:
    //   reference_frame + Σ min(v2, v3, v4)  <  Σ min(v2, v3).
    // Otherwise no reference frame is written and the archive is byte-for-byte the
    // v2/v3 layout it was before — so v4 can only ever shrink the output. Blocks
    // may mix versions freely (decode dispatches on the leading version byte).
    // `--no-rescue` (`rescue = false`) forces the fast v2-only path and skips the
    // global assembly entirely. Ties keep the lower version for determinism.
    let order = params.seq_order as usize;

    // Pass 1: assemble every clustered read into a frozen reference plus per-read
    // placements. The greedy fold is serial, so it runs over a fixed set of
    // windows IN PARALLEL (`assemble_global_windowed`) and the overlap-merge below
    // reclaims the cross-window deduplication — a near-ratio-neutral speedup of
    // what used to be the compress-time floor. Only run when v4 is a candidate
    // (the adaptive `rescue` path).
    let global = if params.rescue && n > 0 {
        let refs_all: Vec<&[u8]> = (0..n).map(cl_read).collect();
        let (reference, places) = pool.install(|| {
            fqxv_reorder::assemble_global_windowed(&refs_all, &cl_anchors, ASSEMBLY_WINDOWS)
        });
        // Overlap-merge the greedy reference: chain contigs whose suffix overlaps
        // another's prefix into fewer, longer super-contigs, store the shared
        // sequence once, and remap placements. Format-transparent (encode_global_block
        // re-derives mismatches against whatever reference it is handed), so the
        // decoder and on-disk layout are unchanged. A single default pass captures
        // ~all the gain (iterating or loosening thresholds adds <=0.1-0.2%); it is
        // strongly scale-dependent — measured -8.4% vs v3 on 4M NovaSeq (the
        // reference shrinks 26.8 -> 23.8 MB), where the plain assembly alone was
        // only -3.1%. The whole-file decision below still gates the reference
        // layout, so this can only ever shrink the archive.
        let (reference, places) =
            pool.install(|| fqxv_reorder::merge_reference(&refs_all, &reference, &places));
        // Code the reference once with our clean-room order-k model, split into
        // fixed blocks compressed in parallel (`encode_reference_frame`). This is
        // fast — it replaces the single-threaded xz pass that used to dominate
        // compress time — and adds no external-compressor dependency.
        let ref_payload = pool.install(|| encode_reference_frame(&reference, order))?;
        Some((reference, places, ref_payload))
    } else {
        None
    };

    // Pass 2: per block, the block-local best (v2/v3) and — when a reference
    // exists — the reference-inclusive best (v2/v3/v4). Both are kept until the
    // whole-file decision below picks one layout for the archive.
    struct BlockChoice {
        block_local: Vec<u8>,
        with_ref: Vec<u8>,
    }
    let choices: Vec<BlockChoice> = pool.install(|| {
        ranges
            .par_iter()
            .map(|&(s, e)| -> Result<BlockChoice> {
                let refs: Vec<&[u8]> = (s..e).map(cl_read).collect();
                let anch = &cl_anchors[s..e];
                let mut block_local = fqxv_reorder::encode_clustered(&refs, anch, order)?;
                if params.rescue {
                    let v3 = fqxv_reorder::encode_clustered_rescue(&refs, anch, order)?;
                    if v3.len() < block_local.len() {
                        block_local = v3;
                    }
                }
                let with_ref = match &global {
                    Some((reference, places, _)) => {
                        let v4 =
                            fqxv_reorder::encode_global_block(&refs, &places[s..e], reference)?;
                        if v4.len() < block_local.len() {
                            v4
                        } else {
                            block_local.clone()
                        }
                    }
                    None => Vec::new(), // unused when there is no reference
                };
                Ok(BlockChoice {
                    block_local,
                    with_ref,
                })
            })
            .collect::<Result<_>>()
    })?;

    // Whole-file decision: adopt the reference layout only if it is strictly
    // smaller than the block-local layout (reference frame included).
    let (use_reference, seq_blocks, ref_payload): (bool, Vec<Vec<u8>>, Vec<u8>) = match global {
        Some((_, _, ref_payload)) => {
            let with_ref_total =
                ref_payload.len() + choices.iter().map(|c| c.with_ref.len()).sum::<usize>();
            let block_local_total = choices.iter().map(|c| c.block_local.len()).sum::<usize>();
            if with_ref_total < block_local_total {
                let blocks = choices.into_iter().map(|c| c.with_ref).collect();
                (true, blocks, ref_payload)
            } else {
                let blocks = choices.into_iter().map(|c| c.block_local).collect();
                (false, blocks, Vec::new())
            }
        }
        None => {
            let blocks = choices.into_iter().map(|c| c.block_local).collect();
            (false, blocks, Vec::new())
        }
    };

    // 5b. Names, then the keep_order decision, then quality.
    //
    // Reorder has two layouts: keep_order codes names/quality in ORIGINAL order
    // and stores a permutation; otherwise they're coded in CLUSTERED order and no
    // permutation is stored (reads emerge clustered). For single-end input we
    // pick ADAPTIVELY: counter-style names (e.g. SRA `.N N`) delta-code to almost
    // nothing in original order, so a permutation is cheaper than a scrambled
    // name stream; random names are the reverse. Grouped input (the permutation
    // reconstructs spots) and an explicit `params.keep_order` force keep_order.

    // Names coded in ORIGINAL order, per block.
    let names_original = || -> Result<Vec<Vec<u8>>> {
        pool.install(|| {
            ranges
                .par_iter()
                .map(|&(s, e)| {
                    let headers: Vec<&[u8]> = (s..e).map(|i| all.header(i)).collect();
                    Ok(fqxv_tokenizer::encode(&headers)?)
                })
                .collect()
        })
    };
    // Names coded in CLUSTERED order, per block.
    let names_clustered = || -> Result<Vec<Vec<u8>>> {
        pool.install(|| {
            ranges
                .par_iter()
                .map(|&(s, e)| {
                    let headers: Vec<&[u8]> = plan.order[s..e]
                        .iter()
                        .map(|&oi| all.header(oi as usize))
                        .collect();
                    Ok(fqxv_tokenizer::encode(&headers)?)
                })
                .collect()
        })
    };
    // Global permutation (byte-plane split → rANS), coded with whichever order is
    // smaller (decode auto-detects). The two regimes differ: on huge inputs
    // order-1 wins (real byte-to-byte correlation, its per-context header
    // amortized); on small-to-medium inputs order-1's ~130 KB header dominates
    // and order-0 wins — picking the smaller keeps keep_order efficient at every
    // size. Perm encode is a small fraction of total, so trying both is cheap.
    let encode_perm = || -> Result<Vec<u8>> {
        let mut planes = vec![0u8; n * 4];
        for (i, &x) in plan.order.iter().enumerate() {
            planes[i] = x as u8;
            planes[n + i] = (x >> 8) as u8;
            planes[2 * n + i] = (x >> 16) as u8;
            planes[3 * n + i] = (x >> 24) as u8;
        }
        let o0 = fqxv_rans::encode(&planes, fqxv_rans::Order::Zero)?;
        let o1 = fqxv_rans::encode(&planes, fqxv_rans::Order::One)?;
        Ok(if o0.len() <= o1.len() { o0 } else { o1 })
    };

    // Discard-order (opt-in, single-end): if the names are purely positional (a
    // counter), regenerate them from a tiny template instead of coding them — no
    // name stream, no permutation. Reorder-lossy (reads are renumbered), gated on
    // `params.regenerate_names`. When the names ARE a regenerable counter we
    // reproduce them exactly; otherwise (Illumina tile/x/y, spliced counters that
    // reset, …) we still renumber — with a fresh synthetic 1..n counter — rather
    // than falling back to storing a full permutation. This is SPRING's behavior:
    // discarding read identity is exactly what the caller opted into, and it drops
    // the ~read-count-sized permutation that otherwise dominates a reordered
    // archive. Only whole-file order-preserving inputs (`keep_order`, grouped)
    // keep the permutation.
    let template = if params.regenerate_names && !(params.keep_order || g > 1) {
        let orig_names: Vec<&[u8]> = (0..n).map(|i| all.header(i)).collect();
        Some(
            fqxv_tokenizer::detect_template(&orig_names)
                .unwrap_or_else(fqxv_tokenizer::NameTemplate::renumber),
        )
    } else {
        None
    };

    let (keep_order, name_blocks, perm_c) = if template.is_some() {
        // Clustered layout, no permutation, empty (regenerated) name blocks.
        (false, vec![Vec::new(); ranges.len()], Vec::new())
    } else if params.keep_order || g > 1 {
        (true, names_original()?, encode_perm()?)
    } else {
        // Adaptive: keep_order iff original-order names + permutation beat the
        // clustered-order name stream. (Quality's order dependence is second-order
        // and ignored here.) Deterministic — sizes don't depend on thread count.
        let orig = names_original()?;
        let clustered = names_clustered()?;
        let perm = encode_perm()?;
        let keep_bytes = orig.iter().map(Vec::len).sum::<usize>() + perm.len();
        let cluster_bytes = clustered.iter().map(Vec::len).sum::<usize>();
        if keep_bytes < cluster_bytes {
            (true, orig, perm)
        } else {
            (false, clustered, Vec::new())
        }
    };

    // Quality in the chosen order: original for keep_order; otherwise clustered,
    // reversed for flipped reads so bytes line up with the reverse-complemented
    // sequence.
    let qual_blocks: Vec<Vec<u8>> = pool.install(|| {
        ranges
            .par_iter()
            .map(|&(s, e)| -> Result<Vec<u8>> {
                if keep_order {
                    Ok(fqxv_fqzcomp::encode(
                        &all.lens[s..e],
                        &all.qual[offs[s]..offs[e]],
                        params.quality_binning,
                    )?)
                } else {
                    let mut cl_lens: Vec<u32> = Vec::with_capacity(e - s);
                    let mut cl_qual: Vec<u8> = Vec::new();
                    for &oi in &plan.order[s..e] {
                        let oi = oi as usize;
                        cl_lens.push(all.lens[oi]);
                        let q = &all.qual[offs[oi]..offs[oi + 1]];
                        if plan.flip[oi] {
                            cl_qual.extend(q.iter().rev());
                        } else {
                            cl_qual.extend_from_slice(q);
                        }
                    }
                    Ok(fqxv_fqzcomp::encode(
                        &cl_lens,
                        &cl_qual,
                        params.quality_binning,
                    )?)
                }
            })
            .collect::<Result<_>>()
    })?;

    let nq_blocks: Vec<(Vec<u8>, Vec<u8>)> = name_blocks.into_iter().zip(qual_blocks).collect();

    // 7. Write: header, then n / flip / perm / seq blocks / name+qual blocks.
    let platform = resolve_platform_block(params.platform, &all);
    let mut w = BufWriter::new(writer);
    let mut flags = FLAG_PLUS_NORMALIZED | FLAG_REORDERED | FLAG_GLOBAL_REORDER;
    if keep_order {
        flags |= FLAG_KEEP_ORDER;
    }
    if template.is_some() {
        flags |= FLAG_REGEN_NAMES;
    }
    if use_reference {
        flags |= FLAG_GLOBAL_REFERENCE;
    }
    // A global reference frame is a whole-file capability a reader must have to
    // reconstruct the referenced blocks, and it's known here before any block is
    // written — so advertise it in the coarse feature word too (the flag bit stays
    // for in-layout dispatch; the feature word is what an older reader checks).
    let required_features = if use_reference {
        crate::feature::GLOBAL_REFERENCE
    } else {
        0
    };
    write_header_prefix(
        &mut w,
        params.seq_order,
        binning_tag(params.quality_binning),
        flags,
        g,
        platform,
        required_features,
    )?;
    w.write_all(&(n as u64).to_le_bytes())?;
    write_framed(&mut w, &flip_bits)?;
    write_framed(&mut w, &perm_c)?;
    // Name-template frame (empty unless regenerating names).
    let tmpl_bytes = template.as_ref().map(|t| t.to_bytes()).unwrap_or_default();
    write_framed(&mut w, &tmpl_bytes)?;
    // Shared global reference frame (only when the v4 layout was chosen); the
    // FLAG_GLOBAL_REFERENCE bit tells the decoder whether to expect it here.
    if use_reference {
        write_framed(&mut w, &ref_payload)?;
    }
    w.write_all(&(ranges.len() as u32).to_le_bytes())?;
    for payload in &seq_blocks {
        write_framed(&mut w, payload)?;
    }
    for (names, qual) in &nq_blocks {
        write_framed(&mut w, names)?;
        write_framed(&mut w, qual)?;
    }

    // Trailing whole-output content digest: fold the reads exactly as decode will
    // emit them (see [`OutputDigest`]). keep-order emits original order/content;
    // otherwise clustered order, original orientation, template-regenerated names.
    // Quality is folded *post-binning* — the values actually stored and
    // reconstructed — so a lossy-quality archive doesn't trip its own check.
    let mut od = OutputDigest::new();
    let binning = params.quality_binning;
    let read_slice = |a: usize| {
        (
            &all.seq[offs[a]..offs[a + 1]],
            &all.qual[offs[a]..offs[a + 1]],
        )
    };
    if keep_order {
        for i in 0..n {
            let (seq, qual) = read_slice(i);
            od.push(all.header(i), seq, &apply_binning(qual, binning));
        }
    } else {
        for j in 0..n {
            let oi = plan.order[j] as usize;
            let (seq, qual) = read_slice(oi);
            let regen;
            let name: &[u8] = if let Some(t) = &template {
                regen = t.regenerate(j);
                &regen
            } else {
                all.header(oi)
            };
            od.push(name, seq, &apply_binning(qual, binning));
        }
    }
    let output_digest = od.finish();
    write_framed(&mut w, &output_digest.to_le_bytes())?;
    w.flush()?;

    // Each framed slice is [4 len][4 crc][bytes]; n_blocks is a bare [4].
    let frame = |len: usize| 4 + CRC_LEN + len;
    let ref_frame = if use_reference {
        frame(ref_payload.len())
    } else {
        0
    };
    let out_bytes = (HEADER_LEN
        + 8
        + frame(flip_bits.len())
        + frame(perm_c.len())
        + ref_frame
        + 4
        + seq_blocks.iter().map(|p| frame(p.len())).sum::<usize>()
        + nq_blocks
            .iter()
            .map(|(nm, q)| frame(nm.len()) + frame(q.len()))
            .sum::<usize>()
        + frame(DIGEST_LEN)) as u64;
    Ok(Stats {
        reads: n as u64,
        blocks: ranges.len() as u64,
        out_bytes,
        group_size: g,
    })
}

/// The decoded whole-file reorder streams, before any un-permutation. `cl_reads`
/// is in clustered order; `names`/`lens`/`quals` are in clustered order without
/// `keep_order` and in original order with it (see [`encode_reordered`]).
pub(crate) struct ReorderStreams {
    n: usize,
    n_blocks: usize,
    flip: Vec<u8>,
    perm_c: Vec<u8>,
    cl_reads: Vec<Vec<u8>>,
    names: Vec<Vec<u8>>,
    lens: Vec<u32>,
    quals: Vec<u8>,
    /// When set (discard-order archives), `names` is empty and each output name
    /// is regenerated from this template at its output position.
    template: Option<fqxv_tokenizer::NameTemplate>,
    /// Whole-output content digest (see [`OutputDigest`]); the decode paths fold
    /// the reads they emit and compare against this.
    output_digest: u64,
}

/// Reference-frame coding method (first byte of the `FLAG_GLOBAL_REFERENCE`
/// frame). The reference is coded once per file with our own clean-room order-k
/// [`fqxv_seq`] model; the default splits it into fixed blocks coded in parallel.
const REF_METHOD_SEQ: u8 = 0; // whole-reference order-k `fqxv_seq` (single pass)
const REF_METHOD_SEQPAR: u8 = 1; // block-parallel order-k `fqxv_seq`
const REF_METHOD_PACK: u8 = 4; // SPRING-style 2-bit-pack + clean-room LZMA on packed
                               // Method tags 2 (LZ77) and 3 (BWT) were prototypes that lost to the order-k model
                               // on the raw bases and have been removed (see `refpack`'s module docs / issue #52).

/// Fixed number of contig blocks the reference is split into for parallel coding.
/// Fixed (never derived from the thread count) so the coded bytes are identical
/// regardless of `--threads`. Enough blocks to keep many cores busy; the
/// per-block context reset is a small ratio cost amortized over ~1M+ bases each
/// on real references.
const REF_SEQ_BLOCKS: usize = 64;

/// Encode the global reference frame, tagged with a leading method byte, picking
/// the smaller of two clean-room codings (never-worse):
///
/// * `REF_METHOD_PACK` — **SPRING's approach**: 2-bit-pack the ACGT consensus (a
///   hard 2 bits/base floor) and LZMA the packed bytes (plus the contig lengths),
///   which captures the long-range near-duplicate-contig repeats the order-k model
///   cannot see. This is the usual winner on real references. Measured −5.9% vs
///   the order-k model on 4M NovaSeq (20.71 vs 22.01 MB), close to SPRING's own
///   2-bit+BSC (~20.4 MB) and fast (~8s on the 22 MB packed stream, 4× smaller
///   than the raw bases).
/// * `REF_METHOD_SEQPAR` — the block-parallel order-k `fqxv_seq` fallback, which
///   wins only on references with little cross-contig repeat (e.g. tiny inputs).
///
/// Both are internally parallel and thread-count-independent, so the chosen bytes
/// are identical regardless of `--threads`. Fully clean-room — no external
/// compressor. (Earlier LZ77/BWT prototypes and a raw-base LZMA lost to the
/// order-k model — feeding byte-domain coders the raw 4-symbol bases underperforms;
/// the 2-bit packing is what makes the LZMA path win. They were removed.)
fn encode_reference_frame(
    reference: &fqxv_reorder::GlobalReference,
    order: usize,
) -> Result<Vec<u8>> {
    // The two codings are independent; overlap them so the never-worse gate costs
    // wall-clock ~max(pack, seqpar), not their sum.
    let (pack, seqpar) = rayon::join(
        || reference.encode_packed(),
        || reference.encode_blocked(order, REF_SEQ_BLOCKS),
    );
    let (method, payload) = match (pack?, seqpar?) {
        (pack, seqpar) if pack.len() <= seqpar.len() => (REF_METHOD_PACK, pack),
        (_, seqpar) => (REF_METHOD_SEQPAR, seqpar),
    };
    let mut out = Vec::with_capacity(1 + payload.len());
    out.push(method);
    out.extend_from_slice(&payload);
    Ok(out)
}

/// Decode a reference frame written by [`encode_reference_frame`].
fn decode_reference_frame(bytes: &[u8]) -> Result<fqxv_reorder::GlobalReference> {
    match bytes.split_first() {
        Some((&REF_METHOD_SEQ, rest)) => Ok(fqxv_reorder::GlobalReference::decode(rest)?),
        Some((&REF_METHOD_SEQPAR, rest)) => {
            Ok(fqxv_reorder::GlobalReference::decode_blocked(rest)?)
        }
        Some((&REF_METHOD_PACK, rest)) => Ok(fqxv_reorder::GlobalReference::decode_packed(rest)?),
        _ => Err(Error::Malformed("reorder reference: bad method tag")),
    }
}

/// Reserve a `Vec` of `cap` capacity, returning an error instead of aborting the
/// process when a corrupt `cap` — e.g. the out-of-range read count read straight
/// from a damaged stream — would demand an absurd allocation. `cap` is only a
/// sizing hint (the vectors grow via `extend` regardless), so this never changes
/// a valid decode; it just fails a bad archive gracefully, mirroring the
/// `try_reserve_exact` guard in [`read_framed`].
fn reserve_checked<T>(cap: usize) -> Result<Vec<T>> {
    let mut v = Vec::new();
    v.try_reserve_exact(cap)
        .map_err(|_| Error::Malformed("reorder read count too large to allocate"))?;
    Ok(v)
}

/// Read and entropy-decode the whole-file reorder layout. `r` is positioned just
/// past the header. Shared by [`decode_reordered_whole`] and
/// [`decode_reordered_split`]. `has_reference` (the `FLAG_GLOBAL_REFERENCE` bit)
/// says whether a shared global reference frame precedes the block count; when
/// set, version-4 sequence blocks are decoded as positions on it.
pub(crate) fn read_reordered_streams<R: Read>(
    mut r: R,
    pool: &rayon::ThreadPool,
    has_reference: bool,
) -> Result<ReorderStreams> {
    let mut n_buf = [0u8; 8];
    r.read_exact(&mut n_buf)?;
    let n = u64::from_le_bytes(n_buf) as usize;
    let flip = read_framed(&mut r, "reorder flip bitmap")?;
    let perm_c = read_framed(&mut r, "reorder permutation")?;
    let tmpl_bytes = read_framed(&mut r, "reorder name template")?;
    let template = if tmpl_bytes.is_empty() {
        None
    } else {
        Some(fqxv_tokenizer::NameTemplate::from_bytes(&tmpl_bytes)?)
    };
    let regen = template.is_some();
    // Shared global reference frame (present iff FLAG_GLOBAL_REFERENCE) — decoded
    // once, then every version-4 block indexes into it.
    let reference = if has_reference {
        let ref_bytes = read_framed(&mut r, "reorder global reference")?;
        Some(decode_reference_frame(&ref_bytes)?)
    } else {
        None
    };
    let mut nb = [0u8; 4];
    r.read_exact(&mut nb)?;
    let n_blocks = u32::from_le_bytes(nb) as usize;

    let mut seq_payloads: Vec<Vec<u8>> = Vec::with_capacity(n_blocks.min(1 << 20));
    for i in 0..n_blocks {
        seq_payloads.push(read_framed(&mut r, &format!("reorder sequence block {i}"))?);
    }
    let mut nq_payloads: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(n_blocks.min(1 << 20));
    for i in 0..n_blocks {
        let names = read_framed(&mut r, &format!("reorder name block {i}"))?;
        let qual = read_framed(&mut r, &format!("reorder quality block {i}"))?;
        nq_payloads.push((names, qual));
    }
    // Trailing whole-output content digest frame.
    let digest_bytes = read_framed(&mut r, "reorder output digest")?;
    let output_digest = u64::from_le_bytes(
        digest_bytes
            .as_slice()
            .try_into()
            .map_err(|_| Error::Malformed("reorder output digest length"))?,
    );

    // Decode both partitions in parallel. Version-4 blocks index the shared
    // reference; version-2/3 blocks ignore it (decode dispatches on the version
    // byte), so an archive may freely mix them.
    let reference_ref = reference.as_ref();
    let seq_dec: Vec<Vec<Vec<u8>>> = pool.install(|| {
        seq_payloads
            .par_iter()
            .map(|p| -> Result<Vec<Vec<u8>>> {
                Ok(fqxv_reorder::decode_clustered_any(p, reference_ref)?)
            })
            .collect::<Result<_>>()
    })?;
    // Per name+quality block: (decoded names, (per-read lengths, quality bytes)).
    type NqBlock = (Vec<Vec<u8>>, (Vec<u32>, Vec<u8>));
    let nq_dec: Vec<NqBlock> = pool.install(|| {
        nq_payloads
            .par_iter()
            .map(|(nm, q)| -> Result<_> {
                // Discard-order archives carry empty name blocks; names are
                // regenerated from the template, so skip name decoding.
                let names = if regen {
                    Vec::new()
                } else {
                    fqxv_tokenizer::decode(nm)?
                };
                Ok((names, fqxv_fqzcomp::decode(q)?))
            })
            .collect::<Result<_>>()
    })?;

    // Flatten the per-block vectors into whole-file streams. `n` is untrusted
    // (read straight from the stream), so size these hints fallibly — a corrupt
    // count must error, not abort the process on a multi-terabyte allocation.
    let mut cl_reads: Vec<Vec<u8>> = reserve_checked(n)?;
    for blk in seq_dec {
        cl_reads.extend(blk);
    }
    let mut names: Vec<Vec<u8>> = reserve_checked(if regen { 0 } else { n })?;
    let mut lens: Vec<u32> = reserve_checked(n)?;
    let mut quals: Vec<u8> = Vec::new();
    for (nm, (ls, qs)) in nq_dec {
        names.extend(nm);
        lens.extend(ls);
        quals.extend(qs);
    }
    let names_ok = if regen {
        names.is_empty()
    } else {
        names.len() == n
    };
    if cl_reads.len() != n || !names_ok || lens.len() != n {
        return Err(Error::Malformed("reordered stream length disagreement"));
    }
    Ok(ReorderStreams {
        n,
        n_blocks,
        flip,
        perm_c,
        cl_reads,
        names,
        lens,
        quals,
        template,
        output_digest,
    })
}

/// Un-permute a `keep_order` reorder archive: place each clustered sequence
/// (un-flipped) at its original position via the stored permutation, yielding the
/// sequences in original order. Consumes `s.cl_reads`.
pub(crate) fn unpermute_sequences(s: &mut ReorderStreams) -> Result<Vec<Vec<u8>>> {
    let n = s.n;
    let perm: Vec<u32> = {
        // The permutation is four byte-planes of `n` u32s, so `n * 4` is not an
        // estimate — it is the only length this can be, and it was already being
        // checked one line down. Checking it BEFORE the decode is what makes it a
        // bound: the length in the stream header is untrusted, and unbounded
        // `decode` allocates whatever it claims. The container CRCs each payload
        // first, but a CRC detects accidents, not a crafted archive.
        let pb = fqxv_rans::decode_bounded(&s.perm_c, n.saturating_mul(4))
            .map_err(|_| Error::Malformed("bad permutation"))?;
        if pb.len() != n * 4 {
            return Err(Error::Malformed("permutation length mismatch"));
        }
        (0..n)
            .map(|i| u32::from_le_bytes([pb[i], pb[n + i], pb[2 * n + i], pb[3 * n + i]]))
            .collect()
    };
    let mut seq_orig: Vec<Vec<u8>> = vec![Vec::new(); n];
    for (j, mut seq) in std::mem::take(&mut s.cl_reads).into_iter().enumerate() {
        if s.flip.get(j / 8).copied().unwrap_or(0) >> (j % 8) & 1 == 1 {
            seq = fqxv_reorder::revcomp(&seq);
        }
        let dest = perm[j] as usize;
        *seq_orig
            .get_mut(dest)
            .ok_or(Error::Malformed("permutation out of range"))? = seq;
    }
    Ok(seq_orig)
}

/// Decode a whole-file globally-clustered reorder archive to interleaved FASTQ on
/// a single writer (see [`compress_reordered_whole`]). `r` is positioned just past
/// the header. `keep_order` (from `FLAG_KEEP_ORDER`) selects the mode: un-permute
/// into original order, or emit in clustered order. `group_size` is recorded in
/// the returned [`Stats`]; grouped archives are always `keep_order`, so their
/// records emerge in original spot-interleaved order. `has_reference` (the
/// `FLAG_GLOBAL_REFERENCE` bit) says whether a shared reference frame is present.
pub(crate) fn decode_reordered_whole<R: Read, W: Write>(
    r: R,
    writer: W,
    threads: usize,
    keep_order: bool,
    group_size: u8,
    has_reference: bool,
) -> Result<Stats> {
    let pool = build_pool(threads)?;
    let mut s = read_reordered_streams(r, &pool, has_reference)?;
    let n = s.n;
    let n_blocks = s.n_blocks;
    let expected_digest = s.output_digest;
    let mut od = OutputDigest::new();
    let mut w = BufWriter::new(writer);
    if keep_order {
        // Un-permute, then emit in original order against the original-order
        // names/quality.
        let seq_orig = unpermute_sequences(&mut s)?;
        let mut qoff = 0usize;
        for i in 0..n {
            let l = s.lens[i] as usize;
            let qual = s
                .quals
                .get(qoff..qoff + l)
                .ok_or(Error::Malformed("quality underrun"))?;
            qoff += l;
            if seq_orig[i].len() != l {
                return Err(Error::Malformed("reordered sequence length mismatch"));
            }
            od.push(&s.names[i], &seq_orig[i], qual);
            let mut rec = Vec::with_capacity(l * 2 + s.names[i].len() + 8);
            write_record(&mut rec, &s.names[i], &seq_orig[i], qual);
            w.write_all(&rec)?;
        }
    } else {
        // Reads emerge in clustered order; names/quality were coded clustered too.
        // Un-flip the reverse-complemented reads (sequence and quality) to restore
        // each record's original content, then emit in clustered order.
        let template = s.template.take();
        let cl_reads = std::mem::take(&mut s.cl_reads);
        let mut qoff = 0usize;
        for (j, mut seq) in cl_reads.into_iter().enumerate() {
            let l = s.lens[j] as usize;
            let mut qual = s
                .quals
                .get(qoff..qoff + l)
                .ok_or(Error::Malformed("quality underrun"))?
                .to_vec();
            qoff += l;
            if seq.len() != l {
                return Err(Error::Malformed("reordered sequence length mismatch"));
            }
            if s.flip.get(j / 8).copied().unwrap_or(0) >> (j % 8) & 1 == 1 {
                seq = fqxv_reorder::revcomp(&seq);
                qual.reverse();
            }
            // Discard-order archives regenerate the name from the template at the
            // output position; otherwise use the clustered-order decoded name.
            let regen_name;
            let name: &[u8] = if let Some(t) = &template {
                regen_name = t.regenerate(j);
                &regen_name
            } else {
                &s.names[j]
            };
            od.push(name, &seq, &qual);
            let mut rec = Vec::with_capacity(l * 2 + name.len() + 8);
            write_record(&mut rec, name, &seq, &qual);
            w.write_all(&rec)?;
        }
    }
    w.flush()?;
    if od.finish() != expected_digest {
        return Err(Error::Corrupt {
            what: "reorder output digest".to_string(),
        });
    }
    Ok(Stats {
        reads: n as u64,
        blocks: n_blocks as u64,
        out_bytes: 0,
        group_size: group_size.max(1),
    })
}

/// Decode a grouped whole-file reorder archive, splitting the reads back into `g`
/// writers by their per-spot member. Only valid for `keep_order` archives (the
/// permutation reconstructs the mate interleaving); the caller guarantees this.
/// Record `i` in restored original order belongs to member `i % g`.
/// `has_reference` (the `FLAG_GLOBAL_REFERENCE` bit) says whether a shared
/// reference frame is present.
pub(crate) fn decode_reordered_split<R: Read, W: Write>(
    r: R,
    writers: &mut [W],
    threads: usize,
    g: usize,
    has_reference: bool,
) -> Result<Stats> {
    let pool = build_pool(threads)?;
    let mut s = read_reordered_streams(r, &pool, has_reference)?;
    let n = s.n;
    let n_blocks = s.n_blocks;
    let expected_digest = s.output_digest;
    let seq_orig = unpermute_sequences(&mut s)?;
    let mut bufs: Vec<BufWriter<&mut W>> = writers.iter_mut().map(BufWriter::new).collect();
    let mut stats = Stats {
        reads: n as u64,
        blocks: n_blocks as u64,
        group_size: g as u8,
        ..Stats::default()
    };
    let mut od = OutputDigest::new();
    let mut qoff = 0usize;
    for i in 0..n {
        let l = s.lens[i] as usize;
        let qual = s
            .quals
            .get(qoff..qoff + l)
            .ok_or(Error::Malformed("quality underrun"))?;
        qoff += l;
        if seq_orig[i].len() != l {
            return Err(Error::Malformed("reordered sequence length mismatch"));
        }
        od.push(&s.names[i], &seq_orig[i], qual);
        let mut rec = Vec::with_capacity(l * 2 + s.names[i].len() + 8);
        write_record(&mut rec, &s.names[i], &seq_orig[i], qual);
        bufs[i % g].write_all(&rec)?;
        stats.out_bytes += rec.len() as u64;
    }
    for b in &mut bufs {
        b.flush()?;
    }
    // Split emits reads in original (i) order across the g writers, matching the
    // keep-order digest folded above.
    if od.finish() != expected_digest {
        return Err(Error::Corrupt {
            what: "reorder output digest".to_string(),
        });
    }
    Ok(stats)
}

/// Rolling xxh3-64 over reads in output order — the whole-file reorder layout's
/// analog of the plain layout's per-block, per-stream [`stream_digests`] (that
/// layout splits reads across seq/name/quality partitions, so there is no single
/// block to digest). Encode
/// folds the reads it *will* emit (original order for keep-order; clustered order,
/// original orientation, with template-regenerated names otherwise); decode folds
/// the reads it *actually* emits and compares. A mismatch means the reorder codec
/// stack (clustering, contig assembly, permutation, flips) round-tripped into
/// wrong output. Per-read name/seq lengths are folded to pin boundaries; the read
/// count is folded last so a short/long read set can't collide. `qual.len()`
/// equals `seq.len()` per read.
pub(crate) struct OutputDigest {
    h: Xxh3,
    n: u64,
}

impl OutputDigest {
    fn new() -> Self {
        OutputDigest {
            h: Xxh3::new(),
            n: 0,
        }
    }
    fn push(&mut self, name: &[u8], seq: &[u8], qual: &[u8]) {
        self.h.update(&(name.len() as u32).to_le_bytes());
        self.h.update(name);
        self.h.update(&(seq.len() as u32).to_le_bytes());
        self.h.update(seq);
        self.h.update(qual);
        self.n += 1;
    }
    fn finish(mut self) -> u64 {
        self.h.update(&self.n.to_le_bytes());
        self.h.digest()
    }
}

/// Fold quality through the active binning table so a content digest covers the
/// *stored* (post-binning) bytes that decode reconstructs — not the original
/// input. Borrows unchanged when lossless. Mirrors the per-block digest's rule.
pub(crate) fn apply_binning(qual: &[u8], binning: QualityBinning) -> Cow<'_, [u8]> {
    match binning {
        QualityBinning::Lossless => Cow::Borrowed(qual),
        b => Cow::Owned(qual.iter().map(|&q| b.apply(q)).collect()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use fqxv_reorder::GlobalReference;

    fn sample_reference() -> GlobalReference {
        // Three ACGT contigs; from_lens_seq validates the lengths sum to seq.len().
        let seq = b"ACGTACGTACGTTTTTGGGGCCCCAAAATTTTACACACACGTGTGTGT".to_vec();
        GlobalReference::from_lens_seq(&[16, 16, 16], seq).unwrap()
    }

    #[test]
    fn reference_frame_roundtrips_both_methods() {
        let r = sample_reference();

        // Default frame round-trips; the never-worse gate picks PACK or SEQPAR.
        let frame = encode_reference_frame(&r, 8).unwrap();
        assert!(matches!(frame[0], REF_METHOD_PACK | REF_METHOD_SEQPAR));
        let back = decode_reference_frame(&frame).unwrap();
        assert_eq!(back.raw_bases(), r.raw_bases());
        assert_eq!(back.contig_lens(), r.contig_lens());

        // The PACK method decodes when forced directly.
        let mut pack_frame = vec![REF_METHOD_PACK];
        pack_frame.extend_from_slice(&r.encode_packed().unwrap());
        let pack_back = decode_reference_frame(&pack_frame).unwrap();
        assert_eq!(pack_back.raw_bases(), r.raw_bases());
        assert_eq!(pack_back.contig_lens(), r.contig_lens());

        // The whole-reference (single-pass) method still decodes.
        let mut seq_frame = vec![REF_METHOD_SEQ];
        seq_frame.extend_from_slice(&r.encode(8, 0, 0).unwrap());
        let seq_back = decode_reference_frame(&seq_frame).unwrap();
        assert_eq!(seq_back.raw_bases(), r.raw_bases());
        assert_eq!(seq_back.contig_lens(), r.contig_lens());

        // Block-parallel output is identical regardless of the block count
        // (fixed contig-index boundaries) — spot-check a couple of block counts
        // decode to the same reference.
        for nb in [1usize, 3, 7] {
            let payload = r.encode_blocked(8, nb).unwrap();
            let rt = GlobalReference::decode_blocked(&payload).unwrap();
            assert_eq!(rt.raw_bases(), r.raw_bases());
            assert_eq!(rt.contig_lens(), r.contig_lens());
        }
    }

    #[test]
    fn reference_frame_rejects_bad_method() {
        assert!(decode_reference_frame(&[]).is_err());
        assert!(decode_reference_frame(&[0xff, 1, 2, 3]).is_err());
    }
}
