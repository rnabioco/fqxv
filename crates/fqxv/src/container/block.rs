//! Per-block coding: `RawBlock`, block compress/decode, and stream framing.

use super::*;
use tracing::trace;

/// One block of parsed FASTQ records. Header text is packed into a single arena
/// (`header_buf` + cumulative `header_ends`) rather than a `Vec` per record — the
/// parse loop is single-threaded and feeds the parallel compressors, so avoiding
/// a per-read allocation keeps that feed from starving the pool.
#[derive(Default)]
pub(crate) struct RawBlock {
    pub(crate) header_buf: Vec<u8>,
    pub(crate) header_ends: Vec<u32>,
    pub(crate) lens: Vec<u32>,
    pub(crate) seq: Vec<u8>,
    pub(crate) qual: Vec<u8>,
}

impl RawBlock {
    pub(crate) fn push(&mut self, name: &[u8], description: &[u8], seq: &[u8], qual: &[u8]) {
        self.header_buf.extend_from_slice(name);
        if !description.is_empty() {
            self.header_buf.push(b' ');
            self.header_buf.extend_from_slice(description);
        }
        self.header_ends.push(self.header_buf.len() as u32);
        self.lens.push(seq.len() as u32);
        self.seq.extend_from_slice(seq);
        self.qual.extend_from_slice(qual);
    }

    /// Append a record whose header text is already in its final (normalized)
    /// form — see [`normalize_header`]. Used by the parallel parser, which builds
    /// the name/description join itself and hands over the finished header bytes.
    pub(crate) fn push_raw(&mut self, header: &[u8], seq: &[u8], qual: &[u8]) {
        self.header_buf.extend_from_slice(header);
        self.header_ends.push(self.header_buf.len() as u32);
        self.lens.push(seq.len() as u32);
        self.seq.extend_from_slice(seq);
        self.qual.extend_from_slice(qual);
    }

    /// Number of records in the block.
    pub(crate) fn n_reads(&self) -> usize {
        self.header_ends.len()
    }

    /// The `i`th record's header bytes.
    pub(crate) fn header(&self, i: usize) -> &[u8] {
        let start = if i == 0 {
            0
        } else {
            self.header_ends[i - 1] as usize
        };
        &self.header_buf[start..self.header_ends[i] as usize]
    }

    /// Borrowed slices for every header, in record order.
    pub(crate) fn header_refs(&self) -> Vec<&[u8]> {
        let mut refs = Vec::with_capacity(self.header_ends.len());
        let mut start = 0usize;
        for &end in &self.header_ends {
            refs.push(&self.header_buf[start..end as usize]);
            start = end as usize;
        }
        refs
    }
}

/// Assemble one output block from the globally-ordered record range `[gs, ge)`,
/// copying each record's header (from its chunk's arena) and sequence/quality
/// (from `buf`) into a fresh [`RawBlock`]. `gstart[c]` is the global index of
/// chunk `c`'s first record.
pub(crate) fn build_block(
    buf: &[u8],
    chunks: &[ChunkParse],
    gstart: &[usize],
    gs: usize,
    ge: usize,
) -> RawBlock {
    let mut blk = RawBlock::default();
    let mut gi = gs;
    // Chunk holding the first record: the last c with gstart[c] <= gi.
    let mut c = gstart.partition_point(|&x| x <= gi) - 1;
    while gi < ge {
        let chunk = &chunks[c];
        let base = gstart[c];
        let local_start = gi - base;
        let take = (ge.min(gstart[c + 1]) - gi) + local_start;
        for local in local_start..take {
            let rec = chunk.recs[local];
            let hdr_start = if local == 0 {
                0
            } else {
                chunk.recs[local - 1].hdr_end as usize
            };
            let header = &chunk.hdr[hdr_start..rec.hdr_end as usize];
            let seq = &buf[rec.seq_off..rec.seq_off + rec.seq_len as usize];
            let qual = &buf[rec.qual_off..rec.qual_off + rec.qual_len as usize];
            blk.push_raw(header, seq, qual);
        }
        gi = ge.min(gstart[c + 1]);
        c += 1;
    }
    blk
}

/// Write a batch's compressed payloads in order, updating `stats` and recording
/// each row group in `index` for the footer.
pub(crate) fn write_blocks<W: Write>(
    w: &mut W,
    blocks: &[RawBlock],
    compressed: Vec<Result<Vec<u8>>>,
    stats: &mut Stats,
    index: &mut FooterIndex,
) -> Result<()> {
    for (b, payload) in blocks.iter().zip(compressed) {
        let payload = payload?;
        index.entries.push((index.offset, b.n_reads() as u32));
        // Record each coded stream's absolute (offset, len, crc) for the footer's
        // column-projection index, before `index.offset` advances past this block.
        index
            .streams
            .push(payload_stream_locs(&payload, index.offset)?);
        // Frame: [4 BLOCK_MAGIC][8 payload_len][4 crc32c(payload)][payload].
        w.write_all(&BLOCK_MAGIC)?;
        w.write_all(&(payload.len() as u64).to_le_bytes())?;
        w.write_all(&crc32c(&payload).to_le_bytes())?;
        w.write_all(&payload)?;
        let framed = (FRAME_HEAD_LEN + payload.len()) as u64;
        index.offset += framed;
        trace!(
            reads = b.n_reads(),
            payload = payload.len(),
            "block written"
        );
        stats.reads += b.n_reads() as u64;
        stats.blocks += 1;
        stats.out_bytes += framed;
    }
    Ok(())
}

/// One xxh3-64 digest per decoded stream (names, sequence, quality) of a block's
/// *decoded canonical form* — the exact bytes `decompress` reconstructs.
pub(crate) struct StreamDigests {
    pub(crate) names: u64,
    pub(crate) seq: u64,
    pub(crate) qual: u64,
}

/// Digest each of a block's three decoded streams independently, so a mismatch
/// localizes *which* stream a codec round-tripped into wrong-but-in-bounds output
/// (corruption the per-payload CRC cannot catch, because the stored bytes were
/// never altered). Computed identically on the encode side (from the post-
/// `QualityBinning` quality — the values actually stored) and the decode side (from
/// the reconstructed streams).
///
/// The digests are over the *stored* (post-binning) form, not the original input,
/// so a lossy archive verifies against what it emits, not against data it never
/// promised to reproduce. Each digest folds in `n_reads` and its stream's per-read
/// lengths (name lengths for names; `lens` for both sequence and quality, which
/// share it) so no byte can silently cross a read — or stream — boundary.
pub(crate) fn stream_digests<'a>(
    n_reads: usize,
    names: impl Iterator<Item = &'a [u8]>,
    lens: &[u32],
    seq: &[u8],
    qual: &[u8],
) -> StreamDigests {
    let mut hn = Xxh3::new();
    hn.update(&(n_reads as u64).to_le_bytes());
    for name in names {
        hn.update(&(name.len() as u32).to_le_bytes());
        hn.update(name);
    }

    // Sequence and quality each pin their per-read boundaries with the shared
    // `lens`, so a byte sliding between reads changes the digest even at constant
    // total length.
    let digest_lens_then = |bytes: &[u8]| {
        let mut h = Xxh3::new();
        h.update(&(n_reads as u64).to_le_bytes());
        for &l in lens {
            h.update(&l.to_le_bytes());
        }
        h.update(bytes);
        h.digest()
    };

    StreamDigests {
        names: hn.digest(),
        seq: digest_lens_then(seq),
        qual: digest_lens_then(qual),
    }
}

/// Sequence-stream codec, recorded as a leading method byte (v4). Short reads use
/// the order-k context model; long reads the cross-read overlap-assembly codec.
pub(crate) const SEQ_METHOD_ORDERK: u8 = 0;
/// See [`SEQ_METHOD_ORDERK`].
pub(crate) const SEQ_METHOD_OVERLAP: u8 = 1;
/// Long-read overlap codec coded against a **shared, whole-file reference**
/// (`fqxv_lroverlap::encode_against`): the block stores only its edit streams; the
/// consensus reference lives once in the archive's `FLAG_GLOBAL_REFERENCE` frame
/// (issue #168). Decoding requires that frame — [`decode_sequence_stream`] fails
/// closed if it is absent.
pub(crate) const SEQ_METHOD_OVERLAP_REF: u8 = 2;
/// Raw-sequence LZMA (`fqxv_reorder::lzma_seq_encode`): a clean-room LZMA over the
/// block's ASCII bases with an ~89 MB match window. At ordinary long-read coverage
/// (a real genome, not 300x of one organism) overlapping reads share long *exact*
/// substrings that neither the within-read order-k model nor the consensus-edit
/// overlap codec captures, but a large-window LZ finds directly — 0.6 vs 1.4
/// b/base on Revio WGS (#197). Self-contained; kept only when it beats the other
/// candidates via [`keep_smaller`].
pub(crate) const SEQ_METHOD_LZMA: u8 = 3;

/// Which index a sketch is being built for.
///
/// Seeding-scheme quality is **coverage-dependent**, and the long-read encoder
/// builds two indexes over very different amounts of data, so they do not want
/// the same scheme. See [`sketch_for`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum SeedContext {
    /// The shared whole-file reference: every read in the file contributes, so
    /// this index sees the input's full coverage.
    WholeFile,
    /// One block's own overlap index, which sees only that block's share of the
    /// coverage — on a file split into `n` blocks, roughly `1/n` of it.
    PerBlock,
}

/// Minimizer sketch for the long-read overlap codec, chosen by the detected
/// platform **and** by how much coverage the index will see.
///
/// PacBio's <1% error rate leaves nearly every k-mer intact, so its sparse
/// `map-hifi` sketch suffices at either coverage; everything else falls back to
/// the dense `map-ont` sketch, the conservative default that also works on HiFi
/// (it only costs index size) whereas HiFi's sparse sketch misses ONT overlaps.
///
/// **Why ONT splits on context.** Closed syncmers conserve anchors better than
/// window minimizers at ~10% error, but only pay off once coverage is deep
/// enough for the extra surviving anchors to find partners; below that their
/// weaker per-window specificity costs more than the conservation gains. The two
/// indexes land on opposite sides of that crossover, measured on `ecoli_ont`
/// (block 1, 268 Mbase):
///
/// | index | syncmer | minimizer |
/// |---|---:|---:|
/// | whole-file reference | **1.280 b/base** | 1.416 |
/// | per-block overlap | 1.517 | **1.243 b/base** |
///
/// A single scheme across both therefore gives up ~18% on one index or ~10% on
/// the other. Using one scheme everywhere (syncmers, issue #177) is what put the
/// ONT archive 2.79 MB behind; splitting recovers it. Both schemes share `w` and
/// `k`, so anchor density (`2/(w + 1)`) and specificity are unchanged — only
/// *which* positions are selected differs.
///
/// The sketch affects ratio and speed only — the block self-describes `(w, k)`,
/// so decode is unaffected, and the keep-the-smaller rule in
/// [`encode_sequence_stream`] means a mis-detected platform can never regress a
/// block below order-k.
pub(crate) fn sketch_for(platform: Platform, ctx: SeedContext) -> fqxv_lroverlap::Sketch {
    match (platform, ctx) {
        (Platform::PacBio, _) => fqxv_lroverlap::Sketch::hifi(),
        (_, SeedContext::WholeFile) => fqxv_lroverlap::Sketch::ont(),
        (_, SeedContext::PerBlock) => fqxv_lroverlap::Sketch {
            scheme: fqxv_lroverlap::SeedScheme::Minimizer,
            ..fqxv_lroverlap::Sketch::ont()
        },
    }
}

fn tag_orderk(coded: Vec<u8>) -> Vec<u8> {
    let mut out = Vec::with_capacity(coded.len() + 1);
    out.push(SEQ_METHOD_ORDERK);
    out.extend_from_slice(&coded);
    out
}

/// Raw-sequence LZMA candidate (see [`SEQ_METHOD_LZMA`]): a method-tagged stream
/// coding the block's ASCII bases with a large-window LZ. It wins where reads
/// carry cross-read redundancy at ordinary coverage; the caller keeps it only when
/// it is smaller than the other candidates.
fn lzma_stream(lens: &[u32], seq: &[u8]) -> Result<Vec<u8>> {
    let coded = fqxv_reorder::lzma_seq_encode(lens, seq)?;
    let mut out = Vec::with_capacity(coded.len() + 1);
    out.push(SEQ_METHOD_LZMA);
    out.extend_from_slice(&coded);
    Ok(out)
}

/// Whether to code the raw-LZMA sequence candidate (#197). `FQXV_SEQ_NO_LZMA`
/// disables it — for A/B measurement of the LZMA lever, and for tests that must
/// exercise a layout LZMA would otherwise win. Off by default (LZMA on), zero-cost
/// when unset.
fn lzma_candidate_enabled() -> bool {
    std::env::var_os("FQXV_SEQ_NO_LZMA").is_none()
}

/// The optional LZMA candidate: `None` (skipping the encode entirely) when
/// disabled, else the coded stream.
fn lzma_candidate(lens: &[u32], seq: &[u8]) -> Result<Option<Vec<u8>>> {
    if lzma_candidate_enabled() {
        Ok(Some(lzma_stream(lens, seq)?))
    } else {
        Ok(None)
    }
}

pub(crate) fn order_k_stream(lens: &[u32], seq: &[u8], params: &Params) -> Result<Vec<u8>> {
    let plain = fqxv_seq::encode(lens, seq, params.seq_order as usize)?;

    // The hashed high-order tier (levels 8+, `seq_hash_order > seq_order`) adds an
    // escape tier above order-k. It captures deep context on repetitive data but is
    // not free: on a low-redundancy library its escapes cost more than they save,
    // which made `-l9`/`--max` produce a *larger* archive than the default (#196).
    // It is auto-detected from the stream on decode, so the tier-on and tier-off
    // encodings decode identically — code both and keep the smaller, the same
    // never-worse rule the overlap-vs-order-k choice already uses. `min(plain,
    // hashed)` also floors this at the plain order-k cost, so enabling the tier can
    // never regress a block.
    let hashed_on = params.seq_hash_bits > 0
        && usize::from(params.seq_hash_order) > usize::from(params.seq_order);
    if !hashed_on {
        return Ok(tag_orderk(plain));
    }
    let hashed = fqxv_seq::encode_hashed(
        lens,
        seq,
        params.seq_order as usize,
        params.seq_hash_order as usize,
        u32::from(params.seq_hash_bits),
    )?;
    Ok(tag_orderk(keep_smaller(hashed, plain)))
}

/// Encode the sequence stream, choosing the codec by a leading method byte.
///
/// Short reads always use the order-k context model. Long reads (mean length over
/// [`REORDER_MAX_MEAN_LEN`]) carry cross-read redundancy the within-read model
/// cannot see, so the overlap-assembly codec is tried — but its advantage
/// depends on within-block coverage, which the per-block budget bounds, so both
/// are coded and the **smaller** kept. The overlap codec then never regresses a
/// block against order-k (low coverage, a large genome, or few reads all fall
/// back automatically), the same "keep the smaller" rule the reference coder and
/// the per-stream entropy coder already use. The choice is per block, so a file
/// of mixed lengths codes each block with whichever fits.
///
/// `platform` selects the overlap codec's minimizer sketch (see [`sketch_for`]);
/// short-read blocks never reach the overlap codec, so it is unused there.
pub(crate) fn encode_sequence_stream(
    lens: &[u32],
    seq: &[u8],
    params: &Params,
    platform: Platform,
) -> Result<Vec<u8>> {
    if !super::reorder::is_long_read(lens) {
        return order_k_stream(lens, seq, params);
    }
    // Three candidates, coded in parallel and kept by the smaller (#197). The
    // overlap codec assembles a reference from THIS block's reads only (a per-block
    // index, see `SeedContext`); order-k is the within-read fallback; and raw LZMA
    // catches the cross-read *exact* matches a real genome at ordinary coverage
    // carries, which the other two miss.
    let opts = fqxv_lroverlap::EncodeOpts {
        sketch: sketch_for(platform, SeedContext::PerBlock),
        ..Default::default()
    };
    let (overlap, (order_k, lzma)) = rayon::join(
        || {
            fqxv_lroverlap::encode(lens, seq, &opts).map(|coded| {
                let mut out = Vec::with_capacity(coded.len() + 1);
                out.push(SEQ_METHOD_OVERLAP);
                out.extend_from_slice(&coded);
                out
            })
        },
        || {
            rayon::join(
                || order_k_stream(lens, seq, params),
                || lzma_candidate(lens, seq),
            )
        },
    );
    let (overlap, order_k, lzma) = (overlap?, order_k?, lzma?);
    // `FQXV_DIAG_SEQ` reports the per-block contest (bytes and bits/base) so
    // seeding/consensus changes can be judged against the real keep-the-smaller
    // outcome, not a proxy. Off by default, zero-cost.
    if std::env::var_os("FQXV_DIAG_SEQ").is_some() {
        let bases: u64 = lens.iter().map(|&l| u64::from(l)).sum::<u64>().max(1);
        let bpb = |n: usize| (n as f64 * 8.0) / bases as f64;
        let lzma_len = lzma.as_ref().map_or(usize::MAX, Vec::len);
        eprintln!(
            "[diag seq] {} reads, {} bases | overlap {} B ({:.3} b/base) | order-k {} B ({:.3} b/base) | lzma {} B ({:.3} b/base)",
            lens.len(),
            bases,
            overlap.len(),
            bpb(overlap.len()),
            order_k.len(),
            bpb(order_k.len()),
            lzma_len,
            bpb(lzma_len),
        );
    }
    let plain = keep_smaller(overlap, order_k);
    Ok(match lzma {
        Some(l) => keep_smaller(l, plain),
        None => plain,
    })
}

/// Code **only** the shared-reference candidate: the smaller of the
/// reference-coded [`SEQ_METHOD_OVERLAP_REF`] stream and order-k.
///
/// The caller has already committed to the reference layout (see the block-0
/// probe in `compress_longread_shared_ref`), so the per-block overlap candidate
/// would be assembled and then thrown away. Skipping it is the difference
/// between one and two full long-read assemblies per block, which on HiFi — where
/// the shared reference wins by ~30x and the plain candidate never has a chance —
/// is the bulk of compress time.
///
/// Order-k is still coded, so a block whose reads do not place on the reference
/// cannot regress below the context model.
pub(crate) fn encode_sequence_stream_shared_only(
    lens: &[u32],
    seq: &[u8],
    params: &Params,
    platform: Platform,
    reference: &fqxv_lroverlap::Reference,
) -> Result<Vec<u8>> {
    let opts_shared = fqxv_lroverlap::EncodeOpts {
        sketch: sketch_for(platform, SeedContext::WholeFile),
        ..Default::default()
    };
    let (against, order_k) = rayon::join(
        || {
            fqxv_lroverlap::encode_against(reference, lens, seq, &opts_shared).map(|coded| {
                let mut out = Vec::with_capacity(coded.len() + 1);
                out.push(SEQ_METHOD_OVERLAP_REF);
                out.extend_from_slice(&coded);
                out
            })
        },
        || order_k_stream(lens, seq, params),
    );
    let (against, order_k) = (against?, order_k?);
    Ok(keep_smaller(against, order_k))
}

/// Code a long-read block's sequence **both ways** and return the two candidate
/// streams the whole-file layout gate chooses between (issue #184):
///
/// - `.0` — the **shared-reference** candidate: the smaller of the reference-coded
///   [`SEQ_METHOD_OVERLAP_REF`] stream and order-k. Valid only in an archive that
///   carries the reference frame.
/// - `.1` — the **plain** candidate: exactly what [`encode_sequence_stream`] would
///   produce for this block, i.e. the smallest of the per-block
///   [`SEQ_METHOD_OVERLAP`] stream, order-k, and raw [`SEQ_METHOD_LZMA`].
///   Self-contained.
///
/// Both are needed because the gate must compare the shared-reference layout
/// against *the layout it would otherwise use*. Gating on order-k alone measures
/// against a bar weaker than the plain fallback (which floors each block at
/// `min(overlap, order-k)`), so a reference that loses to the per-block overlap
/// codec still clears it — the ONT regression in issue #184, where the frame cost
/// 4.37 MB to save 1.58 MB and was adopted anyway.
///
/// Returning the plain candidate rather than just its length lets the caller write
/// the fallback layout without re-coding: the per-block overlap encode is the
/// expensive half of this function, and `block_ranges` is a pure function of the
/// input, so these streams are byte-identical to a fresh
/// [`compress_buffered_plain`] pass.
pub(crate) fn encode_sequence_stream_shared(
    lens: &[u32],
    seq: &[u8],
    params: &Params,
    platform: Platform,
    reference: &fqxv_lroverlap::Reference,
) -> Result<(Vec<u8>, Vec<u8>)> {
    // The two candidates index different amounts of data and take the seeding
    // scheme that suits each: `encode_against` places reads on the whole-file
    // reference (built with the same whole-file sketch by the caller), while
    // `encode` assembles a reference from this block alone.
    let opts_shared = fqxv_lroverlap::EncodeOpts {
        sketch: sketch_for(platform, SeedContext::WholeFile),
        ..Default::default()
    };
    let opts_block = fqxv_lroverlap::EncodeOpts {
        sketch: sketch_for(platform, SeedContext::PerBlock),
        ..Default::default()
    };
    let (against, (overlap, (order_k, lzma))) = rayon::join(
        || {
            fqxv_lroverlap::encode_against(reference, lens, seq, &opts_shared).map(|coded| {
                let mut out = Vec::with_capacity(coded.len() + 1);
                out.push(SEQ_METHOD_OVERLAP_REF);
                out.extend_from_slice(&coded);
                out
            })
        },
        || {
            rayon::join(
                || {
                    fqxv_lroverlap::encode(lens, seq, &opts_block).map(|coded| {
                        let mut out = Vec::with_capacity(coded.len() + 1);
                        out.push(SEQ_METHOD_OVERLAP);
                        out.extend_from_slice(&coded);
                        out
                    })
                },
                || {
                    rayon::join(
                        || order_k_stream(lens, seq, params),
                        || lzma_candidate(lens, seq),
                    )
                },
            )
        },
    );
    let (against, overlap, order_k, lzma) = (against?, overlap?, order_k?, lzma?);
    if std::env::var_os("FQXV_DIAG_SEQ").is_some() {
        let bases: u64 = lens.iter().map(|&l| u64::from(l)).sum::<u64>().max(1);
        let bpb = |n: usize| (n as f64 * 8.0) / bases as f64;
        let lzma_len = lzma.as_ref().map_or(usize::MAX, Vec::len);
        eprintln!(
            "[diag seq/shared] {} reads, {} bases | shared-ref {} B ({:.3} b/base) | per-block overlap {} B ({:.3} b/base) | order-k {} B ({:.3} b/base) | lzma {} B ({:.3} b/base)",
            lens.len(),
            bases,
            against.len(),
            bpb(against.len()),
            overlap.len(),
            bpb(overlap.len()),
            order_k.len(),
            bpb(order_k.len()),
            lzma_len,
            bpb(lzma_len),
        );
    }
    // Shared candidate: the reference-coded stream vs order-k. Plain candidate: the
    // layout the whole-file gate falls back to — the smallest of the per-block
    // overlap codec, order-k, and raw LZMA (#197). Order-k feeds both sides, so it
    // is cloned once (a cheap memcpy beside the encodes above).
    let shared = keep_smaller(against, order_k.clone());
    let mut plain = keep_smaller(overlap, order_k);
    if let Some(l) = lzma {
        plain = keep_smaller(l, plain);
    }
    Ok((shared, plain))
}

/// Decode a method-tagged sequence stream back to `(per-read lengths, bases)`.
///
/// `shared_ref` is the archive's whole-file reference frame (present only when the
/// header's `FLAG_GLOBAL_REFERENCE` bit is set). A [`SEQ_METHOD_OVERLAP_REF`] block
/// is coded against it and cannot be decoded without it, so this **fails closed**
/// (`Error::Malformed`) when the method needs the reference but none was supplied.
pub(crate) fn decode_sequence_stream(
    coded: &[u8],
    shared_ref: Option<&fqxv_lroverlap::Reference>,
) -> Result<(Vec<u32>, Vec<u8>)> {
    let (&method, rest) = coded
        .split_first()
        .ok_or(Error::Malformed("empty sequence stream"))?;
    match method {
        SEQ_METHOD_ORDERK => Ok(fqxv_seq::decode(rest)?),
        SEQ_METHOD_LZMA => Ok(fqxv_reorder::lzma_seq_decode(rest)?),
        SEQ_METHOD_OVERLAP => Ok(fqxv_lroverlap::decode(rest)?),
        SEQ_METHOD_OVERLAP_REF => {
            let reference = shared_ref.ok_or(Error::Malformed(
                "shared-reference sequence block but no reference frame in the archive",
            ))?;
            Ok(fqxv_lroverlap::decode_against(reference, rest)?)
        }
        method => Err(Error::UnsupportedMethod {
            stream: "sequence",
            method,
        }),
    }
}

/// Code one non-reorder block: names (tokenizer), sequence (method-tagged), and
/// quality (fqzcomp), each length-prefixed, behind three leading per-stream
/// [`stream_digests`]. Reorder uses the whole-file path instead.
pub(crate) fn compress_block(b: &RawBlock, params: &Params, platform: Platform) -> Result<Vec<u8>> {
    let header_refs = b.header_refs();
    // The three streams are independent; code them concurrently so a block's
    // wall time is its slowest stream, not their sum. Nested inside the
    // per-block `par_iter`, these joins simply fill cores left idle when there
    // are fewer blocks than threads (and are near-free when every worker is
    // already busy).
    let (names_c, (seq_c, qual_c)) = rayon::join(
        || fqxv_tokenizer::encode(&header_refs),
        || {
            rayon::join(
                || encode_sequence_stream(&b.lens, &b.seq, params, platform),
                // Hand the bases to the quality coder: on long reads it conditions
                // quality on sequence (base + next + homopolymer run); on short
                // reads it ignores them and codes the position context as before.
                || fqxv_fqzcomp::encode_seq(&b.lens, &b.qual, &b.seq, params.quality_binning),
            )
        },
    );
    let (names_c, seq_c, qual_c) = (names_c?, seq_c?, qual_c?);
    assemble_block_payload(b, &names_c, &seq_c, &qual_c, params)
}

/// Code one non-reorder block reusing an already-coded sequence stream: names
/// (tokenizer) and quality (fqzcomp) are coded here, and `precoded_seq` (a
/// method-tagged stream from [`encode_sequence_stream_shared`]) is bundled in
/// as-is. The shared-reference compress path (issue #168) codes the expensive
/// reference-aligned sequence in a first pass and reuses it here, so the sequence
/// is never coded twice.
pub(crate) fn compress_block_with_seq(
    b: &RawBlock,
    params: &Params,
    precoded_seq: &[u8],
) -> Result<Vec<u8>> {
    let header_refs = b.header_refs();
    let (names_c, qual_c) = rayon::join(
        || fqxv_tokenizer::encode(&header_refs),
        || fqxv_fqzcomp::encode_seq(&b.lens, &b.qual, &b.seq, params.quality_binning),
    );
    let (names_c, qual_c) = (names_c?, qual_c?);
    assemble_block_payload(b, &names_c, precoded_seq, &qual_c, params)
}

/// Assemble a block payload from its three already-coded streams: prepend the
/// per-stream decoded-content digests and `n_reads`, then each length-prefixed
/// stream. Shared by the plain ([`compress_block`]) and shared-reference
/// ([`compress_block_shared`]) block builders, which differ only in how they code
/// the sequence stream.
fn assemble_block_payload(
    b: &RawBlock,
    names_c: &[u8],
    seq_c: &[u8],
    qual_c: &[u8],
    params: &Params,
) -> Result<Vec<u8>> {
    // End-to-end round-trip check: digest the block's decoded content (post-binning
    // quality, so lossy archives verify against what they emit) and store it at the
    // head of the payload. Lossless is the common case and borrows without a copy.
    let binned: Cow<[u8]> = match params.quality_binning {
        QualityBinning::Lossless => Cow::Borrowed(&b.qual),
        binning => Cow::Owned(b.qual.iter().map(|&q| binning.apply(q)).collect()),
    };
    let digests = stream_digests(
        b.n_reads(),
        b.header_refs().into_iter(),
        &b.lens,
        &b.seq,
        &binned,
    );

    let mut out =
        Vec::with_capacity(STREAM_DIGESTS_LEN + 16 + names_c.len() + seq_c.len() + qual_c.len());
    out.extend_from_slice(&digests.names.to_le_bytes());
    out.extend_from_slice(&digests.seq.to_le_bytes());
    out.extend_from_slice(&digests.qual.to_le_bytes());
    out.extend_from_slice(&(b.n_reads() as u32).to_le_bytes());
    for stream in [names_c, seq_c, qual_c] {
        // Stream lengths are stored as u32. The MAX_BLOCK_SEQ_BYTES row-group
        // budget keeps every compressed stream well under this, but guard the
        // cast so a future budget change can never silently truncate a length
        // and misframe the block on decode.
        let len = u32::try_from(stream.len())
            .map_err(|_| Error::Malformed("compressed stream exceeds u32 length"))?;
        out.extend_from_slice(&len.to_le_bytes());
        out.extend_from_slice(stream);
    }
    Ok(out)
}

/// Read blocks in batches of `batch`, invoking `f` on each batch.
pub(crate) fn for_each_block_batch<R: Read, F>(r: &mut R, batch: usize, mut f: F) -> Result<()>
where
    F: FnMut(&[Vec<u8>]) -> Result<()>,
{
    let mut block_index = 0u64;
    loop {
        let mut raw_blocks: Vec<Vec<u8>> = Vec::with_capacity(batch);
        for _ in 0..batch {
            match read_block(r, block_index)? {
                Some(block) => {
                    raw_blocks.push(block);
                    block_index += 1;
                }
                None => break,
            }
        }
        if raw_blocks.is_empty() {
            break;
        }
        let full = raw_blocks.len() == batch;
        f(&raw_blocks)?;
        if !full {
            break;
        }
    }
    Ok(())
}

/// Locate a block's three coded streams (names, sequence, quality) as absolute
/// `(offset, len, crc32c)` triples, given the payload bytes and the block frame's
/// absolute offset. Used to build the footer's per-stream projection index.
///
/// The payload is `[24 digests][4 n_reads] ([4 len][bytes])×3`, so each stream's
/// coded bytes start right after its length prefix; the returned offset is
/// absolute (past the frame head), and the CRC is over exactly those `len` bytes —
/// the same slice a remote client fetches, so it can verify a projected stream the
/// decoded-content digests can't cover.
pub(crate) fn payload_stream_locs(payload: &[u8], block_offset: u64) -> Result<[StreamLoc; 3]> {
    let base = block_offset + FRAME_HEAD_LEN as u64;
    let mut c = Cursor::new(payload);
    for _ in 0..3 {
        c.u64()?; // per-stream content digests (names, sequence, quality)
    }
    c.u32()?; // n_reads
    let mut locs = [StreamLoc::default(); 3];
    for loc in &mut locs {
        let len = c.u32()?;
        let off_in_payload = c.pos();
        let bytes = c.take(len as usize)?;
        *loc = StreamLoc {
            offset: base + off_in_payload as u64,
            len,
            crc: crc32c(bytes),
        };
    }
    Ok(locs)
}

/// Decoded block streams: (n_reads, names, per-read lengths, sequence, quality).
pub(crate) type BlockParts = (usize, Vec<Vec<u8>>, Vec<u32>, Vec<u8>, Vec<u8>);

/// Decode a block's streams and slice out each read's (name, seq, qual).
///
/// `shared_ref` is the archive's whole-file reference frame, threaded through to
/// the sequence decoder for [`SEQ_METHOD_OVERLAP_REF`] blocks (issue #168); `None`
/// for archives without a reference frame, which never contain such blocks.
pub(crate) fn decode_block_parts(
    buf: &[u8],
    shared_ref: Option<&fqxv_lroverlap::Reference>,
) -> Result<BlockParts> {
    let mut c = Cursor::new(buf);
    let expected = StreamDigests {
        names: c.u64()?,
        seq: c.u64()?,
        qual: c.u64()?,
    };
    let n_reads = c.u32()? as usize;
    // Slice out the three compressed streams (cheap, sequential), then decode
    // them concurrently — same rationale as the encode side.
    let (names_s, seq_s, qual_s) = (c.slice_u32()?, c.slice_u32()?, c.slice_u32()?);
    // A sequence-conditioned quality stream (long reads) must see the decoded
    // bases, so decode the sequence first and feed it in. Short-read quality is
    // sequence-blind, so keep decoding the two streams in parallel — the common
    // case pays nothing for this. Peek the quality header to tell them apart.
    let seq_first = fqxv_fqzcomp::needs_sequence(qual_s);
    let (names, (seq_r, qual_r)) = rayon::join(
        || fqxv_tokenizer::decode(names_s),
        || {
            if seq_first {
                let seq_r = decode_sequence_stream(seq_s, shared_ref);
                let qual_r = match &seq_r {
                    Ok((_, seq)) => fqxv_fqzcomp::decode_seq(qual_s, seq),
                    // Sequence decode failed; the block errors on `seq_r?` below.
                    // Decode quality with no sequence so the type lines up.
                    Err(_) => fqxv_fqzcomp::decode_seq(qual_s, &[]),
                };
                (seq_r, qual_r)
            } else {
                rayon::join(
                    || decode_sequence_stream(seq_s, shared_ref),
                    || fqxv_fqzcomp::decode_seq(qual_s, &[]),
                )
            }
        },
    );
    let names = names?;
    let (seq_lens, seq) = seq_r?;
    let (_qlens, qual) = qual_r?;
    if names.len() != n_reads || seq_lens.len() != n_reads {
        return Err(Error::Malformed("block stream length disagreement"));
    }
    // End-to-end check: each reconstructed stream must digest to the value the
    // encoder stored. A mismatch here (with the frame CRC intact) means a codec
    // decoded valid bytes into wrong output — the failure mode CRC cannot see —
    // and the per-stream digests name which stream regressed.
    let got = stream_digests(
        n_reads,
        names.iter().map(Vec::as_slice),
        &seq_lens,
        &seq,
        &qual,
    );
    for (ok, what) in [
        (got.names == expected.names, "block names digest"),
        (got.seq == expected.seq, "block sequence digest"),
        (got.qual == expected.qual, "block quality digest"),
    ] {
        if !ok {
            return Err(Error::Corrupt {
                what: what.to_string(),
            });
        }
    }
    Ok((n_reads, names, seq_lens, seq, qual))
}

pub(crate) fn write_record(out: &mut Vec<u8>, name: &[u8], seq: &[u8], qual: &[u8]) {
    out.push(b'@');
    out.extend_from_slice(name);
    out.push(b'\n');
    out.extend_from_slice(seq);
    out.extend_from_slice(b"\n+\n");
    out.extend_from_slice(qual);
    out.push(b'\n');
}

pub(crate) fn decode_block(
    buf: &[u8],
    shared_ref: Option<&fqxv_lroverlap::Reference>,
) -> Result<(u64, Vec<u8>)> {
    let (n_reads, names, lens, seq, qual) = decode_block_parts(buf, shared_ref)?;
    let mut out = Vec::with_capacity(seq.len() * 2 + qual.len());
    let mut off = 0usize;
    for i in 0..n_reads {
        let l = lens[i] as usize;
        // Checked slicing: a block whose per-read lengths overrun the decoded
        // sequence/quality buffers is malformed, not a reason to panic.
        let (s, q) = read_slices(&seq, &qual, off, l)?;
        write_record(&mut out, &names[i], s, q);
        off += l;
    }
    Ok((n_reads as u64, out))
}

/// Bounds-checked `(seq, qual)` slices for one read at `off..off+l`, erroring
/// instead of panicking when corrupted lengths overrun either buffer.
pub(crate) fn read_slices<'a>(
    seq: &'a [u8],
    qual: &'a [u8],
    off: usize,
    l: usize,
) -> Result<(&'a [u8], &'a [u8])> {
    let end = off
        .checked_add(l)
        .ok_or(Error::Malformed("read length overflow"))?;
    let s = seq.get(off..end).ok_or(Error::Malformed(
        "sequence shorter than declared read lengths",
    ))?;
    let q = qual.get(off..end).ok_or(Error::Malformed(
        "quality shorter than declared read lengths",
    ))?;
    Ok((s, q))
}

/// Split a grouped block into `g` FASTQ buffers by local read index mod `g`.
///
/// De-interleaving is a block-local `i % g`, which is correct only because every
/// block holds whole spots and starts on member 0 (enforced at encode by
/// [`block_ranges`] and the streaming loop). That invariant is not otherwise
/// recorded on disk, so a block whose read count is not a multiple of `g` — a
/// regression in the block splitter, or a crafted archive — would silently
/// misroute the trailing partial spot and shift every following block's members.
/// Reject it here so the failure is loud and localized rather than wrong output.
pub(crate) fn decode_block_group(
    buf: &[u8],
    g: usize,
    shared_ref: Option<&fqxv_lroverlap::Reference>,
) -> Result<(u64, Vec<Vec<u8>>)> {
    let (n_reads, names, lens, seq, qual) = decode_block_parts(buf, shared_ref)?;
    if !n_reads.is_multiple_of(g) {
        return Err(Error::Malformed(
            "grouped block read count is not a multiple of the group size",
        ));
    }
    let mut outs = vec![Vec::new(); g];
    let mut off = 0usize;
    for i in 0..n_reads {
        let l = lens[i] as usize;
        let (s, q) = read_slices(&seq, &qual, off, l)?;
        write_record(&mut outs[i % g], &names[i], s, q);
        off += l;
    }
    Ok((n_reads as u64, outs))
}

/// Read one framed, CRC-checked block, or `None` at the terminator / a clean EOF.
/// `index` names the block in any corruption error.
///
/// The frame is `[4 BLOCK_MAGIC][8 payload_len][4 crc32c(payload)][payload]`. A
/// zero-length block (magic + `len == 0`) is the terminator that separates the
/// block region from the footer, so a streaming (non-seekable) decoder stops here
/// without reading into the footer. A clean EOF (no bytes, or a partial marker) is
/// also treated as the end, which keeps truncated pre-footer streams decoding what
/// they can. The CRC is verified before the payload is handed to the entropy
/// decoders so corruption surfaces as a clean [`Error::Corrupt`] rather than
/// garbage output.
pub(crate) fn read_block<R: Read>(r: &mut R, index: u64) -> Result<Option<Vec<u8>>> {
    let mut magic = [0u8; BLOCK_MAGIC.len()];
    match r.read_exact(&mut magic) {
        Ok(()) => {}
        // No marker (or only a partial one) left: a clean end of the block region.
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e.into()),
    }
    if magic != BLOCK_MAGIC {
        return Err(Error::Corrupt {
            what: format!("block {index} (bad sync marker)"),
        });
    }
    let mut len = [0u8; 8];
    r.read_exact(&mut len).map_err(|_| Error::Truncated)?;
    let len = u64::from_le_bytes(len);
    if len == 0 {
        return Ok(None);
    }
    if len > MAX_BLOCK_PAYLOAD {
        return Err(Error::Malformed("block payload length exceeds the maximum"));
    }
    let mut crc = [0u8; CRC_LEN];
    r.read_exact(&mut crc).map_err(|_| Error::Truncated)?;
    let expected = u32::from_le_bytes(crc);
    // Read incrementally rather than `resize(len, 0)` + `read_exact`: `len` is
    // capped above, but a truncated stream claiming the full cap would still
    // zero-fill 2 GB up front before the short read is discovered (#142). `take`
    // bounds the read to the claim, so only the bytes actually present are
    // allocated; a short read surfaces as `Truncated`.
    let mut buf = Vec::new();
    let got = r
        .by_ref()
        .take(len)
        .read_to_end(&mut buf)
        .map_err(|_| Error::Truncated)?;
    if got as u64 != len {
        return Err(Error::Truncated);
    }
    if crc32c(&buf) != expected {
        return Err(Error::Corrupt {
            what: format!("block {index}"),
        });
    }
    Ok(Some(buf))
}
