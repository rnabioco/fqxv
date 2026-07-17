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

/// xxh3-64 over a block's *decoded canonical form*: the exact (name, sequence,
/// quality) bytes `decompress` reconstructs, structured so no byte can silently
/// cross a name/seq/quality boundary. Computed identically on the encode side
/// (from the post-`QualityBinning` quality — the values actually stored) and the
/// decode side (from the reconstructed streams). A mismatch means some codec
/// round-tripped this block into wrong-but-in-bounds output — corruption the
/// per-payload CRC cannot catch, because the stored bytes were never altered.
///
/// The digest is over the *stored* (post-binning) form, not the original input,
/// so a lossy archive verifies against what it emits, not against data it never
/// promised to reproduce. Name/quality lengths are folded in explicitly to pin
/// the per-read boundaries; `qual` shares each read's `lens[i]` with `seq`.
pub(crate) fn content_digest<'a>(
    n_reads: usize,
    names: impl Iterator<Item = &'a [u8]>,
    lens: &[u32],
    seq: &[u8],
    qual: &[u8],
) -> u64 {
    let mut h = Xxh3::new();
    h.update(&(n_reads as u64).to_le_bytes());
    for name in names {
        h.update(&(name.len() as u32).to_le_bytes());
        h.update(name);
    }
    for &l in lens {
        h.update(&l.to_le_bytes());
    }
    h.update(seq);
    h.update(qual);
    h.digest()
}

/// Code one non-reorder block: names (tokenizer), sequence (order-k), and quality
/// (fqzcomp), each length-prefixed, behind a leading [`content_digest`]. Reorder
/// uses the whole-file path instead.
pub(crate) fn compress_block(b: &RawBlock, params: &Params) -> Result<Vec<u8>> {
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
                || {
                    fqxv_seq::encode_hashed(
                        &b.lens,
                        &b.seq,
                        params.seq_order as usize,
                        params.seq_hash_order as usize,
                        u32::from(params.seq_hash_bits),
                    )
                },
                || fqxv_fqzcomp::encode(&b.lens, &b.qual, params.quality_binning),
            )
        },
    );
    let (names_c, seq_c, qual_c) = (names_c?, seq_c?, qual_c?);

    // End-to-end round-trip check: digest the block's decoded content (post-binning
    // quality, so lossy archives verify against what they emit) and store it at the
    // head of the payload. Lossless is the common case and borrows without a copy.
    let binned: Cow<[u8]> = match params.quality_binning {
        QualityBinning::Lossless => Cow::Borrowed(&b.qual),
        binning => Cow::Owned(b.qual.iter().map(|&q| binning.apply(q)).collect()),
    };
    let digest = content_digest(
        b.n_reads(),
        b.header_refs().into_iter(),
        &b.lens,
        &b.seq,
        &binned,
    );

    let mut out = Vec::with_capacity(DIGEST_LEN + 16 + names_c.len() + seq_c.len() + qual_c.len());
    out.extend_from_slice(&digest.to_le_bytes());
    out.extend_from_slice(&(b.n_reads() as u32).to_le_bytes());
    for stream in [&names_c, &seq_c, &qual_c] {
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
/// The payload is `[8 digest][4 n_reads] ([4 len][bytes])×3`, so each stream's
/// coded bytes start right after its length prefix; the returned offset is
/// absolute (past the frame head), and the CRC is over exactly those `len` bytes —
/// the same slice a remote client fetches, so it can verify a projected stream the
/// joint content digest can't cover.
pub(crate) fn payload_stream_locs(payload: &[u8], block_offset: u64) -> Result<[StreamLoc; 3]> {
    let base = block_offset + FRAME_HEAD_LEN as u64;
    let mut c = Cursor::new(payload);
    c.u64()?; // content digest
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
pub(crate) fn decode_block_parts(buf: &[u8]) -> Result<BlockParts> {
    let mut c = Cursor::new(buf);
    let expected_digest = c.u64()?;
    let n_reads = c.u32()? as usize;
    // Slice out the three compressed streams (cheap, sequential), then decode
    // them concurrently — same rationale as the encode side.
    let (names_s, seq_s, qual_s) = (c.slice_u32()?, c.slice_u32()?, c.slice_u32()?);
    let (names, (seq_r, qual_r)) = rayon::join(
        || fqxv_tokenizer::decode(names_s),
        || rayon::join(|| fqxv_seq::decode(seq_s), || fqxv_fqzcomp::decode(qual_s)),
    );
    let names = names?;
    let (seq_lens, seq) = seq_r?;
    let (_qlens, qual) = qual_r?;
    if names.len() != n_reads || seq_lens.len() != n_reads {
        return Err(Error::Malformed("block stream length disagreement"));
    }
    // End-to-end check: the reconstructed content must digest to the value the
    // encoder stored. A mismatch here (with the frame CRC intact) means a codec
    // decoded valid bytes into wrong output — the failure mode CRC cannot see.
    let digest = content_digest(
        n_reads,
        names.iter().map(Vec::as_slice),
        &seq_lens,
        &seq,
        &qual,
    );
    if digest != expected_digest {
        return Err(Error::Corrupt {
            what: "block content digest".to_string(),
        });
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

pub(crate) fn decode_block(buf: &[u8]) -> Result<(u64, Vec<u8>)> {
    let (n_reads, names, lens, seq, qual) = decode_block_parts(buf)?;
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
pub(crate) fn decode_block_group(buf: &[u8], g: usize) -> Result<(u64, Vec<Vec<u8>>)> {
    let (n_reads, names, lens, seq, qual) = decode_block_parts(buf)?;
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
    // Fallible allocation: a corrupted-but-in-range length still shouldn't abort
    // the process with an allocation failure.
    let mut buf = Vec::new();
    buf.try_reserve_exact(len as usize)
        .map_err(|_| Error::Malformed("block payload too large to allocate"))?;
    buf.resize(len as usize, 0);
    r.read_exact(&mut buf).map_err(|_| Error::Truncated)?;
    if crc32c(&buf) != expected {
        return Err(Error::Corrupt {
            what: format!("block {index}"),
        });
    }
    Ok(Some(buf))
}
