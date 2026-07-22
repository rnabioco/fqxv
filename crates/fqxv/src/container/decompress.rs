//! Public decompression entry points and content-stats folding.

use super::*;
use rayon::prelude::*;
use tracing::{debug, info, instrument, trace, warn};

/// Decompress a `.fqxv` stream into interleaved FASTQ on `writer`.
///
/// For grouped archives this yields interleaved output — exactly what aligners
/// that accept interleaved paired reads want (`fqxv decompress x.fqxv | bwa mem -p`).
#[instrument(skip_all, fields(threads))]
pub fn decompress<R: Read, W: Write>(reader: R, writer: W, threads: usize) -> Result<Stats> {
    let pool = build_pool(threads)?;
    let batch = pool.current_num_threads().max(1);
    let mut r = BufReader::new(reader);
    let header = read_header(&mut r)?;
    // Whole-file globally-clustered reorder (both keep-order modes) uses a
    // distinct layout — two stream partitions and, with keep-order, a global
    // permutation — not the per-block loop below.
    if header.flags & FLAG_GLOBAL_REORDER != 0 {
        let keep_order = header.flags & FLAG_KEEP_ORDER != 0;
        let has_reference = header.flags & FLAG_GLOBAL_REFERENCE != 0;
        return decode_reordered_whole(
            r,
            writer,
            threads,
            keep_order,
            header.group_size,
            has_reference,
        );
    }
    // Whole-file shared reference frame (plain layout, issue #168): read it once,
    // right after the header, and thread it into every block's sequence decode.
    let reference = read_reference_frame(&mut r, header.flags)?;
    let reference = reference.as_ref();
    // Sequence-only archives store an empty quality stream and reconstruct FASTA.
    let no_quality = header.required_features & crate::feature::NO_QUALITY != 0;
    let mut w = BufWriter::new(writer);

    debug!(
        threads = pool.current_num_threads(),
        batch,
        backend = ?fqxv_rans::Backend::detect(),
        "decompress pool ready"
    );
    let mut stats = Stats::default();
    for_each_block_batch(&mut r, batch, |raw_blocks| {
        debug!(blocks = raw_blocks.len(), "decoding batch");
        let decoded: Vec<Result<(u64, Vec<u8>)>> = pool.install(|| {
            raw_blocks
                .par_iter()
                .map(|b| decode_block(b, reference, no_quality))
                .collect()
        });
        for d in decoded {
            let (reads, fastq) = d?;
            w.write_all(&fastq)?;
            trace!(reads, bytes = fastq.len(), "block decoded");
            stats.reads += reads;
            stats.blocks += 1;
            stats.out_bytes += fastq.len() as u64;
        }
        Ok(())
    })?;
    w.flush()?;
    Ok(stats)
}

/// A [`Write`] sink that folds a decoded interleaved-FASTQ (or FASTA) stream into
/// [`ContentStats`] on the fly. Fed by [`content_stats`] as the decompressor's
/// output writer, it parses records incrementally — buffering only the current
/// line — so a multi-GB archive is summarized without ever holding the decoded
/// output in memory. FASTQ record lines cycle name → sequence → `+` → quality;
/// a sequence-only (`no_quality`) archive emits FASTA (name → sequence), detected
/// from the first header byte (`>` vs `@`). Only sequence and quality are inspected.
#[derive(Default)]
pub(crate) struct StatsSink {
    stats: ContentStats,
    /// Bytes of the current line seen so far (newline excluded).
    line: Vec<u8>,
    /// Which line of the current record we are on: 0 name, 1 seq, 2 `+`, 3 qual.
    line_no: u8,
    /// Lines per record: 4 for FASTQ, 2 for FASTA. `0` until the first header line
    /// picks it from the leading byte.
    modulus: u8,
    /// Whether any read has been counted (so `min_len` starts from the first).
    seen: bool,
}

impl StatsSink {
    /// Fold one complete line (a trailing `\r` stripped) at the current record
    /// position, then advance to the next record line.
    fn commit_line(&mut self) {
        let mut line = &self.line[..];
        if line.last() == Some(&b'\r') {
            line = &line[..line.len() - 1];
        }
        if self.modulus == 0 {
            // First line of the stream is a record header: `>` marks FASTA
            // (two-line, sequence-only) records, `@` marks four-line FASTQ.
            self.modulus = if line.first() == Some(&b'>') { 2 } else { 4 };
        }
        match self.line_no {
            1 => {
                // Sequence: length spread + base composition.
                let len = line.len() as u32;
                if !self.seen || len < self.stats.min_len {
                    self.stats.min_len = len;
                }
                if !self.seen || len > self.stats.max_len {
                    self.stats.max_len = len;
                }
                self.seen = true;
                self.stats.reads += 1;
                self.stats.bases += line.len() as u64;
                for &b in line {
                    match b {
                        b'A' => self.stats.a += 1,
                        b'C' => self.stats.c += 1,
                        b'G' => self.stats.g += 1,
                        b'T' => self.stats.t += 1,
                        b'N' => self.stats.n += 1,
                        _ => self.stats.other += 1,
                    }
                }
            }
            3 => {
                // Quality: mean + per-Phred histogram.
                for &b in line {
                    let phred = (b.saturating_sub(33) as usize).min(QUAL_MAX - 1);
                    self.stats.qual_sum += phred as u64;
                    self.stats.qual_hist[phred] += 1;
                }
            }
            _ => {} // name / `+` line: ignored
        }
        self.line.clear();
        self.line_no = (self.line_no + 1) % self.modulus;
    }
}

impl Write for StatsSink {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        // Decompressor output arrives in arbitrary-sized chunks; reassemble lines
        // across chunk boundaries by buffering up to each newline.
        let mut rest = buf;
        while let Some(nl) = rest.iter().position(|&b| b == b'\n') {
            self.line.extend_from_slice(&rest[..nl]);
            self.commit_line();
            rest = &rest[nl + 1..];
        }
        self.line.extend_from_slice(rest);
        Ok(buf.len())
    }

    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

/// Compute [`ContentStats`] for an archive by decoding it in full.
///
/// This drives the ordinary [`decompress`] path into a folding sink, so it
/// handles every layout (plain, per-block and whole-file reorder, grouped) with
/// no codec-specific logic and is guaranteed consistent with what `decompress`
/// would actually emit. Cost is O(archive) — a real decode — unlike [`inspect`],
/// which is O(row groups). `threads` matches [`decompress`].
pub fn content_stats<R: Read>(reader: R, threads: usize) -> Result<ContentStats> {
    let mut sink = StatsSink::default();
    decompress(reader, &mut sink, threads)?;
    Ok(sink.stats)
}

/// Decompress a grouped archive, splitting reads back into `G` writers by their
/// per-spot member. `writers.len()` must equal the archive's group size.
#[instrument(skip_all, fields(threads, outputs = writers.len()))]
pub fn decompress_split<R: Read, W: Write>(
    reader: R,
    writers: &mut [W],
    threads: usize,
) -> Result<Stats> {
    let pool = build_pool(threads)?;
    let batch = pool.current_num_threads().max(1);
    let mut r = BufReader::new(reader);
    let header = read_header(&mut r)?;
    let g = header.group_size as usize;
    if writers.len() != g {
        return Err(Error::Malformed(
            "output count does not match archive group size",
        ));
    }
    if header.flags & FLAG_REORDERED != 0 {
        // Grouped global-reorder archives are always `keep_order`, so the
        // permutation restores the original spot interleaving and we can
        // de-interleave into the G writers. Single-end clustered-output archives
        // (no preserved order) cannot be split.
        let global = header.flags & FLAG_GLOBAL_REORDER != 0;
        let keep_order = header.flags & FLAG_KEEP_ORDER != 0;
        if global && keep_order {
            let has_reference = header.flags & FLAG_GLOBAL_REFERENCE != 0;
            return decode_reordered_split(r, writers, threads, g, has_reference);
        }
        return Err(Error::Malformed(
            "reordered archive without preserved order: use decompress, not split",
        ));
    }

    let reference = read_reference_frame(&mut r, header.flags)?;
    let reference = reference.as_ref();
    let no_quality = header.required_features & crate::feature::NO_QUALITY != 0;
    debug!(
        threads = pool.current_num_threads(),
        batch,
        group_size = g,
        backend = ?fqxv_rans::Backend::detect(),
        "decompress-split pool ready"
    );
    let mut stats = Stats::default();
    for_each_block_batch(&mut r, batch, |raw_blocks| {
        debug!(blocks = raw_blocks.len(), "decoding batch");
        let decoded: Vec<Result<(u64, Vec<Vec<u8>>)>> = pool.install(|| {
            raw_blocks
                .par_iter()
                .map(|b| decode_block_group(b, g, reference, no_quality))
                .collect()
        });
        for d in decoded {
            let (reads, parts) = d?;
            for (w, part) in writers.iter_mut().zip(&parts) {
                w.write_all(part)?;
            }
            stats.reads += reads;
            stats.blocks += 1;
            stats.out_bytes += parts.iter().map(|p| p.len() as u64).sum::<u64>();
        }
        Ok(())
    })?;
    for w in writers.iter_mut() {
        w.flush()?;
    }
    Ok(stats)
}

/// Outcome of a best-effort [`decompress_recover`] run.
#[derive(Debug, Default, Clone)]
pub struct Recovery {
    /// Reads and bytes actually recovered (in the intact blocks).
    pub stats: Stats,
    /// Blocks decoded successfully.
    pub blocks_recovered: u64,
    /// Blocks skipped because their CRC failed or they would not decode.
    pub blocks_skipped: u64,
    /// Reads lost to the skipped blocks (from the footer's per-group counts).
    pub reads_lost: u64,
}

/// Decompress as much of a corrupted archive as possible, skipping bad blocks.
///
/// Blocks are independent, and the footer's row-group index gives each block's
/// absolute offset, so a corrupt block can be skipped by seeking straight to the
/// next one — one bad byte costs a single row group, not the whole archive. Each
/// skipped block is logged (with its lost read count) and reported in
/// [`Recovery`]. Output is interleaved FASTQ (as [`decompress`]).
///
/// If the footer itself is unreadable — the common truncated-tail case, which also
/// loses the trailing blocks — recovery falls back to **scanning for the
/// per-block `BLOCK_MAGIC` sync marker**: it needs no index and resynchronizes past
/// a corrupt length prefix or a bad block, so every intact block is still
/// recovered (loss counts aren't available without the footer, so `reads_lost` is
/// 0 in that mode). Only the plain layout is recoverable; the globally-clustered
/// reorder layout is all-or-nothing (its streams are mutually dependent) and
/// returns an error directing the caller to plain [`decompress`].
pub fn decompress_recover<R: Read + Seek, W: Write>(
    reader: R,
    writer: W,
    threads: usize,
) -> Result<Recovery> {
    let pool = build_pool(threads)?;
    let mut r = BufReader::new(reader);
    let header = read_header(&mut r)?;
    if header.flags & FLAG_GLOBAL_REORDER != 0 {
        return Err(Error::Malformed(
            "recover supports only the plain layout; reordered archives decode all-or-nothing",
        ));
    }
    // Read the whole-file shared reference frame (if any) before seeking around, so
    // reference-coded blocks can be recovered. A corrupt frame fails closed here —
    // its blocks then simply won't decode and are skipped like any other bad block.
    let reference = read_reference_frame(&mut r, header.flags)?;
    let reference = reference.as_ref();
    let no_quality = header.required_features & crate::feature::NO_QUALITY != 0;
    // Prefer the footer's row-group index — it carries per-group read counts, so
    // losses can be tallied exactly. If the footer is unreadable (the common
    // truncated-tail case, which also loses the trailing blocks), fall back to
    // scanning for block sync markers: that needs no index and resynchronizes past
    // a corrupt length prefix or a bad block.
    let rec = match read_footer(&mut r) {
        Ok(footer) => recover_via_footer(&mut r, &pool, &footer, writer, reference, no_quality)?,
        Err(footer_err) => {
            debug!(error = %footer_err, "footer unreadable; scanning for block markers");
            recover_via_scan(&mut r, &pool, writer, reference, no_quality)?
        }
    };
    info!(
        recovered = rec.blocks_recovered,
        skipped = rec.blocks_skipped,
        reads_lost = rec.reads_lost,
        "recovery complete"
    );
    Ok(rec)
}

/// Footer-driven recovery: seek to each row group's block via the index and skip
/// any that won't decode. Per-group read counts make loss accounting exact.
fn recover_via_footer<R: Read + Seek, W: Write>(
    r: &mut BufReader<R>,
    pool: &rayon::ThreadPool,
    footer: &Footer,
    writer: W,
    reference: Option<&fqxv_lroverlap::Reference>,
    no_quality: bool,
) -> Result<Recovery> {
    let mut rec = Recovery::default();
    let mut w = BufWriter::new(writer);
    for (i, &(off, read_count)) in footer.groups.iter().enumerate() {
        r.seek(SeekFrom::Start(off))?;
        // read_block checks the marker, bounds the length, and verifies the CRC;
        // any failure — bad marker/CRC, truncation, or a decode error — drops just
        // this group.
        let outcome = match read_block(r, i as u64) {
            Ok(Some(payload)) => pool.install(|| decode_block(&payload, reference, no_quality)),
            Ok(None) => Err(Error::Malformed(
                "row-group offset points at the terminator",
            )),
            Err(e) => Err(e),
        };
        match outcome {
            Ok((reads, fastq)) => {
                w.write_all(&fastq)?;
                rec.stats.reads += reads;
                rec.stats.blocks += 1;
                rec.stats.out_bytes += fastq.len() as u64;
                rec.blocks_recovered += 1;
            }
            Err(e) => {
                warn!(block = i, off, reads_lost = read_count, error = %e, "skipping corrupt block");
                rec.blocks_skipped += 1;
                rec.reads_lost += read_count as u64;
            }
        }
    }
    w.flush()?;
    Ok(rec)
}

/// Footer-independent recovery: scan the block region for `BLOCK_MAGIC` and decode
/// every frame that validates, resynchronizing past corruption or a lost index — so
/// a truncated tail (which takes the footer with it) or a corrupt length prefix no
/// longer strands the blocks after it.
///
/// Reads the block region (header..EOF) into memory to scan it; recovery is an
/// exceptional break-glass operation, so the transient copy is worth the
/// simplicity. Without the footer there are no per-group read counts, so
/// `reads_lost` stays 0 — the recovered counts are still reported, and
/// `blocks_skipped` counts frames whose CRC passed but which would not decode.
fn recover_via_scan<R: Read + Seek, W: Write>(
    r: &mut BufReader<R>,
    pool: &rayon::ThreadPool,
    writer: W,
    reference: Option<&fqxv_lroverlap::Reference>,
    no_quality: bool,
) -> Result<Recovery> {
    r.seek(SeekFrom::Start(HEADER_LEN as u64))?;
    let mut region = Vec::new();
    r.read_to_end(&mut region)?;

    let mut rec = Recovery::default();
    let mut w = BufWriter::new(writer);
    let mut pos = 0usize;
    let mut idx = 0u64;
    while pos < region.len() {
        let Some(rel) = find_block_marker(&region[pos..]) else {
            break;
        };
        let m = pos + rel;
        match read_frame_at(&region, m) {
            Some((len, payload)) => {
                match pool.install(|| decode_block(payload, reference, no_quality)) {
                    Ok((reads, fastq)) => {
                        w.write_all(&fastq)?;
                        rec.stats.reads += reads;
                        rec.stats.blocks += 1;
                        rec.stats.out_bytes += fastq.len() as u64;
                        rec.blocks_recovered += 1;
                    }
                    Err(e) => {
                        warn!(block = idx, off = m, error = %e, "scan: frame CRC ok but won't decode");
                        rec.blocks_skipped += 1;
                    }
                }
                idx += 1;
                pos = m + FRAME_HEAD_LEN + len; // past this frame (validated to fit)
            }
            // Terminator, bad length, short frame, or a false-positive marker:
            // resync one byte on and keep scanning.
            None => pos = m + 1,
        }
    }
    w.flush()?;
    Ok(rec)
}

/// First `BLOCK_MAGIC` occurrence in `buf`, if any.
fn find_block_marker(buf: &[u8]) -> Option<usize> {
    buf.windows(BLOCK_MAGIC.len())
        .position(|w| w == BLOCK_MAGIC)
}

/// If a valid data-block frame starts at `region[m]`, return its payload length
/// and the CRC-verified payload slice. `None` for the terminator (`len == 0`), a
/// bad/oversized length, a frame that runs past the buffer, or a CRC mismatch —
/// i.e. a false-positive marker or a corrupt frame.
fn read_frame_at(region: &[u8], m: usize) -> Option<(usize, &[u8])> {
    let head = region.get(m..m + FRAME_HEAD_LEN)?;
    if head[..BLOCK_MAGIC.len()] != BLOCK_MAGIC {
        return None;
    }
    let len = u64::from_le_bytes(head[4..12].try_into().unwrap());
    if len == 0 || len > MAX_BLOCK_PAYLOAD {
        return None;
    }
    let len = len as usize;
    let expected = u32::from_le_bytes(head[12..16].try_into().unwrap());
    let payload = region.get(m + FRAME_HEAD_LEN..m + FRAME_HEAD_LEN + len)?;
    (crc32c(payload) == expected).then_some((len, payload))
}
