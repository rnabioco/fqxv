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
        let decoded: Vec<Result<(u64, Vec<u8>)>> =
            pool.install(|| raw_blocks.par_iter().map(|b| decode_block(b)).collect());
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

/// A [`Write`] sink that folds a decoded interleaved-FASTQ stream into
/// [`ContentStats`] on the fly. Fed by [`content_stats`] as the decompressor's
/// output writer, it parses the four-line records incrementally — buffering only
/// the current line — so a multi-GB archive is summarized without ever holding
/// the decoded FASTQ in memory. Record lines cycle name → sequence → `+` →
/// quality; only sequence and quality are inspected.
#[derive(Default)]
pub(crate) struct StatsSink {
    stats: ContentStats,
    /// Bytes of the current line seen so far (newline excluded).
    line: Vec<u8>,
    /// Which line of the current record we are on: 0 name, 1 seq, 2 `+`, 3 qual.
    line_no: u8,
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
        self.line_no = (self.line_no + 1) % 4;
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
                .map(|b| decode_block_group(b, g))
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
/// Only the plain layout is recoverable this way; the globally-clustered reorder
/// layout is all-or-nothing (its streams are mutually dependent) and returns an
/// error directing the caller to plain [`decompress`]. If the footer itself is
/// unreadable (e.g. a truncated download), this returns the footer error; callers
/// wanting the intact prefix of such a file can fall back to streaming
/// [`decompress`], which decodes every whole block before the truncation.
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
    let footer = read_footer(&mut r)?;
    let mut rec = Recovery::default();
    let mut w = BufWriter::new(writer);
    for (i, &(off, read_count)) in footer.groups.iter().enumerate() {
        r.seek(SeekFrom::Start(off))?;
        // read_block bounds the length and verifies the block CRC; any failure —
        // bad CRC, truncation, or a decode error below — drops just this group.
        let outcome = match read_block(&mut r, i as u64) {
            Ok(Some(payload)) => pool.install(|| decode_block(&payload)),
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
    info!(
        recovered = rec.blocks_recovered,
        skipped = rec.blocks_skipped,
        reads_lost = rec.reads_lost,
        "recovery complete"
    );
    Ok(rec)
}
