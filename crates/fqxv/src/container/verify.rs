//! Integrity checking: whole-file CRC, per-block scans, and structured reports.

use super::*;
use rayon::prelude::*;

/// Verify an archive's integrity without materializing the decoded FASTQ.
///
/// For the plain layout this checks the footer's own CRC, then re-hashes the
/// archive prefix and compares against the stored whole-file CRC-32C — a single
/// linear pass that catches any corruption in the header, any block payload,
/// framing, or the index. For the globally-clustered reorder layout (which has no
/// footer) it decodes into a sink, so every frame CRC and cross-stream length
/// check is exercised. Returns `Ok(())` iff the archive is intact.
pub fn verify<R: Read + Seek>(reader: R, threads: usize) -> Result<()> {
    let mut r = BufReader::new(reader);
    let header = read_header(&mut r)?;
    if header.flags & FLAG_GLOBAL_REORDER != 0 {
        // No footer/whole-file digest here; decoding drives every frame CRC.
        let keep_order = header.flags & FLAG_KEEP_ORDER != 0;
        let has_reference = header.flags & FLAG_GLOBAL_REFERENCE != 0;
        decode_reordered_whole(
            r,
            io::sink(),
            threads,
            keep_order,
            header.group_size,
            has_reference,
        )?;
        return Ok(());
    }
    let footer = read_footer(&mut r)?;
    r.seek(SeekFrom::Start(0))?;
    if verify_whole_file_crc(&mut r, footer.covered_len, threads)? != footer.whole_file_crc {
        return Err(Error::Corrupt {
            what: "archive (whole-file crc)".to_string(),
        });
    }
    Ok(())
}

/// Full post-write verification: re-read a just-written archive and prove it is
/// both intact on disk *and* decodes cleanly, returning the decoded read count.
///
/// This is the check to run after `compress` before trusting (or deleting) the
/// source: it fully decodes the archive into a sink so every per-block content
/// digest and cross-stream length check is exercised — catching a codec bug that
/// produced a CRC-valid-but-wrong archive, which the structural [`verify`] alone
/// cannot. For the plain layout it *also* runs [`verify`] first, because the
/// whole-file CRC covers the header, inter-block framing, and footer index —
/// bytes a streaming decode never reads. The globally-clustered reorder layout has
/// no whole-file CRC, but a single decode already drives every frame CRC and the
/// trailing output digest, so it is not decoded twice.
///
/// It does not compare against the original FASTQ bytes, so it cannot catch a
/// *parse*-level bug (the content digest is computed from the parsed block); it
/// proves the archive round-trips to the content that was encoded. Being
/// page-cache-coherent, it validates the bytes the encoder emitted rather than the
/// storage medium — latent bit-rot is what a later [`verify`] catches.
///
/// `threads` sizes the rayon pool for the whole-file CRC and the decode exactly as
/// it does for [`decompress`] — `0` means all cores. This is what makes a CLI
/// `--verify` honor the command's `--threads`: the check runs under the same worker
/// budget as the compression it follows rather than silently grabbing every core.
pub fn verify_roundtrip<R: Read + Seek>(mut reader: R, threads: usize) -> Result<u64> {
    let header = read_header(&mut BufReader::new(&mut reader))?;
    reader.seek(SeekFrom::Start(0))?;

    if header.flags & FLAG_GLOBAL_REORDER != 0 {
        // One decode is both the structural and the content check here.
        return Ok(decompress(reader, io::sink(), threads)?.reads);
    }

    // Plain layout: whole-file CRC (structure/framing) then decode (content).
    verify(&mut reader, threads)?;
    reader.seek(SeekFrom::Start(0))?;
    Ok(decompress(reader, io::sink(), threads)?.reads)
}

/// The authoritative number of reads the archive should contain, for detecting a
/// truncated file on decompress.
///
/// For the plain and per-block-reorder layouts this reads the footer's
/// `total_reads` field **with no forward-scan fallback**: a file that lost its
/// trailing blocks also lost the footer/EOF trailer, so [`read_footer`] fails and
/// the truncation surfaces as an error rather than a short, silent success.
/// (Contrast [`inspect`], which deliberately falls back to counting the surviving
/// blocks so a partial file still reports what it holds — the wrong tool for a
/// completeness check.)
///
/// Returns `Ok(None)` for the globally-clustered reorder layout, which carries no
/// footer but is already truncation-safe: its trailing whole-output digest frame is
/// the last thing in the file, so a lost tail fails that digest on decode.
pub fn expected_reads<R: Read + Seek>(reader: R) -> Result<Option<u64>> {
    let mut r = BufReader::new(reader);
    let header = read_header(&mut r)?;
    if header.flags & FLAG_GLOBAL_REORDER != 0 {
        return Ok(None);
    }
    Ok(Some(read_footer(&mut r)?.total_reads))
}

/// CRC-32C of the first `covered` bytes of `r`, computed in parallel.
///
/// The stream is read serially (a single `Read` can't be shared across threads),
/// but the CPU-bound checksum is not: bytes are pulled in batches of fixed-size
/// chunks, every chunk in a batch is hashed on the rayon pool at once, and the
/// per-chunk CRCs are folded back together in order with [`crc32c_combine`]. The
/// result is bit-identical to a single-pass hash for any thread count or
/// chunking, so this stays exactly as strong a check as the serial version — it
/// just stops leaving cores idle while the table-driven CRC runs.
pub(crate) fn verify_whole_file_crc<R: Read>(
    r: &mut R,
    covered: u64,
    threads: usize,
) -> Result<u32> {
    /// Per-chunk CRC granularity: large enough that combine/dispatch overhead is
    /// negligible, small enough to keep a batch's buffers bounded in memory.
    const CHUNK: usize = 1 << 20; // 1 MiB

    // Small archives: one serial pass, no thread-pool spin-up to amortize.
    if covered <= CHUNK as u64 {
        let mut buf = vec![0u8; covered as usize];
        r.read_exact(&mut buf)?;
        return Ok(crc32c(&buf));
    }

    let pool = build_pool(threads)?;
    let batch = pool.current_num_threads().max(1);
    let mut running = crc32c(&[]); // CRC of the empty prefix (0x0000_0000)
    let mut remaining = covered;
    let mut bufs: Vec<Vec<u8>> = Vec::with_capacity(batch);
    while remaining > 0 {
        bufs.clear();
        while bufs.len() < batch && remaining > 0 {
            let want = remaining.min(CHUNK as u64) as usize;
            let mut buf = vec![0u8; want];
            r.read_exact(&mut buf)?;
            remaining -= want as u64;
            bufs.push(buf);
        }
        let crcs: Vec<(u32, usize)> =
            pool.install(|| bufs.par_iter().map(|b| (crc32c(b), b.len())).collect());
        for (crc, len) in crcs {
            running = crc32c_combine(running, crc, len as u64);
        }
    }
    Ok(running)
}

/// Quick integrity check: verify every block's stored CRC-32C via the footer's
/// row-group index, in parallel, without recomputing the whole-file digest.
///
/// This is a deliberately *weaker* check than [`verify`]: it covers each block's
/// coded payload but **not** the 10-byte header, the footer index itself, or the
/// inter-block framing (`[len][crc]`) bytes — corruption confined to those
/// regions slips past. In exchange the blocks are read straight from their footer
/// offsets with positioned reads and checked concurrently, so on a parallel
/// filesystem the many independent reads can outrun [`verify`]'s single serial
/// scan of the whole file. It takes a [`File`] (rather than a generic reader)
/// because that concurrency relies on offset-addressed reads.
///
/// The globally-clustered reorder layout has no per-block footer index, so this
/// transparently falls back to the full decode-driven [`verify`] for that layout
/// (as does any platform without positioned-read support).
pub fn verify_quick(file: &File, threads: usize) -> Result<()> {
    let mut r = BufReader::new(file);
    let header = read_header(&mut r)?;
    if header.flags & FLAG_GLOBAL_REORDER != 0 {
        // No per-block index to check — decoding is the only integrity path.
        r.seek(SeekFrom::Start(0))?;
        return verify(r, threads);
    }
    let footer = read_footer(&mut r)?;
    quick_check_blocks(file, &footer, threads)
}

/// Verify each block's stored CRC in parallel via positioned reads (Unix).
#[cfg(unix)]
pub(crate) fn quick_check_blocks(file: &File, footer: &Footer, threads: usize) -> Result<()> {
    use std::os::unix::fs::FileExt;

    let pool = build_pool(threads)?;
    pool.install(|| {
        footer
            .groups
            .par_iter()
            .enumerate()
            .try_for_each(|(i, &(off, _read_count))| {
                // Block frame on disk: [4 BLOCK_MAGIC][8 payload_len][4 crc32c][payload].
                let mut frame = [0u8; FRAME_HEAD_LEN];
                file.read_exact_at(&mut frame, off)
                    .map_err(|_| Error::Truncated)?;
                if frame[..BLOCK_MAGIC.len()] != BLOCK_MAGIC {
                    return Err(Error::Corrupt {
                        what: format!("block {i} (bad sync marker)"),
                    });
                }
                let len = u64::from_le_bytes(frame[4..12].try_into().unwrap());
                if len == 0 {
                    return Err(Error::Malformed(
                        "row-group offset points at the terminator",
                    ));
                }
                if len > MAX_BLOCK_PAYLOAD {
                    return Err(Error::Malformed("block payload length exceeds the maximum"));
                }
                let expected = u32::from_le_bytes(frame[12..].try_into().unwrap());
                // Fallible allocation: a corrupt-but-in-range length must not abort
                // the process (mirrors `read_block`).
                let mut buf = Vec::new();
                buf.try_reserve_exact(len as usize)
                    .map_err(|_| Error::Malformed("block payload too large to allocate"))?;
                buf.resize(len as usize, 0);
                file.read_exact_at(&mut buf, off + FRAME_HEAD_LEN as u64)
                    .map_err(|_| Error::Truncated)?;
                if crc32c(&buf) != expected {
                    return Err(Error::Corrupt {
                        what: format!("block {i}"),
                    });
                }
                Ok(())
            })
    })
}

/// Platforms without positioned reads can't read blocks concurrently from one
/// handle, so quick verification falls back to the full [`verify`] scan.
#[cfg(not(unix))]
pub(crate) fn quick_check_blocks(file: &File, _footer: &Footer, threads: usize) -> Result<()> {
    verify(BufReader::new(file), threads)
}

/// One named integrity check in a [`VerifyReport`].
#[derive(Debug, Clone)]
pub struct VerifyCheck {
    /// What was checked (e.g. `"header"`, `"footer"`, `"block CRCs"`).
    pub name: String,
    /// Whether it passed.
    pub ok: bool,
    /// Human-readable context (counts, the failure reason, or empty).
    pub detail: String,
}

/// Structured result of [`verify_report`]: the individual integrity checks run
/// against an archive, plus the indices of any blocks whose CRC failed. The
/// overall verdict is [`VerifyReport::passed`].
#[derive(Debug, Clone, Default)]
pub struct VerifyReport {
    /// The checks performed, in report order.
    pub checks: Vec<VerifyCheck>,
    /// Total blocks the archive declares (0 for the reorder layout, which has no
    /// footer index).
    pub blocks_total: u64,
    /// Indices of blocks whose CRC failed (plain layout). In the default check
    /// this is populated only when the whole-file CRC fails and the archive is
    /// scanned to localize the damage; in `quick` mode it is always the scan.
    pub failed_blocks: Vec<u64>,
}

impl VerifyReport {
    /// True iff every check passed (the archive is intact).
    pub fn passed(&self) -> bool {
        self.checks.iter().all(|c| c.ok)
    }

    fn push(&mut self, name: &str, ok: bool, detail: impl Into<String>) {
        self.checks.push(VerifyCheck {
            name: name.to_string(),
            ok,
            detail: detail.into(),
        });
    }
}

/// Verify an archive's integrity and return a per-check [`VerifyReport`] — the
/// structured form behind the CLI's `verify` table (see [`verify`] for the plain
/// boolean and [`verify_quick`] for the block-only check).
///
/// With `quick == false` (the default) it runs the same checks as [`verify`] —
/// header, footer, and the parallel whole-file CRC — and, only if that CRC fails,
/// scans the archive block by block to name the damaged ones in
/// [`VerifyReport::failed_blocks`]. With `quick == true` it runs the weaker
/// per-block check of [`verify_quick`] (each block's stored CRC, in parallel via
/// positioned reads) and skips the whole-file digest. The globally-clustered
/// reorder layout has no footer, so either mode verifies it by decoding into a
/// sink (one `streams (decode)` check).
///
/// Returns `Err` only when the input is not a readable fqxv container (bad
/// magic/version); a recognized-but-corrupt archive comes back as a report whose
/// [`passed`](VerifyReport::passed) is `false`.
pub fn verify_report(file: &File, quick: bool, threads: usize) -> Result<VerifyReport> {
    let mut r = BufReader::new(file);
    let header = read_header(&mut r)?;
    let mut report = VerifyReport::default();

    if header.flags & FLAG_GLOBAL_REORDER != 0 {
        report.push(
            "header",
            true,
            format!(
                "format v{}.{}, global-cluster reorder layout",
                header.major, header.minor
            ),
        );
        // No footer/per-block index; decoding drives every frame CRC.
        r.seek(SeekFrom::Start(0))?;
        match verify(r, threads) {
            Ok(()) => report.push("streams (decode)", true, "all frame CRCs verified"),
            Err(e) => report.push("streams (decode)", false, e.to_string()),
        }
        return Ok(report);
    }

    report.push(
        "header",
        true,
        format!("format v{}.{}, plain layout", header.major, header.minor),
    );

    // Footer (read_footer verifies the footer's own CRC).
    let footer = match read_footer(&mut r) {
        Ok(footer) => footer,
        Err(e) => {
            report.push("footer", false, e.to_string());
            return Ok(report);
        }
    };
    report.blocks_total = footer.groups.len() as u64;
    report.push(
        "footer",
        true,
        format!(
            "{} blocks, {} reads",
            report.blocks_total, footer.total_reads
        ),
    );

    if quick {
        // Weaker, faster: only the per-block payload CRCs (parallel).
        report.failed_blocks = scan_failed_blocks(file, &footer, threads);
        let ok_blocks = report.blocks_total - report.failed_blocks.len() as u64;
        let ok = report.failed_blocks.is_empty();
        let detail = if ok {
            format!("{ok_blocks}/{} intact", report.blocks_total)
        } else {
            format!(
                "{ok_blocks}/{} intact; failed: {}",
                report.blocks_total,
                summarize_indices(&report.failed_blocks)
            )
        };
        report.push("block CRCs", ok, detail);
        return Ok(report);
    }

    // Default: parallel whole-file CRC; localize block by block only on failure.
    r.seek(SeekFrom::Start(0))?;
    let crc_ok =
        verify_whole_file_crc(&mut r, footer.covered_len, threads)? == footer.whole_file_crc;
    if crc_ok {
        report.push(
            "block CRCs",
            true,
            format!("{}/{} intact", report.blocks_total, report.blocks_total),
        );
        report.push("whole-file CRC", true, "");
    } else {
        report.failed_blocks = scan_failed_blocks(file, &footer, threads);
        let ok_blocks = report.blocks_total - report.failed_blocks.len() as u64;
        let detail = if report.failed_blocks.is_empty() {
            format!(
                "{ok_blocks}/{} intact (damage is in the header, framing, or index)",
                report.blocks_total
            )
        } else {
            format!(
                "{ok_blocks}/{} intact; failed: {}",
                report.blocks_total,
                summarize_indices(&report.failed_blocks)
            )
        };
        report.push("block CRCs", report.failed_blocks.is_empty(), detail);
        report.push("whole-file CRC", false, "digest mismatch");
    }

    Ok(report)
}

/// Scan every block's stored CRC and return the indices that fail, in ascending
/// order. Parallel positioned reads on Unix (like [`quick_check_blocks`], but
/// collecting *all* failures rather than stopping at the first); a serial
/// footer-driven scan elsewhere.
#[cfg(unix)]
pub(crate) fn scan_failed_blocks(file: &File, footer: &Footer, threads: usize) -> Vec<u64> {
    let pool = match build_pool(threads) {
        Ok(pool) => pool,
        Err(_) => return scan_failed_blocks_serial(file, footer),
    };
    pool.install(|| {
        footer
            .groups
            .par_iter()
            .enumerate()
            .filter_map(|(i, &(off, _read_count))| {
                let ok = block_crc_ok(file, off).unwrap_or(false);
                (!ok).then_some(i as u64)
            })
            .collect()
    })
}

/// True if the block frame at `off` reads back with a matching stored CRC.
#[cfg(unix)]
pub(crate) fn block_crc_ok(file: &File, off: u64) -> Result<bool> {
    use std::os::unix::fs::FileExt;

    let mut frame = [0u8; FRAME_HEAD_LEN];
    file.read_exact_at(&mut frame, off)
        .map_err(|_| Error::Truncated)?;
    if frame[..BLOCK_MAGIC.len()] != BLOCK_MAGIC {
        return Ok(false);
    }
    let len = u64::from_le_bytes(frame[4..12].try_into().unwrap());
    if len == 0 || len > MAX_BLOCK_PAYLOAD {
        return Ok(false);
    }
    let expected = u32::from_le_bytes(frame[12..].try_into().unwrap());
    let mut buf = Vec::new();
    if buf.try_reserve_exact(len as usize).is_err() {
        return Ok(false);
    }
    buf.resize(len as usize, 0);
    file.read_exact_at(&mut buf, off + FRAME_HEAD_LEN as u64)
        .map_err(|_| Error::Truncated)?;
    Ok(crc32c(&buf) == expected)
}

/// Serial footer-driven block scan (the fallback / non-Unix path).
pub(crate) fn scan_failed_blocks_serial(file: &File, footer: &Footer) -> Vec<u64> {
    let mut r = BufReader::new(file);
    let mut failed = Vec::new();
    for (i, &(off, _read_count)) in footer.groups.iter().enumerate() {
        let ok = r.seek(SeekFrom::Start(off)).is_ok()
            && matches!(read_block(&mut r, i as u64), Ok(Some(_)));
        if !ok {
            failed.push(i as u64);
        }
    }
    failed
}

#[cfg(not(unix))]
pub(crate) fn scan_failed_blocks(file: &File, footer: &Footer, _threads: usize) -> Vec<u64> {
    scan_failed_blocks_serial(file, footer)
}

/// Render a list of failed block indices compactly, capping the enumeration so a
/// pathological archive can't print thousands of numbers (`"3, 91, … (+15)"`).
pub(crate) fn summarize_indices(indices: &[u64]) -> String {
    const CAP: usize = 12;
    let shown: Vec<String> = indices.iter().take(CAP).map(u64::to_string).collect();
    if indices.len() > CAP {
        format!("{}, … (+{})", shown.join(", "), indices.len() - CAP)
    } else {
        shown.join(", ")
    }
}
