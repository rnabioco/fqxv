//! Parallel FASTQ parsing into per-chunk record metadata and row-group ranges.

use super::*;
use rayon::prelude::*;

// --- parallel FASTQ parsing --------------------------------------------------

/// One record's location: its normalized header lives in the owning chunk's
/// header arena (`[prev.hdr_end .. hdr_end)`); its sequence and quality bytes are
/// contiguous ranges of the input buffer (CR/LF already excluded).
#[derive(Clone, Copy)]
pub(crate) struct RecMeta {
    pub(crate) hdr_end: u32,
    pub(crate) seq_off: usize,
    pub(crate) seq_len: u32,
    pub(crate) qual_off: usize,
    pub(crate) qual_len: u32,
}

/// One byte-chunk's parse result: a header arena plus per-record metadata.
pub(crate) struct ChunkParse {
    pub(crate) hdr: Vec<u8>,
    pub(crate) recs: Vec<RecMeta>,
}

/// Read one line from `buf[pos..end]`, returning `(content_off, content_len,
/// next_pos)`. Matches `noodles` line semantics: the trailing `\n` (and a `\r`
/// immediately before it) is stripped, but a `\r` at true end-of-input with no
/// following `\n` is kept.
#[inline]
pub(crate) fn take_line(buf: &[u8], pos: usize, end: usize) -> (usize, usize, usize) {
    match memchr::memchr(b'\n', &buf[pos..end]) {
        Some(k) => {
            let line_end = pos + k;
            let mut len = k;
            if len > 0 && buf[line_end - 1] == b'\r' {
                len -= 1;
            }
            (pos, len, line_end + 1)
        }
        None => (pos, end - pos, end),
    }
}

/// Append a record's raw definition line (the bytes after `@`) to `out`
/// **byte-for-byte** — the header is preserved exactly, including a trailing
/// separator or a tab, per fqxv's losslessness contract.
///
/// (Previously this split on the first space/tab and rejoined with a single space
/// to mirror the `noodles` name/description reader, which silently dropped a
/// trailing separator and rewrote a tab to a space — see #49. All read paths now
/// capture the raw definition instead; see [`read_raw_record`].)
pub(crate) fn normalize_header(out: &mut Vec<u8>, def: &[u8]) {
    out.extend_from_slice(def);
}

/// Strip a trailing `\n` (and a `\r` immediately before it) from a read line.
fn strip_eol(line: &[u8]) -> &[u8] {
    let mut end = line.len();
    if end > 0 && line[end - 1] == b'\n' {
        end -= 1;
        if end > 0 && line[end - 1] == b'\r' {
            end -= 1;
        }
    }
    &line[..end]
}

/// Read one 4-line FASTQ record from `r`, capturing the header **byte-exactly**:
/// `def` gets the definition line's bytes after `@` with no separator
/// normalization, `seq`/`qual` get the sequence and quality lines (all three are
/// cleared first). The `+` separator line is read and dropped (fqxv's one
/// sanctioned deviation). Returns `Ok(false)` at a clean EOF before a record.
///
/// This is the shared reader for the streaming / buffered paths that previously
/// used `noodles`, whose name/description split is lossy for the header separator.
pub(crate) fn read_raw_record<R: BufRead>(
    r: &mut R,
    def: &mut Vec<u8>,
    seq: &mut Vec<u8>,
    qual: &mut Vec<u8>,
) -> Result<bool> {
    def.clear();
    seq.clear();
    qual.clear();
    let mut line = Vec::new();
    if r.read_until(b'\n', &mut line)? == 0 {
        return Ok(false); // clean EOF at a record boundary
    }
    let d = strip_eol(&line);
    if d.first() != Some(&b'@') {
        return Err(Error::Malformed("expected '@' at FASTQ record start"));
    }
    def.extend_from_slice(&d[1..]);

    line.clear();
    r.read_until(b'\n', &mut line)?;
    seq.extend_from_slice(strip_eol(&line));

    line.clear();
    r.read_until(b'\n', &mut line)?; // '+' separator line — dropped

    line.clear();
    r.read_until(b'\n', &mut line)?;
    qual.extend_from_slice(strip_eol(&line));
    // Same per-record invariant as the parallel parser (see `parse_chunk`).
    if seq.len() != qual.len() {
        return Err(Error::RecordLengthMismatch {
            seq: seq.len(),
            qual: qual.len(),
        });
    }
    Ok(true)
}

/// True when `o` begins a well-formed 4-line FASTQ record: `@`-line, sequence,
/// `+`-line, and a quality line of the same length as the sequence. The length
/// check makes a false sync (landing on a quality line that happens to start with
/// `@`) astronomically unlikely, so record boundaries can be found in parallel.
pub(crate) fn is_record_start(buf: &[u8], o: usize) -> bool {
    if buf.get(o) != Some(&b'@') {
        return false;
    }
    let n = buf.len();
    let (_, _, p1) = take_line(buf, o, n);
    if p1 >= n {
        return false;
    }
    let (_, seq_len, p2) = take_line(buf, p1, n);
    if buf.get(p2) != Some(&b'+') {
        return false;
    }
    let (_, _, p3) = take_line(buf, p2, n);
    let (_, qual_len, _) = take_line(buf, p3, n);
    seq_len == qual_len
}

/// Find the first record boundary at or after `from` (a line-starting `@` that
/// passes [`is_record_start`]). Used only to split the buffer into parse chunks,
/// so a conservatively-late boundary is harmless.
pub(crate) fn find_record_start(buf: &[u8], from: usize) -> Option<usize> {
    let mut i = from;
    loop {
        let k = memchr::memchr(b'\n', &buf[i..])?;
        let ls = i + k + 1;
        if ls >= buf.len() {
            return None;
        }
        if buf[ls] == b'@' && is_record_start(buf, ls) {
            return Some(ls);
        }
        i = ls;
    }
}

/// Parse the records wholly contained in `buf[start..end)`. `start` and `end`
/// are record boundaries, so every record's four lines lie inside the range
/// (only the final record of the whole input may lack a trailing `\n`).
pub(crate) fn parse_chunk(buf: &[u8], start: usize, end: usize) -> Result<ChunkParse> {
    let mut hdr = Vec::new();
    let mut recs = Vec::new();
    let mut pos = start;
    while pos < end {
        if buf[pos] != b'@' {
            return Err(Error::Malformed("expected FASTQ record start ('@')"));
        }
        let (def_off, def_len, p1) = take_line(buf, pos, end);
        // Header text is everything after the '@', name/description-normalized.
        normalize_header(&mut hdr, &buf[def_off + 1..def_off + def_len]);
        let hdr_end = hdr.len() as u32;

        let (seq_off, seq_len, p2) = take_line(buf, p1, end);
        if buf.get(p2) != Some(&b'+') {
            return Err(Error::Malformed("expected FASTQ '+' separator line"));
        }
        let (_, _, p3) = take_line(buf, p2, end);
        let (qual_off, qual_len, p4) = take_line(buf, p3, end);
        // A valid FASTQ record has equal sequence and quality lengths. Enforce it
        // per-record: the block-level quality-vs-lens check catches only a net
        // mismatch, so two records that mis-compensate (seq>qual then qual>seq)
        // would otherwise pass and silently slide quality across read boundaries.
        if seq_len != qual_len {
            return Err(Error::RecordLengthMismatch {
                seq: seq_len,
                qual: qual_len,
            });
        }

        recs.push(RecMeta {
            hdr_end,
            seq_off,
            seq_len: seq_len as u32,
            qual_off,
            qual_len: qual_len as u32,
        });
        pos = p4;
    }
    Ok(ChunkParse { hdr, recs })
}

/// Parse the whole input buffer into per-chunk record metadata plus the global
/// record index of each chunk's first record (`gstart`) and the total record
/// count `n`.
///
/// The buffer is split into byte-chunks at record boundaries and parsed in
/// parallel. Callers materialize [`RawBlock`]s lazily from this metadata (via
/// [`build_block`]) so only the blocks being compressed right now are resident,
/// rather than a second full copy of the input. Block contents are re-sliced
/// purely by global record index, so the archive is byte-identical regardless of
/// how many chunks/threads did the parsing — determinism holds by construction.
pub(crate) fn parse_chunks(
    buf: &[u8],
    g: usize,
    pool: &rayon::ThreadPool,
) -> Result<(Vec<ChunkParse>, Vec<usize>, usize)> {
    if buf.is_empty() {
        return Ok((Vec::new(), vec![0], 0));
    }

    // Split into ~8 chunks per worker (for load balance), but never finer than
    // ~1 MiB, and resolve each nominal split to a real record boundary.
    let nthreads = pool.current_num_threads().max(1);
    let min_chunk = 1usize << 20;
    let target = nthreads.saturating_mul(8).max(1);
    let n_chunks = (buf.len() / min_chunk).clamp(1, target);
    let mut bounds = Vec::with_capacity(n_chunks + 1);
    bounds.push(0usize);
    for i in 1..n_chunks {
        let nominal = i * buf.len() / n_chunks;
        if let Some(s) = find_record_start(buf, nominal) {
            if *bounds.last().unwrap() < s {
                bounds.push(s);
            }
        }
    }
    bounds.push(buf.len());

    let chunks: Vec<ChunkParse> = pool.install(|| {
        (0..bounds.len() - 1)
            .into_par_iter()
            .map(|i| parse_chunk(buf, bounds[i], bounds[i + 1]))
            .collect::<Result<Vec<_>>>()
    })?;

    // Global record index of each chunk's first record.
    let mut gstart = Vec::with_capacity(chunks.len() + 1);
    let mut acc = 0usize;
    for ch in &chunks {
        gstart.push(acc);
        acc += ch.recs.len();
    }
    gstart.push(acc);
    let n = acc;
    if g > 1 && !n.is_multiple_of(g) {
        return Err(Error::Malformed(
            "interleaved stream ended mid-spot (record count not a multiple of group size)",
        ));
    }
    Ok((chunks, gstart, n))
}

/// Split the parsed reads into row-group ranges `[gs, ge)`, cutting at whichever
/// comes first — `block_reads` reads or `max_bytes` raw sequence bytes — while
/// keeping whole spots together (every cut lands on a multiple of `g` reads).
/// `block_reads` is assumed already rounded to a multiple of `g`.
///
/// Boundaries depend only on the per-read lengths and the two limits, not on how
/// many parse chunks/threads produced them, so the resulting archive is
/// byte-identical regardless of thread count.
pub(crate) fn block_ranges(
    chunks: &[ChunkParse],
    block_reads: usize,
    max_bytes: usize,
    g: usize,
) -> Vec<(usize, usize)> {
    let mut ranges = Vec::new();
    let mut gs = 0usize; // current row group's first global read index
    let mut gi = 0usize; // reads seen so far
    let mut bytes = 0usize; // raw sequence bytes in the current row group
    for chunk in chunks {
        for rec in &chunk.recs {
            gi += 1;
            bytes += rec.seq_len as usize;
            let in_block = gi - gs;
            if in_block.is_multiple_of(g) && (in_block >= block_reads || bytes >= max_bytes) {
                ranges.push((gs, gi));
                gs = gi;
                bytes = 0;
            }
        }
    }
    if gs < gi {
        ranges.push((gs, gi));
    }
    ranges
}
