//! Public compression entry points and the block-driving machinery.

use super::*;
use rayon::prelude::*;
use tracing::{debug, info, instrument};

/// Default reads per block. Larger blocks populate the sequence model's contexts
/// better (higher ratio) but reduce parallelism and raise memory.
pub(crate) const DEFAULT_BLOCK_READS: usize = 1 << 20;

/// One interleaved spot's records — `(raw_header, sequence, quality)` per member —
/// owned so the platform can be detected before the streaming header is written
/// (see `compress_multi`). The header is the byte-exact definition line.
pub(crate) type PrimedSpot = Vec<(Vec<u8>, Vec<u8>, Vec<u8>)>;

/// Compression parameters.
#[derive(Debug, Clone, Copy)]
pub struct Params {
    /// Sequence context-model order (higher = better ratio, more memory).
    pub seq_order: u8,
    /// Hashed high-order sequence tier order (`0` = disabled). Adds a third escape
    /// tier above `seq_order` that captures deeper context on repetitive data;
    /// gated to top effort levels for its memory cost. Non-reorder path only
    /// (reorder codes just the residual). Decode auto-detects it per block.
    pub seq_hash_order: u8,
    /// Hashed-tier table size in bits: `1 << seq_hash_bits` slots (~`16 << bits`
    /// bytes per active block). Ignored when `seq_hash_order` is 0.
    pub seq_hash_bits: u8,
    /// Reads per block. Blocks are the unit of parallelism and random access;
    /// larger blocks give the order-k sequence model more data to train on.
    pub block_reads: usize,
    /// Quality quantization (lossless by default).
    pub quality_binning: QualityBinning,
    /// Cluster reads (reverse-complement aware) and differentially code the
    /// sequence — captures cross-read duplicate redundancy. Works for grouped
    /// (paired / single-cell) input too; grouped reorder always preserves order.
    pub reorder: bool,
    /// In reorder mode, store a permutation so the original read order is
    /// restored (otherwise reads emerge in clustered order). Forced on for grouped
    /// input (`group_size > 1`), where the permutation reconstructs the spots.
    pub keep_order: bool,
    /// In reorder mode, adaptively use the assembly-aware sequence codecs: each
    /// clustered block is coded with the single-contig (v2), the block-local
    /// literal-rescue (v3, keeps every contig alive and re-attaches would-be
    /// literals via a k-mer-indexed assembly step), and — over a whole-file frozen
    /// global reference (v4, SPRING-style) — the smaller is kept. When the shared
    /// reference nets a whole-file win it is written once and the v4 blocks index
    /// it; otherwise the archive is exactly the v2/v3 layout, so the choice is
    /// never worse than the block-local codecs. Default `true`; set `false` for the
    /// faster v2-only path (which also skips the global assembly). Ignored when
    /// `reorder` is false. Decode dispatches the codec per block from a version
    /// byte, so blocks may mix versions.
    pub rescue: bool,
    /// In single-end reorder mode, if the read names are purely positional (a
    /// counter, e.g. SRA `@RUN.N N`), discard the original order and regenerate
    /// the names from a stored template instead of coding them — no permutation,
    /// no name stream. **Reorder-lossy: reads are renumbered** (sequence/quality
    /// preserved exactly). Ignored unless the names are detected as regenerable;
    /// ignored for grouped input. Opt-in.
    pub regenerate_names: bool,
    /// Worker threads (0 = all available cores); clamped to available cores.
    pub threads: usize,
    /// Sequencing platform to record. `None` (default) auto-detects it from the
    /// leading read names; `Some(_)` forces the recorded value.
    pub platform: Option<Platform>,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            seq_order: 11,
            seq_hash_order: 0,
            seq_hash_bits: 0,
            block_reads: DEFAULT_BLOCK_READS,
            quality_binning: QualityBinning::Lossless,
            reorder: false,
            keep_order: false,
            rescue: true,
            regenerate_names: false,
            threads: 0,
            platform: None,
        }
    }
}

/// Summary of a compress/decompress run.
#[derive(Debug, Default, Clone, Copy)]
pub struct Stats {
    /// Number of reads processed.
    pub reads: u64,
    /// Number of blocks.
    pub blocks: u64,
    /// Bytes written to the output.
    pub out_bytes: u64,
    /// Interleaved members per spot recorded in the archive (1 = single-end,
    /// 2 = paired). Meaningful for compression; 0 from decompression.
    pub group_size: u8,
}

/// Compress single-end FASTQ from `reader` into a `.fqxv` stream.
///
/// Streams the input through the block pipeline rather than buffering it whole
/// (a single-member [`compress_multi`]): the reader parses blocks one at a time
/// and hands them to a pool of compressors, so peak memory is bounded by the
/// blocks in flight, not the file size (#112). Output is byte-identical
/// regardless of thread count.
#[instrument(skip_all, fields(seq_order = params.seq_order, block_reads = params.block_reads, reorder = params.reorder, threads = params.threads))]
pub fn compress<'a, R: Read + Send + 'a, W: Write>(
    reader: R,
    writer: W,
    params: Params,
) -> Result<Stats> {
    if params.reorder {
        return compress_reordered_whole(reader, writer, params, 1);
    }
    compress_multi(
        vec![Box::new(reader) as Box<dyn Read + Send + 'a>],
        writer,
        params,
    )
}

/// How many leading records [`compress_auto`] reads to decide whether a single
/// stream is interleaved paired data. Four spots' worth is plenty to be
/// confident while staying cheap for the common single-end case.
pub(crate) const AUTODETECT_PEEK: usize = 8;

/// Read from `r` into `buf` until it holds `need` complete FASTQ records or the
/// reader is exhausted; returns `true` at EOF. A record is four newline-terminated
/// lines, so this counts to `4 * need` line terminators — unambiguous, unlike the
/// `@`-heuristic boundary finders, which matters because the count routes paired
/// vs single-end detection.
fn read_leading_records<R: Read>(r: &mut R, buf: &mut Vec<u8>, need: usize) -> Result<bool> {
    const STEP: usize = 64 << 10;
    let (mut lines, mut scan) = (0usize, 0usize);
    loop {
        while lines < need * 4 {
            match memchr::memchr(b'\n', &buf[scan..]) {
                Some(k) => {
                    scan += k + 1;
                    lines += 1;
                }
                None => break,
            }
        }
        if lines >= need * 4 {
            return Ok(false);
        }
        let before = buf.len();
        buf.resize(before + STEP, 0);
        let mut got = 0;
        while got < STEP {
            match r.read(&mut buf[before + got..])? {
                0 => {
                    buf.truncate(before + got);
                    return Ok(true);
                }
                k => got += k,
            }
        }
    }
}

/// Split a read name into its mate-independent base and an optional mate marker.
/// Handles the two common conventions: a `/1`|`/2` name suffix, and a mate digit
/// as the first token of the description (`@id 1:N:…` / `@id 2:N:…`).
pub(crate) fn mate_key(rec: &noodles_fastq::Record) -> (&[u8], Option<u8>) {
    let name: &[u8] = rec.name().as_ref();
    if let [base @ .., b'/', m @ (b'1' | b'2')] = name {
        return (base, Some(*m));
    }
    let desc: &[u8] = rec.description().as_ref();
    let marker = desc.first().copied().filter(|c| matches!(c, b'1' | b'2'));
    (name, marker)
}

/// True when `a` and `b` look like the two mates of one spot: same base name and
/// either explicit, differing mate markers (`/1` vs `/2`) or none at all (a bare
/// repeated name, as some interleaved dumps emit).
pub(crate) fn are_mates(a: &noodles_fastq::Record, b: &noodles_fastq::Record) -> bool {
    let (base_a, mate_a) = mate_key(a);
    let (base_b, mate_b) = mate_key(b);
    base_a == base_b
        && match (mate_a, mate_b) {
            (Some(x), Some(y)) => x != y,
            (None, None) => true,
            _ => false,
        }
}

/// Guess the interleaving of a single stream from its leading records. Returns 2
/// only when every peeked pair looks like paired mates; anything ambiguous falls
/// back to 1 (single-end), which is always safe to archive. Single-cell (3-4)
/// interleaving is not auto-detected — pass an explicit group size for that.
pub(crate) fn detect_group_size(peeked: &[noodles_fastq::Record]) -> u8 {
    if peeked.len() < 2 {
        return 1;
    }
    let pairs = peeked.len() / 2;
    for i in 0..pairs {
        if !are_mates(&peeked[2 * i], &peeked[2 * i + 1]) {
            return 1;
        }
    }
    2
}

/// Compress a single FASTQ stream, auto-detecting whether it is interleaved
/// paired data from the leading read names (see `detect_group_size`). This is
/// what the CLI uses by default for a lone input so `sracha get -Z … | fqxv
/// compress -` archives paired downloads with the right spot grouping and no
/// flag. Detection only ever promotes to paired on unambiguous mate names;
/// otherwise it behaves exactly like [`compress`]. `reorder` mode honours the
/// detected grouping too: paired input is globally clustered and a permutation
/// restores the mate interleaving (see `encode_reordered`).
#[instrument(skip_all, fields(seq_order = params.seq_order, block_reads = params.block_reads, reorder = params.reorder, threads = params.threads))]
pub fn compress_auto<'a, R: Read + Send + 'a, W: Write>(
    mut reader: R,
    writer: W,
    params: Params,
) -> Result<Stats> {
    // Read only enough to peek the leading records for layout detection. FASTQ is
    // strictly four lines per record, so `AUTODETECT_PEEK` records is exactly
    // `4 * AUTODETECT_PEEK` complete lines — no whole-file buffer just to decide.
    let mut prefix = Vec::new();
    let eof = read_leading_records(&mut reader, &mut prefix, AUTODETECT_PEEK)?;

    let mut fq = noodles_fastq::io::Reader::new(&prefix[..]);
    let mut peeked: Vec<noodles_fastq::Record> = Vec::with_capacity(AUTODETECT_PEEK);
    for _ in 0..AUTODETECT_PEEK {
        let mut rec = noodles_fastq::Record::default();
        if fq.read_record(&mut rec)? == 0 {
            break;
        }
        peeked.push(rec);
    }
    let g = detect_group_size(&peeked);
    info!(group_size = g, reorder = params.reorder, "detected layout");

    // Single-end, non-reorder: stream the peeked prefix then the rest through the
    // drive path (#112). The other layouts need the whole input — reorder for its
    // global clustering, interleaved for spot regrouping — so complete the buffer.
    if g == 1 && !params.reorder {
        return compress(prefix.as_slice().chain(reader), writer, params);
    }
    if !eof {
        reader.read_to_end(&mut prefix)?;
    }
    if params.reorder {
        return encode_reordered(buffer_records(&prefix)?, writer, params, g);
    }
    compress_buffered(&prefix, writer, params, g)
}

/// Compress a single FASTQ stream whose records are *already* interleaved per
/// spot (`m0₀, m1₀, …, m0₁, m1₁, …`) — e.g. the interleaved paired output of
/// `sracha get -Z`. Equivalent to [`compress_multi`] but from one reader, so a
/// download can be archived in one pass with nothing hitting disk.
///
/// `group_size` is the number of interleaved members per spot (2 = paired). The
/// stream's total record count must be a multiple of `group_size`; a trailing
/// partial spot is an error. Restore with [`decompress_split`], or stream
/// interleaved with [`decompress`].
#[instrument(skip_all, fields(seq_order = params.seq_order, block_reads = params.block_reads, group_size, threads = params.threads))]
pub fn compress_interleaved<R: Read + Send, W: Write>(
    mut reader: R,
    writer: W,
    params: Params,
    group_size: u8,
) -> Result<Stats> {
    let g = group_size.max(1);
    if g == 1 {
        return compress(reader, writer, params);
    }
    let mut buf = Vec::new();
    reader.read_to_end(&mut buf)?;
    if params.reorder {
        let all = buffer_records(&buf)?;
        if all.n_reads() % g as usize != 0 {
            return Err(Error::Malformed(
                "interleaved stream ended mid-spot (record count not a multiple of group size)",
            ));
        }
        return encode_reordered(all, writer, params, g);
    }
    compress_buffered(&buf, writer, params, g)
}

/// Compress `G >= 1` per-spot read files (paired mates, single-cell R1/R2/I1/I2,
/// …) into one `.fqxv` stream, interleaving them.
///
/// Readers are consumed in lockstep; unequal read counts are an error. Restore
/// with [`decompress_split`], or stream interleaved with [`decompress`].
#[instrument(skip_all, fields(seq_order = params.seq_order, block_reads = params.block_reads, inputs = readers.len(), threads = params.threads))]
pub fn compress_multi<'a, W: Write>(
    readers: Vec<Box<dyn Read + Send + 'a>>,
    writer: W,
    params: Params,
) -> Result<Stats> {
    let g = readers.len();
    if g == 0 {
        return Err(Error::Malformed("no input readers"));
    }
    if g > u8::MAX as usize {
        return Err(Error::Malformed("too many interleaved inputs"));
    }
    let mut fqs: Vec<BufReader<Box<dyn Read + Send + 'a>>> =
        readers.into_iter().map(BufReader::new).collect();
    // Per-member scratch buffers, reused across records (raw definition line,
    // sequence, quality). `read_raw_record` keeps the header byte-exact.
    let mut defs: Vec<Vec<u8>> = vec![Vec::new(); g];
    let mut seqs: Vec<Vec<u8>> = vec![Vec::new(); g];
    let mut quals: Vec<Vec<u8>> = vec![Vec::new(); g];

    if params.reorder {
        // Buffer every spot in interleaved order (m0₀, m1₀, …), then globally
        // cluster; the stored permutation restores this spot order on decode, so
        // grouped reorder is order-preserving and de-interleaves cleanly.
        let mut all = RawBlock::default();
        loop {
            for j in 0..g {
                if !read_raw_record(&mut fqs[j], &mut defs[j], &mut seqs[j], &mut quals[j])? {
                    if j == 0 {
                        return encode_reordered(all, writer, params, g as u8);
                    }
                    return Err(Error::Malformed("inputs have unequal read counts"));
                }
            }
            for j in 0..g {
                all.push_raw(&defs[j], &seqs[j], &quals[j]);
            }
        }
    }

    // Keep whole spots together: round the block target down to a multiple of g.
    let block_reads = (params.block_reads / g).max(1) * g;
    // Prime the first spot so the platform can be detected before the header is
    // written (this path streams, so no full buffer exists to peek). The fill
    // closure emits the primed spot before reading any further records, so block
    // boundaries are byte-identical to reading it inline.
    let mut primed: PrimedSpot = Vec::with_capacity(g);
    for j in 0..g {
        if !read_raw_record(&mut fqs[j], &mut defs[j], &mut seqs[j], &mut quals[j])? {
            if j == 0 {
                break; // empty input
            }
            return Err(Error::Malformed("inputs have unequal read counts"));
        }
        primed.push((defs[j].clone(), seqs[j].clone(), quals[j].clone()));
    }
    let refs: Vec<&[u8]> = primed.iter().map(|(d, _, _)| d.as_slice()).collect();
    let platform = params.platform.unwrap_or_else(|| detect_platform(&refs));
    let mut primed = Some(primed).filter(|p| !p.is_empty());
    drive(writer, params, g as u8, platform, |b| {
        // Emit the primed first spot into the first block before reading on.
        if let Some(spot) = primed.take() {
            for (header, seq, qual) in &spot {
                b.push_raw(header, seq, qual);
            }
        }
        // Cut on reads OR the raw-sequence byte budget, whichever comes first;
        // the loop reads whole spots, so a byte cut still lands on a spot
        // boundary. Matches the byte budgeting in `block_ranges`.
        while b.n_reads() < block_reads && b.seq.len() < MAX_BLOCK_SEQ_BYTES {
            // Read one record from each input; member 0 EOF ends cleanly.
            for j in 0..g {
                if !read_raw_record(&mut fqs[j], &mut defs[j], &mut seqs[j], &mut quals[j])? {
                    if j == 0 {
                        return Ok(b.n_reads());
                    }
                    return Err(Error::Malformed("inputs have unequal read counts"));
                }
            }
            for j in 0..g {
                b.push_raw(&defs[j], &seqs[j], &quals[j]);
            }
        }
        Ok(b.n_reads())
    })
}

/// Compress an in-memory FASTQ buffer: parse it in parallel into blocks, then
/// compress the blocks (in parallel) and write them in order. `group_size` is
/// the interleaving already determined by the caller.
pub(crate) fn compress_buffered<W: Write>(
    buf: &[u8],
    writer: W,
    params: Params,
    group_size: u8,
) -> Result<Stats> {
    let g = group_size.max(1) as usize;
    // Keep whole spots together: round the block target down to a multiple of g.
    let block_reads = (params.block_reads.max(1) / g).max(1) * g;
    let pool = build_pool(params.threads)?;
    debug!(
        threads = pool.current_num_threads(),
        block_reads,
        group_size,
        backend = ?fqxv_rans::Backend::detect(),
        "compress pool ready"
    );
    let (chunks, gstart, _n) = parse_chunks(buf, g, &pool)?;
    let platform = resolve_platform_buf(params.platform, buf);

    let mut w = CrcWriter::new(BufWriter::new(writer));
    write_header(&mut w, &params, group_size, platform)?;
    let mut stats = Stats {
        group_size,
        ..Stats::default()
    };
    // Materialize and compress blocks one batch at a time so at most `batch`
    // `RawBlock`s (and their compressed payloads) are ever resident — building
    // every block up front would hold a second full copy of the input alongside
    // `buf`. Each block is a pure function of its global index range, so lazy
    // per-batch building is byte-identical to building them all at once.
    // Byte-budgeted row-group ranges (min of block_reads and a raw-sequence byte
    // cap, on whole-spot boundaries) — a pure function of the read lengths, so
    // determinism holds regardless of thread count.
    let ranges = block_ranges(&chunks, block_reads, MAX_BLOCK_SEQ_BYTES, g);
    let num_blocks = ranges.len();
    let batch = pool.current_num_threads().max(1);
    let mut index = FooterIndex::new();
    for batch_start in (0..num_blocks).step_by(batch) {
        let batch_end = (batch_start + batch).min(num_blocks);
        let (blocks, compressed): (Vec<RawBlock>, Vec<Result<Vec<u8>>>) = pool.install(|| {
            (batch_start..batch_end)
                .into_par_iter()
                .map(|bi| {
                    let (gs, ge) = ranges[bi];
                    let blk = build_block(buf, &chunks, &gstart, gs, ge);
                    let payload = compress_block(&blk, &params, platform);
                    (blk, payload)
                })
                .unzip()
        });
        write_blocks(&mut w, &blocks, compressed, &mut stats, &mut index)?;
    }
    let footer_bytes = write_footer(&mut w, &index, stats.reads)?;
    w.flush()?;
    stats.out_bytes += HEADER_LEN as u64 + footer_bytes;
    Ok(stats)
}

/// Shared block driver: `fill` populates one [`RawBlock`] and returns the number
/// of reads it added (0 at EOF). Blocks are compressed in parallel, written in
/// order.
///
/// Parsing input (the `fill` calls, single-threaded because the FASTQ stream is
/// sequential) runs on a dedicated thread and stays a batch ahead via a bounded
/// channel, so it overlaps the parallel compression of the previous batch
/// instead of alternating with it — the parse phase was otherwise a serial
/// stretch that left cores idle and capped utilization.
pub(crate) fn drive<W, F>(
    writer: W,
    params: Params,
    group_size: u8,
    platform: Platform,
    mut fill: F,
) -> Result<Stats>
where
    W: Write,
    F: FnMut(&mut RawBlock) -> Result<usize> + Send,
{
    let nworkers = resolve_threads(params.threads);
    debug!(
        threads = nworkers,
        backend = ?fqxv_rans::Backend::detect(),
        "compress pipeline ready"
    );
    let mut w = CrcWriter::new(BufWriter::new(writer));
    write_header(&mut w, &params, group_size, platform)?;

    let mut stats = Stats {
        group_size,
        ..Stats::default()
    };
    let mut index = FooterIndex::new();

    // A true block-level pipeline, not a batch barrier. The reader parses blocks
    // one at a time (the FASTQ stream is sequential) and streams them to a pool of
    // `nworkers` compressors; a writer drains completed blocks in index order. So
    // the reader parses block N+k while the workers compress N..N+k-1 and the
    // writer emits everything before N — the serial parse overlaps the parallel
    // compression instead of preceding it. The old code accumulated a whole batch
    // before compressing any of it, which on a file of only a few blocks meant a
    // serial parse phase with every core idle, then a compress phase.
    //
    // Determinism is unaffected: blocks are compressed independently and the
    // writer reorders by index, so the output is byte-identical regardless of how
    // the workers interleave. Bound the channels by `nworkers` so at most ~2x that
    // many blocks are ever resident.
    let (work_tx, work_rx) = std::sync::mpsc::sync_channel::<(usize, RawBlock)>(nworkers);
    let work_rx = std::sync::Arc::new(std::sync::Mutex::new(work_rx));
    #[allow(clippy::type_complexity)]
    let (done_tx, done_rx) =
        std::sync::mpsc::sync_channel::<(usize, RawBlock, Result<Vec<u8>>)>(nworkers);

    std::thread::scope(|scope| -> Result<()> {
        // Reader: serial parse, streaming each block as it is built.
        let reader = scope.spawn(move || -> Result<()> {
            let mut idx = 0usize;
            loop {
                let mut b = RawBlock::default();
                match fill(&mut b)? {
                    0 => break,
                    _ => {
                        if work_tx.send((idx, b)).is_err() {
                            break; // workers gone (downstream error)
                        }
                        idx += 1;
                    }
                }
            }
            Ok(()) // dropping work_tx signals EOF to the workers
        });

        // Workers: pull the next block under a short lock, compress, hand on.
        let params_ref = &params;
        for _ in 0..nworkers {
            let work_rx = std::sync::Arc::clone(&work_rx);
            let done_tx = done_tx.clone();
            scope.spawn(move || {
                loop {
                    let item = work_rx.lock().expect("work lock").recv();
                    match item {
                        Ok((idx, blk)) => {
                            let payload = compress_block(&blk, params_ref, platform);
                            if done_tx.send((idx, blk, payload)).is_err() {
                                break; // writer gone
                            }
                        }
                        Err(_) => break, // reader done and channel drained
                    }
                }
            });
        }
        drop(done_tx); // so `done_rx` ends once every worker has finished

        // Writer: reorder completed blocks by index and emit contiguously.
        let mut result = Ok(());
        let mut pending: std::collections::HashMap<usize, (RawBlock, Result<Vec<u8>>)> =
            std::collections::HashMap::new();
        let mut next = 0usize;
        for (idx, blk, payload) in &done_rx {
            pending.insert(idx, (blk, payload));
            while let Some((blk, payload)) = pending.remove(&next) {
                if let Err(e) = write_blocks(&mut w, &[blk], vec![payload], &mut stats, &mut index)
                {
                    result = Err(e);
                    break;
                }
                next += 1;
            }
            if result.is_err() {
                break;
            }
        }
        // Drain any remaining messages so a worker blocked on `send` cannot
        // deadlock the scope join after a write error.
        drop(done_rx);
        // Reader errors (parse failures) take priority over a downstream error,
        // matching the buffered path, which sees the parse error first.
        match reader.join().expect("reader thread panicked") {
            Err(e) => Err(e),
            Ok(()) => result,
        }
    })?;

    let footer_bytes = write_footer(&mut w, &index, stats.reads)?;
    w.flush()?;
    stats.out_bytes += HEADER_LEN as u64 + footer_bytes;
    Ok(stats)
}
