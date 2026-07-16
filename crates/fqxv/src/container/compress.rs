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
/// The whole input is read into memory and parsed in parallel (see
/// [`parse_chunks`]) before the blocks are compressed — the serial FASTQ parse
/// was otherwise the dominant single-threaded cost and left most cores idle.
/// Output is byte-identical regardless of thread count.
#[instrument(skip_all, fields(seq_order = params.seq_order, block_reads = params.block_reads, reorder = params.reorder, threads = params.threads))]
pub fn compress<R: Read + Send, W: Write>(
    mut reader: R,
    writer: W,
    params: Params,
) -> Result<Stats> {
    if params.reorder {
        return compress_reordered_whole(reader, writer, params, 1);
    }
    let mut buf = Vec::new();
    reader.read_to_end(&mut buf)?;
    compress_buffered(&buf, writer, params, 1)
}

/// How many leading records [`compress_auto`] reads to decide whether a single
/// stream is interleaved paired data. Four spots' worth is plenty to be
/// confident while staying cheap for the common single-end case.
pub(crate) const AUTODETECT_PEEK: usize = 8;

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
/// paired data from the leading read names (see [`detect_group_size`]). This is
/// what the CLI uses by default for a lone input so `sracha get -Z … | fqxv
/// compress -` archives paired downloads with the right spot grouping and no
/// flag. Detection only ever promotes to paired on unambiguous mate names;
/// otherwise it behaves exactly like [`compress`]. `reorder` mode honours the
/// detected grouping too: paired input is globally clustered and a permutation
/// restores the mate interleaving (see [`encode_reordered`]).
#[instrument(skip_all, fields(seq_order = params.seq_order, block_reads = params.block_reads, reorder = params.reorder, threads = params.threads))]
pub fn compress_auto<R: Read + Send, W: Write>(
    mut reader: R,
    writer: W,
    params: Params,
) -> Result<Stats> {
    let mut buf = Vec::new();
    reader.read_to_end(&mut buf)?;

    // Peek the leading records to decide the layout. `&[u8]` is `BufRead`, so the
    // noodles reader parses straight out of the buffer with no extra copy.
    let mut fq = noodles_fastq::io::Reader::new(&buf[..]);
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
    if params.reorder {
        return encode_reordered(buffer_records(&buf)?, writer, params, g);
    }
    compress_buffered(&buf, writer, params, g)
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
                    let payload = compress_block(&blk, &params);
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
    let pool = build_pool(params.threads)?;
    let batch = pool.current_num_threads().max(1);
    debug!(
        threads = pool.current_num_threads(),
        batch,
        backend = ?fqxv_rans::Backend::detect(),
        "compress pool ready"
    );
    let mut w = CrcWriter::new(BufWriter::new(writer));
    write_header(&mut w, &params, group_size, platform)?;

    let mut stats = Stats {
        group_size,
        ..Stats::default()
    };
    // One batch of buffering: the reader parses the next batch while this thread
    // compresses and writes the current one.
    let mut index = FooterIndex::new();
    let (tx, rx) = std::sync::mpsc::sync_channel::<Result<Vec<RawBlock>>>(1);
    std::thread::scope(|scope| -> Result<()> {
        let reader = scope.spawn(move || {
            loop {
                let mut blocks: Vec<RawBlock> = Vec::with_capacity(batch);
                let mut eof = false;
                for _ in 0..batch {
                    let mut b = RawBlock::default();
                    match fill(&mut b) {
                        Ok(0) => {
                            eof = true;
                            break;
                        }
                        Ok(_) => blocks.push(b),
                        // Surface the parse error to the consumer, then stop.
                        Err(e) => {
                            let _ = tx.send(Err(e));
                            return;
                        }
                    }
                }
                if !blocks.is_empty() && tx.send(Ok(blocks)).is_err() {
                    return; // consumer went away (write/compress error)
                }
                if eof {
                    return;
                }
            }
        });

        // Consume batches; funnel every error through `result` so the receiver
        // is always dropped before the scope joins the reader (a reader blocked
        // on `send` would otherwise deadlock the join).
        let mut result = Ok(());
        for msg in &rx {
            let blocks = match msg {
                Ok(blocks) => blocks,
                Err(e) => {
                    result = Err(e);
                    break;
                }
            };
            debug!(blocks = blocks.len(), "compressing batch");
            let compressed: Vec<Result<Vec<u8>>> = pool.install(|| {
                blocks
                    .par_iter()
                    .map(|b| compress_block(b, &params))
                    .collect()
            });
            if let Err(e) = write_blocks(&mut w, &blocks, compressed, &mut stats, &mut index) {
                result = Err(e);
                break;
            }
        }
        drop(rx);
        reader.join().expect("reader thread panicked");
        result
    })?;

    let footer_bytes = write_footer(&mut w, &index, stats.reads)?;
    w.flush()?;
    stats.out_bytes += HEADER_LEN as u64 + footer_bytes;
    Ok(stats)
}
