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
    /// Quality quantization (lossless by default). Ignored when [`Self::no_quality`]
    /// is set — quality is not coded at all in that mode.
    pub quality_binning: QualityBinning,
    /// Discard quality entirely: store names + sequence only, and reconstruct
    /// **FASTA** on decompress. The single largest size lever (quality is the
    /// majority of most archives) but explicitly lossy — the original FASTQ cannot
    /// be recovered. Gated by [`crate::feature::NO_QUALITY`] so older readers refuse
    /// the archive rather than mis-decode it. Not supported with [`Self::reorder`]
    /// (the whole-file reorder layout keeps quality); rejected there.
    pub no_quality: bool,
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
    /// Alignment band half-width for the multi-reference tiler (Nanopore long-read
    /// blocks only — the only path that runs it). Wider bands recover more of each
    /// read's drift against its neighbour at more alignment time. Ratio/speed only;
    /// the block self-describes, so decode is unaffected. Default 256 (the codec's
    /// own default); the CLI raises it at the top effort levels.
    pub tile_band: usize,
    /// Best-of-N reference fan-out for the multi-reference tiler (Nanopore blocks
    /// only). At ONT coverage many earlier reads overlap a span with independent
    /// error patterns; trying `tile_max_refs` of them and keeping the cheapest edit
    /// script is the dominant ONT sequence-ratio lever. Ratio/speed only. Default 1
    /// (greedy single reference); the CLI raises it at the top effort levels, `--max`
    /// to the CoLoRd-parity operating point.
    pub tile_max_refs: usize,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            seq_order: 11,
            seq_hash_order: 0,
            seq_hash_bits: 0,
            block_reads: DEFAULT_BLOCK_READS,
            quality_binning: QualityBinning::Lossless,
            no_quality: false,
            reorder: false,
            keep_order: false,
            rescue: true,
            regenerate_names: false,
            threads: 0,
            platform: None,
            tile_band: 256,
            tile_max_refs: 1,
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
    // Mean sequence length over the peeked records: long reads (nanopore/PacBio)
    // want the buffered shared-reference layout (issue #168), so they can't take
    // the streaming single-end shortcut below.
    let long_read = {
        let (sum, n) = peeked.iter().fold((0u64, 0u64), |(s, n), r| {
            (s + r.sequence().len() as u64, n + 1)
        });
        n > 0 && sum / n > REORDER_MAX_MEAN_LEN
    };
    info!(
        group_size = g,
        reorder = params.reorder,
        long_read,
        "detected layout"
    );

    // Single-end, non-reorder: stream the peeked prefix then the rest through the
    // drive path (#112) instead of buffering the whole input. The other layouts
    // need the whole input — reorder for its global clustering, interleaved for
    // spot regrouping, long-read for its whole-file shared reference (#168) — so
    // they complete the buffer.
    //
    // Explicit Nanopore is the exception among long reads: #211 disabled its
    // whole-file shared reference (the consensus is too noisy to pay off), so the
    // only reason to buffer is gone. Streaming produces a byte-identical archive —
    // `compress_multi` cuts blocks at the same `MAX_BLOCK_SEQ_BYTES` budget as the
    // buffered `block_ranges`, and with no shared reference each Nanopore block
    // codes with the same plain per-block codec either way — while holding one
    // block instead of the whole file (the ~25% peak-heap sink, #225). It must be
    // *explicit* Nanopore: the streaming path detects platform from read names,
    // which is `Unknown` for SRA-style names, so an auto-detected ONT set still
    // takes the buffered path that content-classifies it.
    let stream_ok = !long_read || params.platform == Some(Platform::Nanopore);
    if g == 1 && !params.reorder && stream_ok {
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

/// Mean sequence length over the first few FASTQ records of a buffer, for routing
/// the long-read shared-reference path without a full parse. Mirrors
/// [`mean_read_len`](super::reorder::mean_read_len) over a small sample; returns 0
/// for empty or unparseable input (so it is never treated as long-read).
fn sample_mean_read_len(buf: &[u8]) -> u64 {
    const SAMPLE: u64 = 256;
    let mut fq = noodles_fastq::io::Reader::new(buf);
    let mut rec = noodles_fastq::Record::default();
    let (mut sum, mut n) = (0u64, 0u64);
    while n < SAMPLE {
        match fq.read_record(&mut rec) {
            Ok(0) | Err(_) => break,
            Ok(_) => {
                sum += rec.sequence().len() as u64;
                n += 1;
            }
        }
    }
    sum.checked_div(n).unwrap_or(0)
}

/// Compress an in-memory FASTQ buffer. Long-read, non-reorder input is routed to
/// the shared whole-file reference layout ([`compress_longread_shared_ref`], issue
/// #168); everything else uses the plain per-block layout
/// ([`compress_buffered_plain`]).
pub(crate) fn compress_buffered<W: Write>(
    buf: &[u8],
    writer: W,
    params: Params,
    group_size: u8,
) -> Result<Stats> {
    // Long reads carry cross-read redundancy the within-read model can't see; a
    // single whole-file consensus reference, coded against by every block, captures
    // it once instead of re-storing a reference per block. Short reads keep the
    // plain layout (the reference buys nothing and the streaming path is cheaper).
    if !params.reorder && sample_mean_read_len(buf) > super::reorder::REORDER_MAX_MEAN_LEN {
        // The shared whole-file reference (#168) pays off on low-error long reads
        // (PacBio HiFi): a clean consensus is cheap to store once and cheap to code
        // against. On high-error Nanopore the consensus is noisy, so the reference
        // frame costs more than coding against it saves and the whole-file gate
        // ALWAYS rejects it (#184) — but only after building the whole-file reference
        // AND coding every block against it (`encode_against`), a second full
        // long-read assembly per block that profiling put at ~45% of ONT compress
        // CPU, all thrown away. Skip the shared path on Nanopore and go straight to
        // the plain layout it would fall back to regardless: **byte-identical**
        // output at ~2.3x the speed (measured 7:20 -> 3:11 on a 600 MB E. coli ONT
        // file). Mirrors the #206 LZMA platform gate. (Peak RSS is ~neutral: the
        // plain fallback fans more blocks out concurrently, which offsets the freed
        // second assembly; candidate sequentialization reins the per-block peak in.)
        if resolve_platform_buf(params.platform, buf) != Platform::Nanopore {
            return compress_longread_shared_ref(buf, writer, params, group_size);
        }
    }
    compress_buffered_plain(buf, writer, params, group_size)
}

/// Plain per-block layout: parse the buffer in parallel into blocks, then compress
/// the blocks (in parallel) and write them in order. `group_size` is the
/// interleaving already determined by the caller. This is the layout for short
/// reads and the fallback when the long-read shared reference does not pay off.
pub(crate) fn compress_buffered_plain<W: Write>(
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
    // Byte-budgeted row-group ranges (min of block_reads and a raw-sequence byte
    // cap, on whole-spot boundaries) — a pure function of the read lengths, so
    // determinism holds regardless of thread count.
    let ranges = block_ranges(&chunks, block_reads, MAX_BLOCK_SEQ_BYTES, g);
    write_plain_layout(
        writer, buf, &chunks, &gstart, &ranges, &params, group_size, platform, &pool, None, None,
    )
}

/// Write the plain per-block layout: header, blocks in order, footer.
///
/// `precoded_seq`, when present, supplies each block's already-coded sequence
/// stream (indexed by block) and only names and quality are coded here. The
/// long-read shared-reference path uses it to write the fallback layout without
/// re-running the per-block overlap encode it already did in pass 1.
///
/// `precoded_ns`, when present, supplies each block's already-coded **names AND
/// sequence** streams and only quality is coded here. The single-end reorder
/// never-worse gate uses it to write the plain candidate it decided to keep
/// without re-coding names/sequence. It takes precedence over `precoded_seq`;
/// passing both `None` codes all three streams per block, the ordinary short-read
/// path.
#[allow(clippy::too_many_arguments)]
pub(crate) fn write_plain_layout<W: Write>(
    writer: W,
    buf: &[u8],
    chunks: &[ChunkParse],
    gstart: &[usize],
    ranges: &[(usize, usize)],
    params: &Params,
    group_size: u8,
    platform: Platform,
    pool: &rayon::ThreadPool,
    precoded_seq: Option<&[Vec<u8>]>,
    precoded_ns: Option<&[(Vec<u8>, Vec<u8>)]>,
) -> Result<Stats> {
    let mut w = CrcWriter::new(BufWriter::new(writer));
    write_header(&mut w, params, group_size, platform)?;
    let mut stats = Stats {
        group_size,
        ..Stats::default()
    };
    // Materialize and compress blocks one batch at a time so at most `batch`
    // `RawBlock`s (and their compressed payloads) are ever resident — building
    // every block up front would hold a second full copy of the input alongside
    // `buf`. Each block is a pure function of its global index range, so lazy
    // per-batch building is byte-identical to building them all at once.
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
                    let blk = build_block(buf, chunks, gstart, gs, ge);
                    let payload = match (precoded_ns, precoded_seq) {
                        (Some(ns), _) => {
                            compress_block_with_names_seq(&blk, params, &ns[bi].0, &ns[bi].1)
                        }
                        (None, Some(seqs)) => compress_block_with_seq(&blk, params, &seqs[bi]),
                        (None, None) => compress_block(&blk, params, platform),
                    };
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

/// Whole-file never-worse gate for the long-read shared reference.
///
/// Adopt the reference layout only when the frame plus the reference-coded
/// sequence beats **the plain layout it would otherwise fall back to** —
/// `plain_total` being the sum, per block, of the smaller of the overlap codec and
/// order-k.
///
/// The distinction is the whole point (issue #184). Gating against the order-k
/// total alone is a weaker bar than the fallback actually achieves, so a reference
/// that loses to the per-block overlap codec still clears it. That shipped as a
/// real ONT regression: the frame cost 4.37 MB to save 1.58 MB, inflating the
/// archive by ~2.8 MB while comfortably beating order-k. Ties do not adopt — an
/// equal-size archive plus a reference frame is pure overhead.
///
/// Pulled out as a named function because the arithmetic is what regressed and the
/// end-to-end behaviour cannot be pinned by a small fixture: whether the frame pays
/// for itself depends on how well the whole-file assembly collapses, which only
/// diverges from the per-block assemblies at real coverage and read counts. It is
/// the shared-reference wording of the general [`adopt_over`] rule (#203).
pub(crate) fn adopt_shared_reference(
    ref_frame: usize,
    shared_total: usize,
    plain_total: usize,
) -> bool {
    adopt_over(shared_total, ref_frame, plain_total)
}

/// How many times the *predicted* whole-file saving must exceed the reference
/// frame before the block-0 probe is trusted to stand in for the exact gate.
///
/// The probe extrapolates one block's margin to the whole file, so this is the
/// headroom for blocks that are less favourable than the first. On `ecoli_hifi`
/// the predicted saving is ~5.1x the frame, so the shortcut fires with room to
/// spare; on `ecoli_ont` block 0's plain candidate is *smaller* than its shared
/// one, so the margin is negative and the exact gate runs.
const SHORTCUT_FRAME_MARGIN: u128 = 3;

/// Whether the block-0 probe justifies skipping the plain candidate for the
/// remaining blocks.
///
/// Coding both layouts makes the whole-file gate exact, but it costs a second
/// long-read assembly per block — roughly a 2x compress-time penalty. On input
/// where the reference wins overwhelmingly (HiFi: 0.067 vs 0.102 b/base, a ~30x
/// margin over order-k) that second assembly is pure waste: the plain candidate
/// never had a chance and is coded only to be dropped.
///
/// So extrapolate block 0's per-base margin over the whole file and require the
/// predicted saving to clear the frame by [`SHORTCUT_FRAME_MARGIN`]:
///
/// ```text
/// (plain_0 - shared_0) * total_bases  >  MARGIN * ref_frame * bases_0
/// ```
///
/// Both sides are exact integer arithmetic in `u128`, so the decision is a pure
/// function of the input and stays thread-count invariant. A negative margin
/// (the plain candidate already wins on block 0) never shortcuts.
///
/// This trades the *exact* never-worse guarantee for a predicted one, and only
/// in the regime where the two layouts are nowhere near each other. The downside
/// is bounded: each block's shared candidate is still `min(reference-coded,
/// order-k)`, so a mispredicted shortcut can never push a block below the
/// context model — it can only forgo a per-block overlap win that the probe said
/// was far out of reach.
pub(crate) fn shortcut_to_shared_layout(
    shared_0: usize,
    plain_0: usize,
    bases_0: u64,
    total_bases: u64,
    ref_frame: usize,
) -> bool {
    let Some(margin) = plain_0.checked_sub(shared_0) else {
        return false; // the plain layout already wins block 0
    };
    if bases_0 == 0 || margin == 0 {
        return false;
    }
    let predicted = margin as u128 * u128::from(total_bases);
    let required = SHORTCUT_FRAME_MARGIN * ref_frame as u128 * u128::from(bases_0);
    predicted > required
}

/// Long-read plain layout with a **shared whole-file reference** (issue #168).
///
/// A per-block overlap codec re-assembles and re-stores the same consensus
/// reference in every block; at high coverage a file split into several 256 MiB
/// blocks stores ~one copy of the genome per block. This path assembles the
/// consensus **once** over the whole file, stores it in a single reference frame
/// between the header and the first block, and codes every block's reads against
/// that frozen frame ([`SEQ_METHOD_OVERLAP_REF`]). Placement is per-read against an
/// immutable frame, so a read codes identically regardless of which block holds it
/// — blocks stay 256 MiB, parallel, and independently decodable.
///
/// Two passes over the buffered input:
/// 1. build the reference and code each block's sequence **both ways** — against
///    the shared reference and with the plain per-block codec
///    ([`encode_sequence_stream_shared`]) — holding the small streams;
/// 2. whole-file never-worse gate — adopt the reference layout only when
///    `reference frame + Σ shared sequence` beats `Σ plain sequence`, the layout it
///    would otherwise use; else write the plain layout, reusing the pass-1 plain
///    streams rather than re-coding them.
///
/// Coding both ways is what makes the gate exact (issue #184); the loser's streams
/// are dropped, and neither branch codes the sequence twice.
///
/// On the reference layout, pass 2 codes each block's names and quality and reuses
/// the pass-1 sequence, writing blocks in order with bounded memory.
fn compress_longread_shared_ref<W: Write>(
    buf: &[u8],
    writer: W,
    params: Params,
    group_size: u8,
) -> Result<Stats> {
    let g = group_size.max(1) as usize;
    let block_reads = (params.block_reads.max(1) / g).max(1) * g;
    let pool = build_pool(params.threads)?;
    let (chunks, gstart, _n) = parse_chunks(buf, g, &pool)?;
    let platform = resolve_platform_buf(params.platform, buf);
    let ranges = block_ranges(&chunks, block_reads, MAX_BLOCK_SEQ_BYTES, g);
    let num_blocks = ranges.len();
    let batch = pool.current_num_threads().max(1);

    // Gather the whole file's read lengths and bases (sequence only — the reference
    // needs no quality) in read order, from the parsed record offsets into `buf`.
    let total_bases: usize = chunks
        .iter()
        .flat_map(|c| c.recs.iter())
        .map(|r| r.seq_len as usize)
        .sum();
    let mut all_lens: Vec<u32> = Vec::new();
    let mut all_seq: Vec<u8> = Vec::with_capacity(total_bases);
    for chunk in &chunks {
        for rec in &chunk.recs {
            all_lens.push(rec.seq_len);
            all_seq.extend_from_slice(&buf[rec.seq_off..rec.seq_off + rec.seq_len as usize]);
        }
    }

    // Assemble the shared reference once over every read — a whole-file index, so
    // it takes the whole-file seeding scheme, matching what `encode_against` uses
    // to place reads on it.
    let opts = fqxv_lroverlap::EncodeOpts {
        sketch: sketch_for(platform, SeedContext::WholeFile),
        ..fqxv_lroverlap::EncodeOpts::default()
    };
    let reference = pool.install(|| fqxv_lroverlap::build_reference(&all_lens, &all_seq, &opts))?;
    // Free the whole-file sequence buffer before the block passes allocate.
    drop(all_seq);
    drop(all_lens);

    // No usable reference (no shared locus, e.g. amplicon-free or tiny input): the
    // shared layout can only add overhead, so use the plain per-block layout.
    if reference.is_empty() {
        debug!("shared reference assembled no contigs; using plain layout");
        return compress_buffered_plain(buf, writer, params, group_size);
    }
    let ref_frame = reference.encode()?;

    // Pass 1, probe: code the FIRST block both ways — against the shared reference
    // and with the plain per-block codec. Block 0 is a pure function of the input,
    // so this decision is thread-count invariant.
    let blk0 = build_block(buf, &chunks, &gstart, ranges[0].0, ranges[0].1);
    let bases_0: u64 = blk0.lens.iter().map(|&l| u64::from(l)).sum();
    let (shared_0, plain_0) = pool.install(|| {
        encode_sequence_stream_shared(&blk0.lens, &blk0.seq, &params, platform, &reference)
    })?;
    drop(blk0);

    // When the reference wins the probe by a wide enough margin, skip the plain
    // candidate for the remaining blocks: it costs a second full long-read
    // assembly per block and would only be discarded. See
    // [`shortcut_to_shared_layout`].
    let shortcut = num_blocks > 1
        && shortcut_to_shared_layout(
            shared_0.len(),
            plain_0.len(),
            bases_0,
            total_bases as u64,
            ref_frame.len(),
        );

    let shared_seq: Vec<Vec<u8>> = if shortcut {
        info!(
            shared_0 = shared_0.len(),
            plain_0 = plain_0.len(),
            ref_frame = ref_frame.len(),
            "reference wins the block-0 probe by a wide margin; skipping the plain candidate"
        );
        let mut v = Vec::with_capacity(num_blocks);
        v.push(shared_0);
        v.extend(pool.install(|| {
            ranges[1..]
                .par_iter()
                .map(|&(gs, ge)| -> Result<Vec<u8>> {
                    let blk = build_block(buf, &chunks, &gstart, gs, ge);
                    encode_sequence_stream_shared_only(
                        &blk.lens, &blk.seq, &params, platform, &reference,
                    )
                })
                .collect::<Result<Vec<_>>>()
        })?);
        v
    } else {
        // Exact whole-file never-worse gate (issue #184): the reference frame plus
        // the reference-coded sequence must beat **the plain layout it would
        // otherwise use** — which floors each block at `min(per-block overlap,
        // order-k)`, not at order-k. Gating on order-k alone is a weaker bar than
        // the fallback actually achieves, so a reference that loses to the
        // per-block overlap codec still cleared it: on ONT the frame cost 4.37 MB
        // to save 1.58 MB and was adopted anyway, inflating the archive by ~2.8 MB.
        // Both layouts are coded here, so the comparison is exact rather than
        // predicted; the loser's streams are simply dropped.
        let rest: Vec<(Vec<u8>, Vec<u8>)> = pool.install(|| {
            ranges[1..]
                .par_iter()
                .map(|&(gs, ge)| -> Result<(Vec<u8>, Vec<u8>)> {
                    let blk = build_block(buf, &chunks, &gstart, gs, ge);
                    encode_sequence_stream_shared(
                        &blk.lens, &blk.seq, &params, platform, &reference,
                    )
                })
                .collect::<Result<_>>()
        })?;
        let shared_total = shared_0.len() + rest.iter().map(|(s, _)| s.len()).sum::<usize>();
        let plain_total = plain_0.len() + rest.iter().map(|(_, p)| p.len()).sum::<usize>();

        if !adopt_shared_reference(ref_frame.len(), shared_total, plain_total) {
            info!(
                ref_frame = ref_frame.len(),
                shared_total, plain_total, "shared reference does not pay off; using plain layout"
            );
            // Reuse the pass-1 plain streams: they are byte-identical to a fresh
            // `compress_buffered_plain` pass (same ranges, same per-block choice),
            // so this only skips re-running the expensive per-block overlap encode.
            let mut plain_seq = Vec::with_capacity(num_blocks);
            plain_seq.push(plain_0);
            plain_seq.extend(rest.into_iter().map(|(_, p)| p));
            return write_plain_layout(
                writer,
                buf,
                &chunks,
                &gstart,
                &ranges,
                &params,
                group_size,
                platform,
                &pool,
                Some(&plain_seq),
                None,
            );
        }
        let mut v = Vec::with_capacity(num_blocks);
        v.push(shared_0);
        v.extend(rest.into_iter().map(|(s, _)| s));
        v
    };
    info!(
        contigs = reference.len(),
        ref_bases = reference.total_bases(),
        ref_frame = ref_frame.len(),
        shortcut,
        "shared reference adopted"
    );

    // Pass 2: write header, the reference frame, then blocks (names + quality coded
    // here, reusing the pass-1 sequence) in order.
    let mut w = CrcWriter::new(BufWriter::new(writer));
    let flags = FLAG_PLUS_NORMALIZED | FLAG_GLOBAL_REFERENCE;
    // Carry the sequence-only feature bit alongside the shared-reference bit so a
    // long-read `no_quality` archive is gated the same way the plain layout is.
    let mut required_features = crate::feature::GLOBAL_REFERENCE;
    if params.no_quality {
        required_features |= crate::feature::NO_QUALITY;
    }
    write_header_prefix(
        &mut w,
        params.seq_order,
        binning_tag(params.quality_binning),
        flags,
        group_size,
        platform,
        required_features,
    )?;
    write_framed(&mut w, &ref_frame)?;
    // Framed slice on disk is [4 len][4 crc][bytes]; blocks begin past it.
    let ref_frame_bytes = (4 + CRC_LEN + ref_frame.len()) as u64;

    let mut stats = Stats {
        group_size,
        ..Stats::default()
    };
    let mut index = FooterIndex::new_at(HEADER_LEN as u64 + ref_frame_bytes);
    for batch_start in (0..num_blocks).step_by(batch) {
        let batch_end = (batch_start + batch).min(num_blocks);
        let (blocks, compressed): (Vec<RawBlock>, Vec<Result<Vec<u8>>>) = pool.install(|| {
            (batch_start..batch_end)
                .into_par_iter()
                .map(|bi| {
                    let (gs, ge) = ranges[bi];
                    let blk = build_block(buf, &chunks, &gstart, gs, ge);
                    let payload = compress_block_with_seq(&blk, &params, &shared_seq[bi]);
                    (blk, payload)
                })
                .unzip()
        });
        write_blocks(&mut w, &blocks, compressed, &mut stats, &mut index)?;
    }
    let footer_bytes = write_footer(&mut w, &index, stats.reads)?;
    w.flush()?;
    stats.out_bytes += HEADER_LEN as u64 + ref_frame_bytes + footer_bytes;
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
