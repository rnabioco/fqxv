//! The `.fqxv` container: a header followed by independent, parallel-codable
//! blocks.
//!
//! ```text
//! [4] magic "FQXV"
//! [2] format version (LE)
//! [1] sequence context order (k)
//! [1] quality binning tag
//! [1] flags (bit0: '+' normalized; bit1: reordered; bit2: order preserved)
//! [1] group size G (reads interleaved per spot: 1 single-end, 2 paired,
//!                   3-4 single-cell R1/R2/I1[/I2], ...)
//! repeated until EOF:
//!   [8] block payload length (LE)
//!   [ ] block payload
//! block payload (plain):
//!   [4] n_reads (LE)
//!   [4] names_len (LE)  [ ] names   (fqxv-tokenizer)
//!   [4] seq_len   (LE)  [ ] seq     (fqxv-seq)
//!   [4] qual_len  (LE)  [ ] qual    (fqxv-fqzcomp)
//! block payload (reordered): as plain, but after n_reads:
//!   [ceil(n/8)] flip bitmap        (reads stored reverse-complemented)
//!   [4] perm_len (LE) [ ] perm     (byte-plane-split u32 permutation, rANS'd;
//!                                   empty unless order is kept)
//!   seq is always fqxv-reorder clustered/differential. With order kept, names
//!   and qual are coded in ORIGINAL order (the permutation reunites each
//!   clustered sequence with its read, so their models aren't scrambled);
//!   without it, names and qual follow the clustered order.
//!
//! `--reorder --keep-order` instead uses a whole-file, globally-clustered layout
//! (flag bit3), SPRING-style: all reads are clustered in one pass, then the
//! clustered sequence and the original-order names/quality are each coded in
//! independent moderate blocks that fan out across cores. Clustering is global,
//! so block size is free to be moderate for parallelism without hurting ratio.
//! Layout after the header: `[8] n  [ ] flip  [ ] perm  [4] n_blocks
//! [seq block]*n  [ [names][qual] ]*n` (each `[ ]` is a `[u32 len][bytes]` frame;
//! seq blocks are clustered order, name/qual blocks original order).
//! ```
//!
//! When `G > 1`, reads are interleaved per spot (`m0₀, m1₀, …, m0₁, m1₁, …`).
//! Blocks always hold whole spots and start on member 0, so a block splits back
//! into the `G` files by local read index mod `G`. Interleaving lets the name
//! tokenizer collapse the near-identical mate names and keeps reads from one
//! spot adjacent for the sequence model. [`decompress`] streams interleaved
//! FASTQ (pipe to an aligner); [`decompress_split`] restores the `G` files.

use std::io::{self, BufReader, BufWriter, Read, Write};

use rayon::prelude::*;
use tracing::{debug, info, instrument, trace};

use crate::{Error, Result, FORMAT_VERSION, MAGIC};
use fqxv_fqzcomp::QualityBinning;

/// Default reads per block. Larger blocks populate the sequence model's contexts
/// better (higher ratio) but reduce parallelism and raise memory.
const DEFAULT_BLOCK_READS: usize = 1 << 20;
const HEADER_LEN: usize = 10;
const FLAG_PLUS_NORMALIZED: u8 = 0x01;
const FLAG_REORDERED: u8 = 0x02;
const FLAG_KEEP_ORDER: u8 = 0x04;
/// Whole-file, globally-clustered reorder layout (see `compress_reordered_whole`)
/// as opposed to the older per-block reorder blocks.
const FLAG_GLOBAL_REORDER: u8 = 0x08;
/// Minimizer length for clustering in reorder mode.
const REORDER_K: usize = 15;
/// Reads per block in whole-file (global-cluster) reorder mode. Moderate, so the
/// sequence and name/quality blocks fan out across cores. Clustering is global,
/// so block size no longer trades against ratio — only parallelism and the
/// per-block model reset (cheap, since clustered duplicates collapse to MATCH).
const REORDER_BLOCK_READS: usize = 1 << 18;

/// Compression parameters.
#[derive(Debug, Clone, Copy)]
pub struct Params {
    /// Sequence context-model order (higher = better ratio, more memory).
    pub seq_order: u8,
    /// Reads per block. Blocks are the unit of parallelism and random access;
    /// larger blocks give the order-k sequence model more data to train on.
    pub block_reads: usize,
    /// Quality quantization (lossless by default).
    pub quality_binning: QualityBinning,
    /// Cluster reads (reverse-complement aware) and differentially code the
    /// sequence — captures cross-read duplicate redundancy. Single-end only.
    pub reorder: bool,
    /// In reorder mode, store a permutation so the original read order is
    /// restored (otherwise reads emerge in clustered order).
    pub keep_order: bool,
    /// Worker threads (0 = all available cores); clamped to available cores.
    pub threads: usize,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            seq_order: 11,
            block_reads: DEFAULT_BLOCK_READS,
            quality_binning: QualityBinning::Lossless,
            reorder: false,
            keep_order: false,
            threads: 0,
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

/// Container header + per-stream size summary, from [`inspect`] / [`peek`].
#[derive(Debug, Default, Clone)]
pub struct Info {
    /// Sequence context order.
    pub seq_order: u8,
    /// Quality binning tag (0 = lossless).
    pub quality_binning: u8,
    /// Whether the `+` line was normalized.
    pub plus_normalized: bool,
    /// Reads interleaved per spot (1 = single-end, 2 = paired, 3-4 = single-cell).
    pub group_size: u8,
    /// Whether reads were clustered/reordered.
    pub reordered: bool,
    /// Number of blocks (0 from [`peek`]).
    pub blocks: u64,
    /// Total reads (0 from [`peek`]).
    pub reads: u64,
    /// Compressed bytes in the name stream.
    pub names_bytes: u64,
    /// Compressed sequence bytes.
    pub seq_bytes: u64,
    /// Compressed quality bytes.
    pub qual_bytes: u64,
}

/// One block of parsed FASTQ records. Header text is packed into a single arena
/// (`header_buf` + cumulative `header_ends`) rather than a `Vec` per record — the
/// parse loop is single-threaded and feeds the parallel compressors, so avoiding
/// a per-read allocation keeps that feed from starving the pool.
#[derive(Default)]
struct RawBlock {
    header_buf: Vec<u8>,
    header_ends: Vec<u32>,
    lens: Vec<u32>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl RawBlock {
    fn push(&mut self, name: &[u8], description: &[u8], seq: &[u8], qual: &[u8]) {
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
    fn push_raw(&mut self, header: &[u8], seq: &[u8], qual: &[u8]) {
        self.header_buf.extend_from_slice(header);
        self.header_ends.push(self.header_buf.len() as u32);
        self.lens.push(seq.len() as u32);
        self.seq.extend_from_slice(seq);
        self.qual.extend_from_slice(qual);
    }

    /// Number of records in the block.
    fn n_reads(&self) -> usize {
        self.header_ends.len()
    }

    /// The `i`th record's header bytes.
    fn header(&self, i: usize) -> &[u8] {
        let start = if i == 0 {
            0
        } else {
            self.header_ends[i - 1] as usize
        };
        &self.header_buf[start..self.header_ends[i] as usize]
    }

    /// Borrowed slices for every header, in record order.
    fn header_refs(&self) -> Vec<&[u8]> {
        let mut refs = Vec::with_capacity(self.header_ends.len());
        let mut start = 0usize;
        for &end in &self.header_ends {
            refs.push(&self.header_buf[start..end as usize]);
            start = end as usize;
        }
        refs
    }
}

/// Compress single-end FASTQ from `reader` into a `.fqxv` stream.
///
/// The whole input is read into memory and parsed in parallel (see
/// [`parse_blocks`]) before the blocks are compressed — the serial FASTQ parse
/// was otherwise the dominant single-threaded cost and left most cores idle.
/// Output is byte-identical regardless of thread count.
#[instrument(skip_all, fields(seq_order = params.seq_order, block_reads = params.block_reads, reorder = params.reorder, threads = params.threads))]
pub fn compress<R: Read + Send, W: Write>(
    mut reader: R,
    writer: W,
    params: Params,
) -> Result<Stats> {
    if params.reorder && params.keep_order {
        return compress_reordered_whole(reader, writer, params);
    }
    let mut buf = Vec::new();
    reader.read_to_end(&mut buf)?;
    compress_buffered(&buf, writer, params, 1)
}

/// How many leading records [`compress_auto`] reads to decide whether a single
/// stream is interleaved paired data. Four spots' worth is plenty to be
/// confident while staying cheap for the common single-end case.
const AUTODETECT_PEEK: usize = 8;

/// Split a read name into its mate-independent base and an optional mate marker.
/// Handles the two common conventions: a `/1`|`/2` name suffix, and a mate digit
/// as the first token of the description (`@id 1:N:…` / `@id 2:N:…`).
fn mate_key(rec: &noodles_fastq::Record) -> (&[u8], Option<u8>) {
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
fn are_mates(a: &noodles_fastq::Record, b: &noodles_fastq::Record) -> bool {
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
fn detect_group_size(peeked: &[noodles_fastq::Record]) -> u8 {
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
/// otherwise it behaves exactly like [`compress`]. In `reorder` mode the stream
/// is always treated as single-end (reorder is single-end only).
#[instrument(skip_all, fields(seq_order = params.seq_order, block_reads = params.block_reads, reorder = params.reorder, threads = params.threads))]
pub fn compress_auto<R: Read + Send, W: Write>(
    mut reader: R,
    writer: W,
    params: Params,
) -> Result<Stats> {
    if params.reorder && params.keep_order {
        // Global-cluster reorder buffers all reads itself; skip peek/autodetect
        // (reorder is single-end, so there is no group size to detect).
        return compress_reordered_whole(reader, writer, params);
    }
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
    let g = if params.reorder {
        1
    } else {
        detect_group_size(&peeked)
    };
    info!(group_size = g, "detected layout");
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
    if params.reorder {
        return Err(Error::Malformed(
            "reorder is single-end only (would break spot grouping)",
        ));
    }
    let mut buf = Vec::new();
    reader.read_to_end(&mut buf)?;
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
    if params.reorder && g > 1 {
        return Err(Error::Malformed(
            "reorder is single-end only (would break spot grouping)",
        ));
    }
    let mut fqs: Vec<_> = readers
        .into_iter()
        .map(|r| noodles_fastq::io::Reader::new(BufReader::new(r)))
        .collect();
    let mut recs: Vec<_> = (0..g).map(|_| noodles_fastq::Record::default()).collect();

    // Keep whole spots together: round the block target down to a multiple of g.
    let block_reads = (params.block_reads / g).max(1) * g;
    drive(writer, params, g as u8, |b| {
        while b.n_reads() < block_reads {
            // Read one record from each input; member 0 EOF ends cleanly.
            let mut got = 0;
            for j in 0..g {
                if fqs[j].read_record(&mut recs[j])? == 0 {
                    if j == 0 {
                        return Ok(b.n_reads());
                    }
                    return Err(Error::Malformed("inputs have unequal read counts"));
                }
                got += 1;
            }
            debug_assert_eq!(got, g);
            for r in &recs {
                b.push(r.name(), r.description(), r.sequence(), r.quality_scores());
            }
        }
        Ok(b.n_reads())
    })
}

// --- parallel FASTQ parsing --------------------------------------------------

/// One record's location: its normalized header lives in the owning chunk's
/// header arena (`[prev.hdr_end .. hdr_end)`); its sequence and quality bytes are
/// contiguous ranges of the input buffer (CR/LF already excluded).
#[derive(Clone, Copy)]
struct RecMeta {
    hdr_end: u32,
    seq_off: usize,
    seq_len: u32,
    qual_off: usize,
    qual_len: u32,
}

/// One byte-chunk's parse result: a header arena plus per-record metadata.
struct ChunkParse {
    hdr: Vec<u8>,
    recs: Vec<RecMeta>,
}

/// Read one line from `buf[pos..end]`, returning `(content_off, content_len,
/// next_pos)`. Matches `noodles` line semantics: the trailing `\n` (and a `\r`
/// immediately before it) is stripped, but a `\r` at true end-of-input with no
/// following `\n` is kept.
#[inline]
fn take_line(buf: &[u8], pos: usize, end: usize) -> (usize, usize, usize) {
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

/// Build a record's header exactly as [`RawBlock::push`] would from the
/// `noodles` name/description split: name is the bytes before the first space or
/// tab; the separator becomes a single space; an empty description is dropped.
fn normalize_header(out: &mut Vec<u8>, def: &[u8]) {
    match def.iter().position(|&b| b == b' ' || b == b'\t') {
        None => out.extend_from_slice(def),
        Some(i) => {
            out.extend_from_slice(&def[..i]);
            let desc = &def[i + 1..];
            if !desc.is_empty() {
                out.push(b' ');
                out.extend_from_slice(desc);
            }
        }
    }
}

/// True when `o` begins a well-formed 4-line FASTQ record: `@`-line, sequence,
/// `+`-line, and a quality line of the same length as the sequence. The length
/// check makes a false sync (landing on a quality line that happens to start with
/// `@`) astronomically unlikely, so record boundaries can be found in parallel.
fn is_record_start(buf: &[u8], o: usize) -> bool {
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
fn find_record_start(buf: &[u8], from: usize) -> Option<usize> {
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
fn parse_chunk(buf: &[u8], start: usize, end: usize) -> Result<ChunkParse> {
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

/// Assemble one output block from the globally-ordered record range `[gs, ge)`,
/// copying each record's header (from its chunk's arena) and sequence/quality
/// (from `buf`) into a fresh [`RawBlock`]. `gstart[c]` is the global index of
/// chunk `c`'s first record.
fn build_block(
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

/// Parse the whole input buffer into ordered [`RawBlock`]s of up to `block_reads`
/// records each (a multiple of `g`).
///
/// The buffer is split into byte-chunks at record boundaries and parsed in
/// parallel, then re-sliced into blocks purely by global record index — so block
/// contents (and thus the archive) are byte-identical regardless of how many
/// chunks/threads did the parsing. Determinism holds by construction.
fn parse_blocks(
    buf: &[u8],
    g: usize,
    block_reads: usize,
    pool: &rayon::ThreadPool,
) -> Result<Vec<RawBlock>> {
    if buf.is_empty() {
        return Ok(Vec::new());
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
    if n == 0 {
        return Ok(Vec::new());
    }

    let num_blocks = n.div_ceil(block_reads);
    let blocks: Vec<RawBlock> = pool.install(|| {
        (0..num_blocks)
            .into_par_iter()
            .map(|b| {
                let gs = b * block_reads;
                let ge = ((b + 1) * block_reads).min(n);
                build_block(buf, &chunks, &gstart, gs, ge)
            })
            .collect()
    });
    Ok(blocks)
}

/// Write the container header.
fn write_header<W: Write>(w: &mut W, params: &Params, group_size: u8) -> Result<()> {
    let mut flags = FLAG_PLUS_NORMALIZED;
    if params.reorder {
        flags |= FLAG_REORDERED;
        if params.keep_order {
            flags |= FLAG_KEEP_ORDER;
        }
    }
    w.write_all(&MAGIC)?;
    w.write_all(&FORMAT_VERSION.to_le_bytes())?;
    w.write_all(&[
        params.seq_order,
        binning_tag(params.quality_binning),
        flags,
        group_size,
    ])?;
    Ok(())
}

/// Compress an in-memory FASTQ buffer: parse it in parallel into blocks, then
/// compress the blocks (in parallel) and write them in order. `group_size` is
/// the interleaving already determined by the caller.
fn compress_buffered<W: Write>(
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
    let blocks = parse_blocks(buf, g, block_reads, &pool)?;

    let mut w = BufWriter::new(writer);
    write_header(&mut w, &params, group_size)?;
    let mut stats = Stats {
        group_size,
        ..Stats::default()
    };
    // Compress in batches so at most `batch` compressed payloads are buffered
    // before being written in block order.
    let batch = pool.current_num_threads().max(1);
    for chunk in blocks.chunks(batch) {
        let compressed: Vec<Result<Vec<u8>>> = pool.install(|| {
            chunk
                .par_iter()
                .map(|b| compress_block(b, &params))
                .collect()
        });
        write_blocks(&mut w, chunk, compressed, &mut stats)?;
    }
    w.flush()?;
    stats.out_bytes += HEADER_LEN as u64;
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
fn drive<W, F>(writer: W, params: Params, group_size: u8, mut fill: F) -> Result<Stats>
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
    let mut w = BufWriter::new(writer);
    write_header(&mut w, &params, group_size)?;

    let mut stats = Stats {
        group_size,
        ..Stats::default()
    };
    // One batch of buffering: the reader parses the next batch while this thread
    // compresses and writes the current one.
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
            if let Err(e) = write_blocks(&mut w, &blocks, compressed, &mut stats) {
                result = Err(e);
                break;
            }
        }
        drop(rx);
        reader.join().expect("reader thread panicked");
        result
    })?;

    w.flush()?;
    stats.out_bytes += HEADER_LEN as u64;
    Ok(stats)
}

/// Write a batch's compressed payloads in order, updating `stats`.
fn write_blocks<W: Write>(
    w: &mut W,
    blocks: &[RawBlock],
    compressed: Vec<Result<Vec<u8>>>,
    stats: &mut Stats,
) -> Result<()> {
    for (b, payload) in blocks.iter().zip(compressed) {
        let payload = payload?;
        w.write_all(&(payload.len() as u64).to_le_bytes())?;
        w.write_all(&payload)?;
        trace!(
            reads = b.n_reads(),
            payload = payload.len(),
            "block written"
        );
        stats.reads += b.n_reads() as u64;
        stats.blocks += 1;
        stats.out_bytes += 8 + payload.len() as u64;
    }
    Ok(())
}

/// Write `[u32 len][bytes]`.
fn write_framed<W: Write>(w: &mut W, bytes: &[u8]) -> Result<()> {
    w.write_all(&(bytes.len() as u32).to_le_bytes())?;
    w.write_all(bytes)?;
    Ok(())
}

/// Read a `[u32 len][bytes]` frame, guarding the length allocation.
fn read_framed<R: Read>(r: &mut R) -> Result<Vec<u8>> {
    let mut lb = [0u8; 4];
    r.read_exact(&mut lb)?;
    let len = u32::from_le_bytes(lb) as usize;
    let mut buf = Vec::new();
    buf.try_reserve_exact(len)
        .map_err(|_| Error::Malformed("framed slice too large to allocate"))?;
    buf.resize(len, 0);
    r.read_exact(&mut buf)?;
    Ok(buf)
}

/// Whole-file reorder with GLOBAL clustering (SPRING-style). Buffers all reads,
/// clusters them once, then codes the clustered sequence and the ORIGINAL-order
/// names+quality in independent moderate blocks that fan out across cores. A
/// global permutation (byte-plane rANS) restores original order. This is the
/// `--reorder --keep-order` path; clustering is global so block size is free to
/// be moderate for parallelism without hurting ratio.
fn compress_reordered_whole<R: Read + Send, W: Write>(
    reader: R,
    writer: W,
    params: Params,
) -> Result<Stats> {
    let pool = build_pool(params.threads)?;

    // 1. Buffer every read.
    let mut all = RawBlock::default();
    {
        let mut fq = noodles_fastq::io::Reader::new(BufReader::new(reader));
        let mut rec = noodles_fastq::Record::default();
        while fq.read_record(&mut rec)? != 0 {
            all.push(
                rec.name(),
                rec.description(),
                rec.sequence(),
                rec.quality_scores(),
            );
        }
    }
    let n = all.n_reads();

    // Cumulative byte offsets into the concatenated seq/qual.
    let mut offs = Vec::with_capacity(n + 1);
    let mut acc = 0usize;
    for &l in &all.lens {
        offs.push(acc);
        acc += l as usize;
    }
    offs.push(acc);

    // 2. One global clustering pass over every read.
    let plan = pool.install(|| fqxv_reorder::plan(&all.lens, &all.seq, REORDER_K));

    // 3. Clustered, oriented sequences (parallel copy/revcomp) + flip bitmap.
    let cl_reads: Vec<Vec<u8>> = pool.install(|| {
        plan.order
            .par_iter()
            .map(|&oi| {
                let oi = oi as usize;
                let s = &all.seq[offs[oi]..offs[oi + 1]];
                if plan.flip[oi] {
                    fqxv_reorder::revcomp(s)
                } else {
                    s.to_vec()
                }
            })
            .collect()
    });
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
    let seq_blocks: Vec<Vec<u8>> = pool.install(|| {
        ranges
            .par_iter()
            .map(|&(s, e)| -> Result<Vec<u8>> {
                let refs: Vec<&[u8]> = cl_reads[s..e].iter().map(Vec::as_slice).collect();
                Ok(fqxv_reorder::encode_clustered(
                    &refs,
                    &cl_anchors[s..e],
                    params.seq_order as usize,
                )?)
            })
            .collect::<Result<_>>()
    })?;

    // 5b. Names + quality — ORIGINAL order, per block, in parallel.
    let nq_blocks: Vec<(Vec<u8>, Vec<u8>)> = pool.install(|| {
        ranges
            .par_iter()
            .map(|&(s, e)| -> Result<(Vec<u8>, Vec<u8>)> {
                let headers: Vec<&[u8]> = (s..e).map(|i| all.header(i)).collect();
                let names = fqxv_tokenizer::encode(&headers)?;
                let qual = fqxv_fqzcomp::encode(
                    &all.lens[s..e],
                    &all.qual[offs[s]..offs[e]],
                    params.quality_binning,
                )?;
                Ok((names, qual))
            })
            .collect::<Result<_>>()
    })?;

    // 6. Global permutation (byte-plane split → rANS).
    let mut planes = vec![0u8; n * 4];
    for (i, &x) in plan.order.iter().enumerate() {
        planes[i] = x as u8;
        planes[n + i] = (x >> 8) as u8;
        planes[2 * n + i] = (x >> 16) as u8;
        planes[3 * n + i] = (x >> 24) as u8;
    }
    let perm_c = fqxv_rans::encode(&planes, fqxv_rans::Order::One)?;

    // 7. Write: header, then n / flip / perm / seq blocks / name+qual blocks.
    let mut w = BufWriter::new(writer);
    let flags = FLAG_PLUS_NORMALIZED | FLAG_REORDERED | FLAG_KEEP_ORDER | FLAG_GLOBAL_REORDER;
    w.write_all(&MAGIC)?;
    w.write_all(&FORMAT_VERSION.to_le_bytes())?;
    w.write_all(&[
        params.seq_order,
        binning_tag(params.quality_binning),
        flags,
        1,
    ])?;
    w.write_all(&(n as u64).to_le_bytes())?;
    write_framed(&mut w, &flip_bits)?;
    write_framed(&mut w, &perm_c)?;
    w.write_all(&(ranges.len() as u32).to_le_bytes())?;
    for payload in &seq_blocks {
        write_framed(&mut w, payload)?;
    }
    for (names, qual) in &nq_blocks {
        write_framed(&mut w, names)?;
        write_framed(&mut w, qual)?;
    }
    w.flush()?;

    let out_bytes = (HEADER_LEN
        + 8
        + 4
        + flip_bits.len()
        + 4
        + perm_c.len()
        + 4
        + seq_blocks.iter().map(|p| 4 + p.len()).sum::<usize>()
        + nq_blocks
            .iter()
            .map(|(nm, q)| 8 + nm.len() + q.len())
            .sum::<usize>()) as u64;
    Ok(Stats {
        reads: n as u64,
        blocks: ranges.len() as u64,
        out_bytes,
        group_size: 1,
    })
}

/// Decode a whole-file globally-clustered reorder archive (see
/// [`compress_reordered_whole`]). `r` is positioned just past the header.
fn decode_reordered_whole<R: Read, W: Write>(mut r: R, writer: W, threads: usize) -> Result<Stats> {
    let pool = build_pool(threads)?;
    let mut n_buf = [0u8; 8];
    r.read_exact(&mut n_buf)?;
    let n = u64::from_le_bytes(n_buf) as usize;
    let flip = read_framed(&mut r)?;
    let perm_c = read_framed(&mut r)?;
    let mut nb = [0u8; 4];
    r.read_exact(&mut nb)?;
    let n_blocks = u32::from_le_bytes(nb) as usize;

    let mut seq_payloads: Vec<Vec<u8>> = Vec::with_capacity(n_blocks.min(1 << 20));
    for _ in 0..n_blocks {
        seq_payloads.push(read_framed(&mut r)?);
    }
    let mut nq_payloads: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(n_blocks.min(1 << 20));
    for _ in 0..n_blocks {
        let names = read_framed(&mut r)?;
        let qual = read_framed(&mut r)?;
        nq_payloads.push((names, qual));
    }

    // Global permutation: clustered position j -> original read index.
    let perm: Vec<u32> = {
        let pb = fqxv_rans::decode(&perm_c).map_err(|_| Error::Malformed("bad permutation"))?;
        if pb.len() != n * 4 {
            return Err(Error::Malformed("permutation length mismatch"));
        }
        (0..n)
            .map(|i| u32::from_le_bytes([pb[i], pb[n + i], pb[2 * n + i], pb[3 * n + i]]))
            .collect()
    };

    // Decode both partitions in parallel.
    let seq_dec: Vec<Vec<Vec<u8>>> = pool.install(|| {
        seq_payloads
            .par_iter()
            .map(|p| -> Result<Vec<Vec<u8>>> { Ok(fqxv_reorder::decode_clustered(p)?) })
            .collect::<Result<_>>()
    })?;
    // Per name+quality block: (decoded names, (per-read lengths, quality bytes)).
    type NqBlock = (Vec<Vec<u8>>, (Vec<u32>, Vec<u8>));
    let nq_dec: Vec<NqBlock> = pool.install(|| {
        nq_payloads
            .par_iter()
            .map(|(nm, q)| -> Result<_> {
                Ok((fqxv_tokenizer::decode(nm)?, fqxv_fqzcomp::decode(q)?))
            })
            .collect::<Result<_>>()
    })?;

    // Flatten: clustered reads (clustered order); names/lens/quals (original).
    let mut cl_reads: Vec<Vec<u8>> = Vec::with_capacity(n);
    for blk in seq_dec {
        cl_reads.extend(blk);
    }
    let mut names: Vec<Vec<u8>> = Vec::with_capacity(n);
    let mut lens: Vec<u32> = Vec::with_capacity(n);
    let mut quals: Vec<u8> = Vec::new();
    for (nm, (ls, qs)) in nq_dec {
        names.extend(nm);
        lens.extend(ls);
        quals.extend(qs);
    }
    if cl_reads.len() != n || names.len() != n || lens.len() != n {
        return Err(Error::Malformed("reordered stream length disagreement"));
    }

    // Place each clustered sequence at its original position (un-flipped).
    let mut seq_orig: Vec<Vec<u8>> = vec![Vec::new(); n];
    for (j, mut s) in cl_reads.into_iter().enumerate() {
        if flip.get(j / 8).copied().unwrap_or(0) >> (j % 8) & 1 == 1 {
            s = fqxv_reorder::revcomp(&s);
        }
        let dest = perm[j] as usize;
        *seq_orig
            .get_mut(dest)
            .ok_or(Error::Malformed("permutation out of range"))? = s;
    }

    // Emit records in original order.
    let mut w = BufWriter::new(writer);
    let mut qoff = 0usize;
    for i in 0..n {
        let l = lens[i] as usize;
        let qual = quals
            .get(qoff..qoff + l)
            .ok_or(Error::Malformed("quality underrun"))?;
        qoff += l;
        if seq_orig[i].len() != l {
            return Err(Error::Malformed("reordered sequence length mismatch"));
        }
        let mut rec = Vec::with_capacity(l * 2 + names[i].len() + 8);
        write_record(&mut rec, &names[i], &seq_orig[i], qual);
        w.write_all(&rec)?;
    }
    w.flush()?;
    Ok(Stats {
        reads: n as u64,
        blocks: n_blocks as u64,
        out_bytes: 0,
        group_size: 1,
    })
}

fn compress_block(b: &RawBlock, params: &Params) -> Result<Vec<u8>> {
    if params.reorder {
        compress_block_reordered(b, params)
    } else {
        compress_block_plain(b, params)
    }
}

fn compress_block_plain(b: &RawBlock, params: &Params) -> Result<Vec<u8>> {
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
                || fqxv_seq::encode(&b.lens, &b.seq, params.seq_order as usize),
                || fqxv_fqzcomp::encode(&b.lens, &b.qual, params.quality_binning),
            )
        },
    );
    let (names_c, seq_c, qual_c) = (names_c?, seq_c?, qual_c?);

    let mut out = Vec::with_capacity(16 + names_c.len() + seq_c.len() + qual_c.len());
    out.extend_from_slice(&(b.n_reads() as u32).to_le_bytes());
    for stream in [&names_c, &seq_c, &qual_c] {
        out.extend_from_slice(&(stream.len() as u32).to_le_bytes());
        out.extend_from_slice(stream);
    }
    Ok(out)
}

/// Reorder mode: cluster the block's reads, reverse-complement the flipped ones,
/// then code names (tokenizer), sequence (clustered differential), and quality
/// (fqzcomp) in the clustered order. A flip bitmap (always) and a permutation
/// (only with `keep_order`) let decode restore the reads.
fn compress_block_reordered(b: &RawBlock, params: &Params) -> Result<Vec<u8>> {
    let n = b.n_reads();
    let mut offs = Vec::with_capacity(n + 1);
    let mut acc = 0usize;
    for &l in &b.lens {
        offs.push(acc);
        acc += l as usize;
    }
    offs.push(acc);

    let plan = fqxv_reorder::plan(&b.lens, &b.seq, REORDER_K);

    // Clustered, oriented sequences + flip bitmap — always needed for the
    // differential sequence coder. Flipping (reverse-complement) is a sequence-
    // only trick to make an RC-duplicate byte-identical to its partner.
    let mut r_reads: Vec<Vec<u8>> = Vec::with_capacity(n);
    let mut flip_bits = vec![0u8; n.div_ceil(8)];
    for (j, &oi) in plan.order.iter().enumerate() {
        let oi = oi as usize;
        let s = &b.seq[offs[oi]..offs[oi + 1]];
        if plan.flip[oi] {
            flip_bits[j / 8] |= 1 << (j % 8);
            r_reads.push(fqxv_reorder::revcomp(s));
        } else {
            r_reads.push(s.to_vec());
        }
    }
    let read_refs: Vec<&[u8]> = r_reads.iter().map(Vec::as_slice).collect();
    let r_anchors: Vec<u32> = plan
        .order
        .iter()
        .map(|&oi| plan.anchor[oi as usize])
        .collect();

    // Names and quality: with `keep_order`, code them in ORIGINAL read order —
    // the permutation reunites each clustered sequence with its read, so the
    // tokenizer keeps its positional-delta structure and the quality model its
    // per-position structure (clustering by sequence would scramble both).
    // Without `keep_order`, reads emerge clustered, so names/quality follow.
    let (names_c, seq_c, qual_c) = if params.keep_order {
        let headers = b.header_refs();
        let (n_c, (s_c, q_c)) = rayon::join(
            || fqxv_tokenizer::encode(&headers),
            || {
                rayon::join(
                    || {
                        fqxv_reorder::encode_clustered(
                            &read_refs,
                            &r_anchors,
                            params.seq_order as usize,
                        )
                    },
                    || fqxv_fqzcomp::encode(&b.lens, &b.qual, params.quality_binning),
                )
            },
        );
        (n_c?, s_c?, q_c?)
    } else {
        let mut r_headers: Vec<&[u8]> = Vec::with_capacity(n);
        let mut r_lens: Vec<u32> = Vec::with_capacity(n);
        let mut r_qual: Vec<u8> = Vec::with_capacity(b.qual.len());
        for &oi in &plan.order {
            let oi = oi as usize;
            r_headers.push(b.header(oi));
            r_lens.push(b.lens[oi]);
            let q = &b.qual[offs[oi]..offs[oi + 1]];
            if plan.flip[oi] {
                let mut rq = q.to_vec();
                rq.reverse();
                r_qual.extend_from_slice(&rq);
            } else {
                r_qual.extend_from_slice(q);
            }
        }
        let (n_c, (s_c, q_c)) = rayon::join(
            || fqxv_tokenizer::encode(&r_headers),
            || {
                rayon::join(
                    || {
                        fqxv_reorder::encode_clustered(
                            &read_refs,
                            &r_anchors,
                            params.seq_order as usize,
                        )
                    },
                    || fqxv_fqzcomp::encode(&r_lens, &r_qual, params.quality_binning),
                )
            },
        );
        (n_c?, s_c?, q_c?)
    };
    let perm_c = if params.keep_order {
        // Byte-plane (SoA) split before entropy coding: all byte-0s, then
        // byte-1s, 2s, 3s. For a permutation of n < 2^24 reads the upper planes
        // are near-constant, so per-plane rANS captures them for almost nothing;
        // interleaved LE bytes, by contrast, are near-uniform and barely
        // compress (they used to erase much of reorder's sequence saving).
        let mut planes = vec![0u8; n * 4];
        for (i, &x) in plan.order.iter().enumerate() {
            planes[i] = x as u8;
            planes[n + i] = (x >> 8) as u8;
            planes[2 * n + i] = (x >> 16) as u8;
            planes[3 * n + i] = (x >> 24) as u8;
        }
        fqxv_rans::encode(&planes, fqxv_rans::Order::One)?
    } else {
        Vec::new()
    };

    let mut out = Vec::new();
    out.extend_from_slice(&(n as u32).to_le_bytes());
    out.extend_from_slice(&flip_bits);
    out.extend_from_slice(&(perm_c.len() as u32).to_le_bytes());
    out.extend_from_slice(&perm_c);
    for stream in [&names_c, &seq_c, &qual_c] {
        out.extend_from_slice(&(stream.len() as u32).to_le_bytes());
        out.extend_from_slice(stream);
    }
    Ok(out)
}

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
    // Whole-file globally-clustered reorder uses a distinct layout (two stream
    // partitions + a global permutation), not the per-block block loop.
    if header.flags & FLAG_GLOBAL_REORDER != 0 {
        return decode_reordered_whole(r, writer, threads);
    }
    let reordered = header.flags & FLAG_REORDERED != 0;
    let keep_order = header.flags & FLAG_KEEP_ORDER != 0;
    let mut w = BufWriter::new(writer);

    debug!(
        threads = pool.current_num_threads(),
        batch,
        reordered,
        backend = ?fqxv_rans::Backend::detect(),
        "decompress pool ready"
    );
    let mut stats = Stats::default();
    for_each_block_batch(&mut r, batch, |raw_blocks| {
        debug!(blocks = raw_blocks.len(), "decoding batch");
        let decoded: Vec<Result<(u64, Vec<u8>)>> = pool.install(|| {
            raw_blocks
                .par_iter()
                .map(|b| decode_block(b, reordered, keep_order))
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
    if header.flags & FLAG_REORDERED != 0 {
        return Err(Error::Malformed(
            "reordered archive: use decompress, not split",
        ));
    }
    let g = header.group_size as usize;
    if writers.len() != g {
        return Err(Error::Malformed(
            "output count does not match archive group size",
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

/// Read blocks in batches of `batch`, invoking `f` on each batch.
fn for_each_block_batch<R: Read, F>(r: &mut R, batch: usize, mut f: F) -> Result<()>
where
    F: FnMut(&[Vec<u8>]) -> Result<()>,
{
    loop {
        let mut raw_blocks: Vec<Vec<u8>> = Vec::with_capacity(batch);
        for _ in 0..batch {
            match read_block(r)? {
                Some(block) => raw_blocks.push(block),
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

/// Decoded block streams: (n_reads, names, per-read lengths, sequence, quality).
type BlockParts = (usize, Vec<Vec<u8>>, Vec<u32>, Vec<u8>, Vec<u8>);

/// Decode a block's streams and slice out each read's (name, seq, qual).
fn decode_block_parts(buf: &[u8]) -> Result<BlockParts> {
    let mut c = Cursor::new(buf);
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
    Ok((n_reads, names, seq_lens, seq, qual))
}

fn write_record(out: &mut Vec<u8>, name: &[u8], seq: &[u8], qual: &[u8]) {
    out.push(b'@');
    out.extend_from_slice(name);
    out.push(b'\n');
    out.extend_from_slice(seq);
    out.extend_from_slice(b"\n+\n");
    out.extend_from_slice(qual);
    out.push(b'\n');
}

fn decode_block(buf: &[u8], reordered: bool, keep_order: bool) -> Result<(u64, Vec<u8>)> {
    if reordered {
        return decode_block_reordered(buf, keep_order);
    }
    let (n_reads, names, lens, seq, qual) = decode_block_parts(buf)?;
    let mut out = Vec::with_capacity(seq.len() * 2 + qual.len());
    let mut off = 0usize;
    for i in 0..n_reads {
        let l = lens[i] as usize;
        write_record(&mut out, &names[i], &seq[off..off + l], &qual[off..off + l]);
        off += l;
    }
    Ok((n_reads as u64, out))
}

/// Decode a reorder-mode block: undo the clustering (un-flip reverse-complemented
/// reads, and un-permute when the order was preserved).
fn decode_block_reordered(buf: &[u8], keep_order: bool) -> Result<(u64, Vec<u8>)> {
    let mut c = Cursor::new(buf);
    let n = c.u32()? as usize;
    let flip_bits = c.take(n.div_ceil(8))?.to_vec();
    let perm_c = c.slice_u32()?;
    let (names_s, reads_s, qual_s) = (c.slice_u32()?, c.slice_u32()?, c.slice_u32()?);
    let (names, (reads_r, qual_r)) = rayon::join(
        || fqxv_tokenizer::decode(names_s),
        || {
            rayon::join(
                || fqxv_reorder::decode_clustered(reads_s),
                || fqxv_fqzcomp::decode(qual_s),
            )
        },
    );
    let names = names?;
    let reads = reads_r?;
    let (r_lens, r_qual) = qual_r?;
    if names.len() != n || reads.len() != n || r_lens.len() != n {
        return Err(Error::Malformed(
            "reordered block stream length disagreement",
        ));
    }

    let perm: Vec<u32> = if keep_order {
        let pb = fqxv_rans::decode(perm_c).map_err(|_| Error::Malformed("bad permutation"))?;
        if pb.len() != n * 4 {
            return Err(Error::Malformed("permutation length mismatch"));
        }
        // Reassemble from the four byte planes written by the SoA split above.
        (0..n)
            .map(|i| u32::from_le_bytes([pb[i], pb[n + i], pb[2 * n + i], pb[3 * n + i]]))
            .collect()
    } else {
        Vec::new()
    };

    let mut out = Vec::new();
    if keep_order {
        // Names and quality are in original order. Place each clustered sequence
        // at its original position (un-flipped) via the permutation, then emit
        // records in original order.
        let mut seq_orig: Vec<Vec<u8>> = vec![Vec::new(); n];
        for j in 0..n {
            let mut s = reads[j].clone();
            if flip_bits[j / 8] >> (j % 8) & 1 == 1 {
                s = fqxv_reorder::revcomp(&s);
            }
            let dest = perm[j] as usize;
            *seq_orig
                .get_mut(dest)
                .ok_or(Error::Malformed("permutation out of range"))? = s;
        }
        let mut qoff = 0usize;
        for i in 0..n {
            let l = r_lens[i] as usize;
            let qual = r_qual
                .get(qoff..qoff + l)
                .ok_or(Error::Malformed("quality underrun"))?;
            qoff += l;
            if seq_orig[i].len() != l {
                return Err(Error::Malformed("reordered sequence length mismatch"));
            }
            write_record(&mut out, &names[i], &seq_orig[i], qual);
        }
    } else {
        // Reads emerge in clustered order; names/quality were coded clustered too.
        let mut qoff = 0usize;
        for j in 0..n {
            let l = r_lens[j] as usize;
            let mut seq = reads[j].clone();
            let mut qual = r_qual
                .get(qoff..qoff + l)
                .ok_or(Error::Malformed("quality underrun"))?
                .to_vec();
            qoff += l;
            if flip_bits[j / 8] >> (j % 8) & 1 == 1 {
                seq = fqxv_reorder::revcomp(&seq);
                qual.reverse();
            }
            write_record(&mut out, &names[j], &seq, &qual);
        }
    }
    Ok((n as u64, out))
}

/// Split a grouped block into `g` FASTQ buffers by local read index mod `g`.
fn decode_block_group(buf: &[u8], g: usize) -> Result<(u64, Vec<Vec<u8>>)> {
    let (n_reads, names, lens, seq, qual) = decode_block_parts(buf)?;
    let mut outs = vec![Vec::new(); g];
    let mut off = 0usize;
    for i in 0..n_reads {
        let l = lens[i] as usize;
        write_record(
            &mut outs[i % g],
            &names[i],
            &seq[off..off + l],
            &qual[off..off + l],
        );
        off += l;
    }
    Ok((n_reads as u64, outs))
}

/// Read only the header (cheap) — for discovering the group size before opening
/// split outputs.
pub fn peek<R: Read>(reader: R) -> Result<Info> {
    let mut r = reader;
    let header = read_header(&mut r)?;
    Ok(Info {
        seq_order: header.seq_order,
        quality_binning: header.quality_binning,
        plus_normalized: header.flags & FLAG_PLUS_NORMALIZED != 0,
        group_size: header.group_size,
        reordered: header.flags & FLAG_REORDERED != 0,
        ..Info::default()
    })
}

/// Read a `[u32 len][bytes]` frame's length and skip its bytes without
/// allocating them (for metadata-only scans).
fn skip_framed<R: Read>(r: &mut R) -> Result<usize> {
    let mut lb = [0u8; 4];
    r.read_exact(&mut lb)?;
    let len = u32::from_le_bytes(lb) as usize;
    io::copy(&mut r.by_ref().take(len as u64), &mut io::sink())?;
    Ok(len)
}

/// Read the header and per-stream sizes without decoding block payloads.
pub fn inspect<R: Read>(reader: R) -> Result<Info> {
    let mut r = BufReader::new(reader);
    let header = read_header(&mut r)?;
    let mut info = Info {
        seq_order: header.seq_order,
        quality_binning: header.quality_binning,
        plus_normalized: header.flags & FLAG_PLUS_NORMALIZED != 0,
        group_size: header.group_size,
        reordered: header.flags & FLAG_REORDERED != 0,
        ..Info::default()
    };
    // Whole-file global-cluster layout: [u64 n][flip][perm][u32 n_blocks]
    // [seq blocks][name+qual blocks]. Permutation overhead is charged to seq.
    if header.flags & FLAG_GLOBAL_REORDER != 0 {
        let mut n8 = [0u8; 8];
        r.read_exact(&mut n8)?;
        info.reads = u64::from_le_bytes(n8);
        skip_framed(&mut r)?; // flip bitmap
        info.seq_bytes += skip_framed(&mut r)? as u64; // permutation
        let mut nb = [0u8; 4];
        r.read_exact(&mut nb)?;
        let n_blocks = u32::from_le_bytes(nb) as usize;
        info.blocks = n_blocks as u64;
        for _ in 0..n_blocks {
            info.seq_bytes += skip_framed(&mut r)? as u64;
        }
        for _ in 0..n_blocks {
            info.names_bytes += skip_framed(&mut r)? as u64;
            info.qual_bytes += skip_framed(&mut r)? as u64;
        }
        return Ok(info);
    }
    while let Some(block) = read_block(&mut r)? {
        let mut c = Cursor::new(&block);
        let n = c.u32()? as usize;
        info.reads += n as u64;
        if info.reordered {
            // Reordered blocks carry a flip bitmap and a permutation stream
            // before the three coded streams.
            c.take(n.div_ceil(8))?;
            c.slice_u32()?;
        }
        info.names_bytes += c.slice_u32()?.len() as u64;
        info.seq_bytes += c.slice_u32()?.len() as u64;
        info.qual_bytes += c.slice_u32()?.len() as u64;
        info.blocks += 1;
    }
    Ok(info)
}

// --- header / block framing --------------------------------------------------

struct Header {
    seq_order: u8,
    quality_binning: u8,
    flags: u8,
    group_size: u8,
}

fn read_header<R: Read>(r: &mut R) -> Result<Header> {
    let mut buf = [0u8; HEADER_LEN];
    r.read_exact(&mut buf)?;
    if buf[..4] != MAGIC {
        return Err(Error::BadMagic);
    }
    let ver = u16::from_le_bytes([buf[4], buf[5]]);
    if ver != FORMAT_VERSION {
        return Err(Error::UnsupportedVersion(ver));
    }
    let group_size = buf[9].max(1);
    Ok(Header {
        seq_order: buf[6],
        quality_binning: buf[7],
        flags: buf[8],
        group_size,
    })
}

/// Read one length-prefixed block, or `None` at a clean EOF.
fn read_block<R: Read>(r: &mut R) -> Result<Option<Vec<u8>>> {
    let mut len = [0u8; 8];
    match r.read_exact(&mut len) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e.into()),
    }
    let len = u64::from_le_bytes(len) as usize;
    let mut buf = vec![0u8; len];
    r.read_exact(&mut buf).map_err(|_| Error::Truncated)?;
    Ok(Some(buf))
}

fn build_pool(threads: usize) -> Result<rayon::ThreadPool> {
    // Resolve the effective worker count: 0 means "all available cores", and any
    // explicit request is clamped to what physically exists so we never
    // over-subscribe the pool.
    let available = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let n = if threads == 0 {
        available
    } else {
        threads.min(available)
    };
    rayon::ThreadPoolBuilder::new()
        .num_threads(n)
        .build()
        .map_err(|e| Error::Io(io::Error::other(e.to_string())))
}

fn binning_tag(b: QualityBinning) -> u8 {
    match b {
        QualityBinning::Lossless => 0,
        QualityBinning::Bin8 => 1,
        QualityBinning::Bin4 => 2,
        QualityBinning::Bin2 => 3,
    }
}

struct Cursor<'a> {
    buf: &'a [u8],
    pos: usize,
}

impl<'a> Cursor<'a> {
    fn new(buf: &'a [u8]) -> Self {
        Cursor { buf, pos: 0 }
    }
    fn u32(&mut self) -> Result<u32> {
        let end = self.pos + 4;
        let s = self.buf.get(self.pos..end).ok_or(Error::Truncated)?;
        self.pos = end;
        Ok(u32::from_le_bytes(s.try_into().unwrap()))
    }
    fn slice_u32(&mut self) -> Result<&'a [u8]> {
        let n = self.u32()? as usize;
        self.take(n)
    }
    fn take(&mut self, n: usize) -> Result<&'a [u8]> {
        let end = self.pos + n;
        let s = self.buf.get(self.pos..end).ok_or(Error::Truncated)?;
        self.pos = end;
        Ok(s)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE: &[u8] = b"\
@SRR1.1 INST:1:FC:1:1101:1000:2000 length=8\n\
ACGTACGT\n\
+SRR1.1 INST:1:FC:1:1101:1000:2000 length=8\n\
IIIIFFF#\n\
@SRR1.2 INST:1:FC:1:1101:1005:2050 length=8\n\
NNGGCCTA\n\
+\n\
###IIIFF\n";

    fn compress_bytes(input: &[u8], params: Params) -> Vec<u8> {
        let mut out = Vec::new();
        compress(input, &mut out, params).expect("compress");
        out
    }

    #[test]
    fn roundtrip_normalizes_plus() {
        let archive = compress_bytes(SAMPLE, Params::default());
        let mut fastq = Vec::new();
        decompress(&archive[..], &mut fastq, 1).expect("decompress");
        let expected = b"\
@SRR1.1 INST:1:FC:1:1101:1000:2000 length=8\n\
ACGTACGT\n\
+\n\
IIIIFFF#\n\
@SRR1.2 INST:1:FC:1:1101:1005:2050 length=8\n\
NNGGCCTA\n\
+\n\
###IIIFF\n";
        assert_eq!(fastq, expected);
    }

    fn make_reads(tag: &str, n: usize) -> Vec<u8> {
        let mut v = Vec::new();
        for i in 0..n {
            v.extend_from_slice(format!("@r.{i} {tag}\nACGT\n+\nIIII\n").as_bytes());
        }
        v
    }

    #[test]
    fn paired_roundtrip_splits() {
        let r1 = make_reads("a", 3);
        let r2 = make_reads("b", 3);
        let mut archive = Vec::new();
        let readers: Vec<Box<dyn Read + Send>> =
            vec![Box::new(&r1[..]) as Box<dyn Read + Send>, Box::new(&r2[..])];
        let s = compress_multi(readers, &mut archive, Params::default()).unwrap();
        assert_eq!(s.reads, 6);
        assert_eq!(peek(&archive[..]).unwrap().group_size, 2);

        let (mut o1, mut o2) = (Vec::new(), Vec::new());
        {
            let mut outs: Vec<&mut Vec<u8>> = vec![&mut o1, &mut o2];
            decompress_split(&archive[..], &mut outs, 1).unwrap();
        }
        assert_eq!(o1, r1);
        assert_eq!(o2, r2);
    }

    #[test]
    fn single_cell_four_way_roundtrip() {
        // 10x-style: R1, R2, I1, I2.
        let files: Vec<Vec<u8>> = ["R1", "R2", "I1", "I2"]
            .iter()
            .map(|t| make_reads(t, 5))
            .collect();
        let mut archive = Vec::new();
        let readers: Vec<Box<dyn Read + Send>> = files
            .iter()
            .map(|f| Box::new(&f[..]) as Box<dyn Read + Send>)
            .collect();
        compress_multi(readers, &mut archive, Params::default()).unwrap();
        assert_eq!(peek(&archive[..]).unwrap().group_size, 4);

        let mut outs: Vec<Vec<u8>> = vec![Vec::new(); 4];
        decompress_split(&archive[..], &mut outs, 1).unwrap();
        assert_eq!(outs, files);
    }

    #[test]
    fn grouped_archive_streams_interleaved() {
        let r1 = b"@r.1 a\nACGT\n+\nIIII\n";
        let r2 = b"@r.1 b\nGGGG\n+\n####\n";
        let mut archive = Vec::new();
        let readers: Vec<Box<dyn Read + Send>> =
            vec![Box::new(&r1[..]) as Box<dyn Read + Send>, Box::new(&r2[..])];
        compress_multi(readers, &mut archive, Params::default()).unwrap();
        let mut out = Vec::new();
        decompress(&archive[..], &mut out, 1).unwrap();
        assert_eq!(out, b"@r.1 a\nACGT\n+\nIIII\n@r.1 b\nGGGG\n+\n####\n");
    }

    // Two paired spots, mates interleaved on one stream with /1 /2 names.
    const INTERLEAVED: &[u8] = b"\
@s1/1\nAAAA\n+\nIIII\n\
@s1/2\nTTTT\n+\nFFFF\n\
@s2/1\nCCCC\n+\nIIII\n\
@s2/2\nGGGG\n+\nFFFF\n";

    #[test]
    fn interleaved_stream_forces_pairing_and_splits() {
        let mut archive = Vec::new();
        let s = compress_interleaved(INTERLEAVED, &mut archive, Params::default(), 2).unwrap();
        assert_eq!(s.reads, 4);
        assert_eq!(s.group_size, 2);
        assert_eq!(peek(&archive[..]).unwrap().group_size, 2);

        let (mut o1, mut o2) = (Vec::new(), Vec::new());
        {
            let mut outs: Vec<&mut Vec<u8>> = vec![&mut o1, &mut o2];
            decompress_split(&archive[..], &mut outs, 1).unwrap();
        }
        assert_eq!(o1, b"@s1/1\nAAAA\n+\nIIII\n@s2/1\nCCCC\n+\nIIII\n");
        assert_eq!(o2, b"@s1/2\nTTTT\n+\nFFFF\n@s2/2\nGGGG\n+\nFFFF\n");
    }

    #[test]
    fn interleaved_odd_count_errors() {
        let mut truncated = INTERLEAVED.to_vec();
        truncated.extend_from_slice(b"@s3/1\nACGT\n+\nIIII\n"); // dangling mate
        let err = compress_interleaved(&truncated[..], &mut Vec::new(), Params::default(), 2);
        assert!(matches!(err, Err(Error::Malformed(_))));
    }

    #[test]
    fn auto_detects_interleaved_pairing() {
        let mut archive = Vec::new();
        let s = compress_auto(INTERLEAVED, &mut archive, Params::default()).unwrap();
        assert_eq!(
            s.group_size, 2,
            "paired /1 /2 names should auto-detect as paired"
        );
        assert_eq!(peek(&archive[..]).unwrap().group_size, 2);
    }

    #[test]
    fn auto_leaves_single_end_ungrouped() {
        // Distinct, unpaired names must not be mistaken for mates.
        let single = make_reads("x", 6);
        let mut archive = Vec::new();
        let s = compress_auto(&single[..], &mut archive, Params::default()).unwrap();
        assert_eq!(s.group_size, 1);

        let mut out = Vec::new();
        decompress(&archive[..], &mut out, 1).unwrap();
        assert_eq!(out, single);
    }

    #[test]
    fn unequal_mate_counts_error() {
        let r1 = make_reads("a", 2);
        let r2 = make_reads("b", 1);
        let readers: Vec<Box<dyn Read + Send>> =
            vec![Box::new(&r1[..]) as Box<dyn Read + Send>, Box::new(&r2[..])];
        let err = compress_multi(readers, &mut Vec::new(), Params::default());
        assert!(matches!(err, Err(Error::Malformed(_))));
    }

    #[test]
    fn split_count_mismatch_errors() {
        let archive = compress_bytes(SAMPLE, Params::default()); // group_size 1
        let mut outs: Vec<Vec<u8>> = vec![Vec::new(); 2];
        let err = decompress_split(&archive[..], &mut outs, 1);
        assert!(matches!(err, Err(Error::Malformed(_))));
    }

    #[test]
    fn inspect_reports_streams() {
        let archive = compress_bytes(SAMPLE, Params::default());
        let info = inspect(&archive[..]).expect("inspect");
        assert_eq!(info.reads, 2);
        assert_eq!(info.blocks, 1);
        assert_eq!(info.group_size, 1);
        assert!(info.plus_normalized);
        assert!(info.names_bytes > 0 && info.seq_bytes > 0 && info.qual_bytes > 0);
    }

    // Concatenate every Nth record line (line 4 = quality, line 2 = sequence)
    // across a FASTQ byte stream, in order.
    fn record_line(fastq: &[u8], which: usize) -> Vec<u8> {
        fastq
            .split(|&b| b == b'\n')
            .enumerate()
            .filter(|(i, l)| i % 4 == which && !l.is_empty())
            .flat_map(|(_, l)| l.iter().copied())
            .collect()
    }

    #[test]
    fn lossy_binning_roundtrips_and_reports_tag() {
        for (bin, tag) in [
            (QualityBinning::Bin8, 1u8),
            (QualityBinning::Bin4, 2),
            (QualityBinning::Bin2, 3),
        ] {
            let params = Params {
                quality_binning: bin,
                ..Params::default()
            };
            let archive = compress_bytes(SAMPLE, params);

            // The header tag round-trips through inspect.
            assert_eq!(
                inspect(&archive[..]).expect("inspect").quality_binning,
                tag,
                "info tag for {bin:?}"
            );

            let mut fastq = Vec::new();
            decompress(&archive[..], &mut fastq, 1).expect("decompress");

            // Lossy contract: recovered qualities equal the input qualities passed
            // through the same bin table; bases survive exactly.
            let want: Vec<u8> = record_line(SAMPLE, 3)
                .iter()
                .map(|&b| bin.apply(b))
                .collect();
            assert_eq!(record_line(&fastq, 3), want, "binned qualities for {bin:?}");
            assert_eq!(
                record_line(&fastq, 1),
                record_line(SAMPLE, 1),
                "bases must be exact for {bin:?}"
            );
        }
    }

    #[test]
    fn lossless_default_reports_zero_tag() {
        let archive = compress_bytes(SAMPLE, Params::default());
        assert_eq!(inspect(&archive[..]).unwrap().quality_binning, 0);
    }

    fn dup_rich_input(keep_order_marker: char) -> Vec<u8> {
        // Duplicate-rich single-end reads, including a reverse-complement pair so
        // clustering flips a read (exercises the un-flip path).
        let a = b"ACGTTTGACCGATTGCAACGT";
        let ra = fqxv_reorder::revcomp(a);
        let mut input = Vec::new();
        for i in 0..40u32 {
            let s = match i % 3 {
                0 => a.to_vec(),
                1 => ra.clone(),
                _ => b"TTTTGGGGCCCCAAAATTTTG".to_vec(),
            };
            input.extend_from_slice(format!("@read.{i} {keep_order_marker}\n").as_bytes());
            input.extend_from_slice(&s);
            input.extend_from_slice(format!("\n+\n{}\n", "I".repeat(s.len())).as_bytes());
        }
        input
    }

    fn record_set(fastq: &[u8]) -> Vec<Vec<u8>> {
        let lines: Vec<&[u8]> = fastq.split(|&b| b == b'\n').collect();
        let mut recs: Vec<Vec<u8>> = lines
            .chunks(4)
            .filter(|c| c.len() == 4)
            .map(|c| c.join(&b"\n"[..]))
            .collect();
        recs.sort();
        recs
    }

    #[test]
    fn reorder_keep_order_is_byte_exact() {
        let input = dup_rich_input('d');
        let params = Params {
            reorder: true,
            keep_order: true,
            ..Params::default()
        };
        let mut archive = Vec::new();
        compress(&input[..], &mut archive, params).unwrap();
        assert_eq!(inspect(&archive[..]).unwrap().group_size, 1);
        let mut out = Vec::new();
        decompress(&archive[..], &mut out, 1).unwrap();
        assert_eq!(out, input, "reorder --keep-order must be byte-exact");
    }

    #[test]
    fn reorder_free_preserves_records_as_a_set() {
        let input = dup_rich_input('e');
        let params = Params {
            reorder: true,
            keep_order: false,
            ..Params::default()
        };
        let mut archive = Vec::new();
        compress(&input[..], &mut archive, params).unwrap();
        let mut out = Vec::new();
        decompress(&archive[..], &mut out, 1).unwrap();
        assert_eq!(record_set(&out), record_set(&input));
    }

    #[test]
    fn reorder_rejects_paired() {
        let r: &[u8] = b"@r.1 a\nACGT\n+\nIIII\n";
        let params = Params {
            reorder: true,
            ..Params::default()
        };
        let readers: Vec<Box<dyn Read + Send>> = vec![Box::new(r), Box::new(r)];
        let err = compress_multi(readers, &mut Vec::new(), params);
        assert!(matches!(err, Err(Error::Malformed(_))));
    }

    #[test]
    fn empty_input() {
        let archive = compress_bytes(b"", Params::default());
        let mut fastq = Vec::new();
        let stats = decompress(&archive[..], &mut fastq, 1).unwrap();
        assert_eq!(stats.reads, 0);
        assert!(fastq.is_empty());
    }

    #[test]
    fn bad_magic() {
        let err = decompress(&b"not an fqxv file at all"[..], &mut Vec::new(), 1);
        assert!(matches!(err, Err(Error::BadMagic)));
    }
}
