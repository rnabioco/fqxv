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
//! repeated until the terminator:
//!   [8] block payload length (LE, nonzero)
//!   [4] CRC-32C of the payload (LE) -- verified before decode, so corruption is
//!       caught and localized to one block instead of decoded into garbage
//!   [ ] block payload
//! [8] 0  (zero-length terminator block: a streaming, non-seekable decoder
//!         stops here; seekable readers jump to the footer via the trailer)
//! footer (row-group index — lets `inspect`/random access seek, not scan):
//!   [4] n_row_groups (LE)
//!   per row group: [8] byte_offset (LE, points at the group's length field)
//!                  [4] read_count  (LE)
//!   [8] total_reads (LE)
//!   [4] whole_file_crc (LE)  -- CRC-32C of the archive from byte 0 through the
//!       total_reads field; a one-pass end-to-end integrity check (`verify`)
//!   [4] footer_crc (LE)      -- CRC-32C of the footer body above, checked before
//!       any offset in the index is trusted
//! trailer (fixed, at EOF):
//!   [8] footer_offset (LE)   -- seek straight to the footer
//!   [4] magic "FQXF"
//! block payload:
//!   [4] n_reads (LE)
//!   [4] names_len (LE)  [ ] names   (fqxv-tokenizer)
//!   [4] seq_len   (LE)  [ ] seq     (fqxv-seq)
//!   [4] qual_len  (LE)  [ ] qual    (fqxv-fqzcomp)
//!
//! `--reorder` uses a distinct whole-file, globally-clustered layout (flag bit3),
//! SPRING-style: all reads are clustered in one pass, then the clustered sequence
//! and the names/quality are each coded in independent moderate blocks that fan
//! out across cores. Clustering is global, so block size is free to be moderate
//! for parallelism without hurting ratio. Both `--reorder` modes share this one
//! path — with `--keep-order` (flag bit2) names/quality are coded in ORIGINAL
//! order and a permutation restores it; without it they are coded in CLUSTERED
//! order and no permutation is written. Grouped (paired / single-cell, `G > 1`)
//! input reorders too: the reads are clustered ignoring mate structure, but the
//! permutation reconstructs the original spot interleaving, so `keep_order` is
//! forced on and the archive de-interleaves cleanly on `decompress_split`.
//! Layout after the header:
//! `[8] n  [ ] flip  [ ] perm  [4] n_blocks  [seq block]*n  [ [names][qual] ]*n`
//! (each `[ ]` is a `[u32 len][u32 crc32c][bytes]` frame, CRC-verified on decode;
//! `perm` is empty without keep-order).
//! This layout is self-describing and carries no footer/terminator — decode
//! dispatches on flag bit3 before ever reading a block, so the block-region
//! terminator and footer index above apply only to the plain layout.
//! ```
//!
//! When `G > 1`, reads are interleaved per spot (`m0₀, m1₀, …, m0₁, m1₁, …`).
//! Blocks always hold whole spots and start on member 0, so a block splits back
//! into the `G` files by local read index mod `G`. Interleaving lets the name
//! tokenizer collapse the near-identical mate names and keeps reads from one
//! spot adjacent for the sequence model. [`decompress`] streams interleaved
//! FASTQ (pipe to an aligner); [`decompress_split`] restores the `G` files.

use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};

use rayon::prelude::*;
use tracing::{debug, info, instrument, trace, warn};

use crate::crc::{crc32c, CrcWriter, Hasher};
use crate::{Error, Result, FORMAT_VERSION, MAGIC};
use fqxv_fqzcomp::QualityBinning;

/// Default reads per block. Larger blocks populate the sequence model's contexts
/// better (higher ratio) but reduce parallelism and raise memory.
const DEFAULT_BLOCK_READS: usize = 1 << 20;
/// Raw-sequence byte budget per row group. A group is cut at whichever comes
/// first — `block_reads` reads or this many raw sequence bytes — so long-read
/// (nanopore-style) data does not collapse into one enormous row group that
/// destroys parallelism and random-access granularity and could overflow the
/// `u32` per-stream compressed length. For fixed short reads the read count is
/// the binding limit and this never triggers.
const MAX_BLOCK_SEQ_BYTES: usize = 256 << 20;
const HEADER_LEN: usize = 10;
/// Bytes of CRC-32C appended after a frame's length field (plain block frames)
/// or after a `[u32 len]` framed slice (reorder layout).
const CRC_LEN: usize = 4;
/// Upper bound on a single block payload's declared length. A block holds at most
/// `block_reads` reads and `MAX_BLOCK_SEQ_BYTES` of raw sequence, and the three
/// compressed streams are each smaller than their raw input in the common case,
/// so a real payload is comfortably under this. It exists only so a corrupted
/// length field can't drive a multi-exabyte allocation before the CRC is even
/// checked; anything larger is rejected as malformed.
const MAX_BLOCK_PAYLOAD: u64 = (MAX_BLOCK_SEQ_BYTES as u64) * 8;
/// Magic at the very end of a v1 archive, just after the `[8] footer_offset`
/// back-pointer, so a reader can confirm it found a real footer.
const FOOTER_MAGIC: [u8; 4] = *b"FQXF";
/// Bytes in the fixed EOF trailer: `[8] footer_offset` + `[4] FOOTER_MAGIC`.
const TRAILER_LEN: usize = 12;
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
    /// sequence — captures cross-read duplicate redundancy. Works for grouped
    /// (paired / single-cell) input too; grouped reorder always preserves order.
    pub reorder: bool,
    /// In reorder mode, store a permutation so the original read order is
    /// restored (otherwise reads emerge in clustered order). Forced on for grouped
    /// input (`group_size > 1`), where the permutation reconstructs the spots.
    pub keep_order: bool,
    /// In reorder mode, use the literal-rescue sequence codec: keep every contig
    /// alive and re-attach would-be literals to any contig they overlap (a
    /// k-mer-indexed assembly step). Smaller sequence stream on deep data, at a
    /// higher encode cost. Ignored when `reorder` is false. Decode auto-detects.
    pub rescue: bool,
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
            rescue: false,
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
    /// Whether original read order is restored on decompress (a permutation is
    /// stored). Always true for non-reordered archives; for reordered archives it
    /// reflects the keep-order layout choice.
    pub keep_order: bool,
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
    let mut fqs: Vec<_> = readers
        .into_iter()
        .map(|r| noodles_fastq::io::Reader::new(BufReader::new(r)))
        .collect();
    let mut recs: Vec<_> = (0..g).map(|_| noodles_fastq::Record::default()).collect();

    if params.reorder {
        // Buffer every spot in interleaved order (m0₀, m1₀, …), then globally
        // cluster; the stored permutation restores this spot order on decode, so
        // grouped reorder is order-preserving and de-interleaves cleanly.
        let mut all = RawBlock::default();
        loop {
            for j in 0..g {
                if fqs[j].read_record(&mut recs[j])? == 0 {
                    if j == 0 {
                        return encode_reordered(all, writer, params, g as u8);
                    }
                    return Err(Error::Malformed("inputs have unequal read counts"));
                }
            }
            for r in &recs {
                all.push(r.name(), r.description(), r.sequence(), r.quality_scores());
            }
        }
    }

    // Keep whole spots together: round the block target down to a multiple of g.
    let block_reads = (params.block_reads / g).max(1) * g;
    drive(writer, params, g as u8, |b| {
        // Cut on reads OR the raw-sequence byte budget, whichever comes first;
        // the loop reads whole spots, so a byte cut still lands on a spot
        // boundary. Matches the byte budgeting in `block_ranges`.
        while b.n_reads() < block_reads && b.seq.len() < MAX_BLOCK_SEQ_BYTES {
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
fn parse_chunks(
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
fn block_ranges(
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

/// Write the container header.
fn write_header<W: Write>(w: &mut W, params: &Params, group_size: u8) -> Result<()> {
    // The block layout is always non-reorder — reorder (both keep-order modes)
    // uses the whole-file path, which writes its own header.
    debug_assert!(!params.reorder);
    let flags = FLAG_PLUS_NORMALIZED;
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
    let (chunks, gstart, _n) = parse_chunks(buf, g, &pool)?;

    let mut w = CrcWriter::new(BufWriter::new(writer));
    write_header(&mut w, &params, group_size)?;
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
    let mut w = CrcWriter::new(BufWriter::new(writer));
    write_header(&mut w, &params, group_size)?;

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

/// Row-group index accumulated as blocks are written, serialized into the v1
/// footer. `offset` tracks the current byte position from the start of the file,
/// so each entry records where its row group's `[8] length` field begins.
struct FooterIndex {
    entries: Vec<(u64, u32)>,
    offset: u64,
}

impl FooterIndex {
    fn new() -> Self {
        // Blocks begin right after the fixed header.
        FooterIndex {
            entries: Vec::new(),
            offset: HEADER_LEN as u64,
        }
    }
}

/// Close the block region and append the footer + EOF trailer, returning the
/// number of bytes written (terminator + footer + trailer).
///
/// A zero-length block terminates the region so a streaming, non-seekable
/// decoder stops before the footer; seekable readers ([`inspect`], random
/// access) instead seek to the footer via the trailer's back-pointer. The footer
/// carries two checksums: `whole_file_crc` (the running CRC of every byte written
/// so far, read from the [`CrcWriter`] tee) and `footer_crc` (over the footer
/// body itself), so a reader can both trust the index and detect archive-wide
/// corruption.
fn write_footer<W: Write>(
    w: &mut CrcWriter<W>,
    index: &FooterIndex,
    total_reads: u64,
) -> Result<u64> {
    // Zero-length terminator block (fed through the tee like everything else).
    w.write_all(&0u64.to_le_bytes())?;
    let footer_offset = index.offset + 8;

    let mut body = Vec::with_capacity(4 + index.entries.len() * 12 + 8 + FOOTER_CRC_TAIL);
    body.extend_from_slice(&(index.entries.len() as u32).to_le_bytes());
    for &(off, read_count) in &index.entries {
        body.extend_from_slice(&off.to_le_bytes());
        body.extend_from_slice(&read_count.to_le_bytes());
    }
    body.extend_from_slice(&total_reads.to_le_bytes());
    // Feed the body-so-far through the tee; the tee's CRC is now the digest of
    // everything preceding the whole_file_crc field — exactly what it records.
    w.write_all(&body)?;
    let whole_file_crc = w.crc();
    body.extend_from_slice(&whole_file_crc.to_le_bytes());
    // footer_crc covers the footer body up to but not including itself.
    let footer_crc = crc32c(&body);
    w.write_all(&whole_file_crc.to_le_bytes())?;
    w.write_all(&footer_crc.to_le_bytes())?;

    w.write_all(&footer_offset.to_le_bytes())?;
    w.write_all(&FOOTER_MAGIC)?;
    // body = n_groups..whole_file_crc; +CRC_LEN for footer_crc, +8 terminator.
    Ok(8 + body.len() as u64 + CRC_LEN as u64 + TRAILER_LEN as u64)
}

/// Write a batch's compressed payloads in order, updating `stats` and recording
/// each row group in `index` for the footer.
fn write_blocks<W: Write>(
    w: &mut W,
    blocks: &[RawBlock],
    compressed: Vec<Result<Vec<u8>>>,
    stats: &mut Stats,
    index: &mut FooterIndex,
) -> Result<()> {
    for (b, payload) in blocks.iter().zip(compressed) {
        let payload = payload?;
        index.entries.push((index.offset, b.n_reads() as u32));
        // Frame: [8 payload_len][4 crc32c(payload)][payload].
        w.write_all(&(payload.len() as u64).to_le_bytes())?;
        w.write_all(&crc32c(&payload).to_le_bytes())?;
        w.write_all(&payload)?;
        let framed = (8 + CRC_LEN + payload.len()) as u64;
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

/// Write `[u32 len][u32 crc32c][bytes]`.
fn write_framed<W: Write>(w: &mut W, bytes: &[u8]) -> Result<()> {
    w.write_all(&(bytes.len() as u32).to_le_bytes())?;
    w.write_all(&crc32c(bytes).to_le_bytes())?;
    w.write_all(bytes)?;
    Ok(())
}

/// Read a `[u32 len][u32 crc32c][bytes]` frame, guarding the length allocation and
/// verifying the CRC. `what` names the frame in any corruption error.
fn read_framed<R: Read>(r: &mut R, what: &str) -> Result<Vec<u8>> {
    let mut lb = [0u8; 4];
    r.read_exact(&mut lb)?;
    let len = u32::from_le_bytes(lb) as usize;
    let mut cb = [0u8; CRC_LEN];
    r.read_exact(&mut cb)?;
    let expected = u32::from_le_bytes(cb);
    let mut buf = Vec::new();
    buf.try_reserve_exact(len)
        .map_err(|_| Error::Malformed("framed slice too large to allocate"))?;
    buf.resize(len, 0);
    r.read_exact(&mut buf)?;
    if crc32c(&buf) != expected {
        return Err(Error::Corrupt {
            what: what.to_string(),
        });
    }
    Ok(buf)
}

/// Buffer every record of a single (possibly interleaved) FASTQ stream into one
/// [`RawBlock`], preserving input order. Used by the reorder path, which needs the
/// whole file resident before it can cluster globally.
fn buffer_records(buf: &[u8]) -> Result<RawBlock> {
    let mut all = RawBlock::default();
    let mut fq = noodles_fastq::io::Reader::new(buf);
    let mut rec = noodles_fastq::Record::default();
    while fq.read_record(&mut rec)? != 0 {
        all.push(
            rec.name(),
            rec.description(),
            rec.sequence(),
            rec.quality_scores(),
        );
    }
    Ok(all)
}

/// Buffer a single reader and hand off to [`encode_reordered`] (single-end when
/// `group_size == 1`, or an already-interleaved stream for `group_size > 1`).
fn compress_reordered_whole<R: Read + Send, W: Write>(
    reader: R,
    writer: W,
    params: Params,
    group_size: u8,
) -> Result<Stats> {
    // Buffer every read; input order is preserved, so an interleaved stream stays
    // interleaved and the permutation can restore that spot order on decode.
    let mut all = RawBlock::default();
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
    encode_reordered(all, writer, params, group_size)
}

/// Globally cluster the buffered reads (SPRING-style) and write the whole-file
/// reorder archive: cluster once, then code the clustered sequence in independent
/// moderate blocks that fan out across cores (clustering is global, so block size
/// trades only against parallelism, not ratio). Two modes:
///
/// - `keep_order`: names+quality are coded in ORIGINAL order and a global
///   permutation (byte-plane rANS) restores it, so reads come back byte-exact in
///   input order.
/// - without `keep_order`: names+quality are coded in CLUSTERED order and NO
///   permutation is written, so decode emits reads in clustered order (records
///   preserved as a set) — smaller, but order is not restorable.
///
/// `group_size` is the mate interleaving (1 = single-end). When `group_size > 1`
/// the original spot order *is* the mate interleaving, so the permutation is
/// required to reconstruct it: `keep_order` is forced on regardless of
/// `params.keep_order`, making grouped reorder order-preserving.
fn encode_reordered<W: Write>(
    all: RawBlock,
    writer: W,
    params: Params,
    group_size: u8,
) -> Result<Stats> {
    let g = group_size.max(1);
    let pool = build_pool(params.threads)?;
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
                let anch = &cl_anchors[s..e];
                let order = params.seq_order as usize;
                Ok(if params.rescue {
                    fqxv_reorder::encode_clustered_rescue(&refs, anch, order)?
                } else {
                    fqxv_reorder::encode_clustered(&refs, anch, order)?
                })
            })
            .collect::<Result<_>>()
    })?;

    // 5b. Names, then the keep_order decision, then quality.
    //
    // Reorder has two layouts: keep_order codes names/quality in ORIGINAL order
    // and stores a permutation; otherwise they're coded in CLUSTERED order and no
    // permutation is stored (reads emerge clustered). For single-end input we
    // pick ADAPTIVELY: counter-style names (e.g. SRA `.N N`) delta-code to almost
    // nothing in original order, so a permutation is cheaper than a scrambled
    // name stream; random names are the reverse. Grouped input (the permutation
    // reconstructs spots) and an explicit `params.keep_order` force keep_order.

    // Names coded in ORIGINAL order, per block.
    let names_original = || -> Result<Vec<Vec<u8>>> {
        pool.install(|| {
            ranges
                .par_iter()
                .map(|&(s, e)| {
                    let headers: Vec<&[u8]> = (s..e).map(|i| all.header(i)).collect();
                    Ok(fqxv_tokenizer::encode(&headers)?)
                })
                .collect()
        })
    };
    // Names coded in CLUSTERED order, per block.
    let names_clustered = || -> Result<Vec<Vec<u8>>> {
        pool.install(|| {
            ranges
                .par_iter()
                .map(|&(s, e)| {
                    let headers: Vec<&[u8]> =
                        plan.order[s..e].iter().map(|&oi| all.header(oi as usize)).collect();
                    Ok(fqxv_tokenizer::encode(&headers)?)
                })
                .collect()
        })
    };
    // Global permutation (byte-plane split → rANS), coded with whichever order is
    // smaller (decode auto-detects). The two regimes differ: on huge inputs
    // order-1 wins (real byte-to-byte correlation, its per-context header
    // amortized); on small-to-medium inputs order-1's ~130 KB header dominates
    // and order-0 wins — picking the smaller keeps keep_order efficient at every
    // size. Perm encode is a small fraction of total, so trying both is cheap.
    let encode_perm = || -> Result<Vec<u8>> {
        let mut planes = vec![0u8; n * 4];
        for (i, &x) in plan.order.iter().enumerate() {
            planes[i] = x as u8;
            planes[n + i] = (x >> 8) as u8;
            planes[2 * n + i] = (x >> 16) as u8;
            planes[3 * n + i] = (x >> 24) as u8;
        }
        let o0 = fqxv_rans::encode(&planes, fqxv_rans::Order::Zero)?;
        let o1 = fqxv_rans::encode(&planes, fqxv_rans::Order::One)?;
        Ok(if o0.len() <= o1.len() { o0 } else { o1 })
    };

    let (keep_order, name_blocks, perm_c) = if params.keep_order || g > 1 {
        (true, names_original()?, encode_perm()?)
    } else {
        // Adaptive: keep_order iff original-order names + permutation beat the
        // clustered-order name stream. (Quality's order dependence is second-order
        // and ignored here.) Deterministic — sizes don't depend on thread count.
        let orig = names_original()?;
        let clustered = names_clustered()?;
        let perm = encode_perm()?;
        let keep_bytes = orig.iter().map(Vec::len).sum::<usize>() + perm.len();
        let cluster_bytes = clustered.iter().map(Vec::len).sum::<usize>();
        if keep_bytes < cluster_bytes {
            (true, orig, perm)
        } else {
            (false, clustered, Vec::new())
        }
    };

    // Quality in the chosen order: original for keep_order; otherwise clustered,
    // reversed for flipped reads so bytes line up with the reverse-complemented
    // sequence.
    let qual_blocks: Vec<Vec<u8>> = pool.install(|| {
        ranges
            .par_iter()
            .map(|&(s, e)| -> Result<Vec<u8>> {
                if keep_order {
                    Ok(fqxv_fqzcomp::encode(
                        &all.lens[s..e],
                        &all.qual[offs[s]..offs[e]],
                        params.quality_binning,
                    )?)
                } else {
                    let mut cl_lens: Vec<u32> = Vec::with_capacity(e - s);
                    let mut cl_qual: Vec<u8> = Vec::new();
                    for &oi in &plan.order[s..e] {
                        let oi = oi as usize;
                        cl_lens.push(all.lens[oi]);
                        let q = &all.qual[offs[oi]..offs[oi + 1]];
                        if plan.flip[oi] {
                            cl_qual.extend(q.iter().rev());
                        } else {
                            cl_qual.extend_from_slice(q);
                        }
                    }
                    Ok(fqxv_fqzcomp::encode(&cl_lens, &cl_qual, params.quality_binning)?)
                }
            })
            .collect::<Result<_>>()
    })?;

    let nq_blocks: Vec<(Vec<u8>, Vec<u8>)> = name_blocks.into_iter().zip(qual_blocks).collect();

    // 7. Write: header, then n / flip / perm / seq blocks / name+qual blocks.
    let mut w = BufWriter::new(writer);
    let mut flags = FLAG_PLUS_NORMALIZED | FLAG_REORDERED | FLAG_GLOBAL_REORDER;
    if keep_order {
        flags |= FLAG_KEEP_ORDER;
    }
    w.write_all(&MAGIC)?;
    w.write_all(&FORMAT_VERSION.to_le_bytes())?;
    w.write_all(&[
        params.seq_order,
        binning_tag(params.quality_binning),
        flags,
        g,
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

    // Each framed slice is [4 len][4 crc][bytes]; n_blocks is a bare [4].
    let frame = |len: usize| 4 + CRC_LEN + len;
    let out_bytes = (HEADER_LEN
        + 8
        + frame(flip_bits.len())
        + frame(perm_c.len())
        + 4
        + seq_blocks.iter().map(|p| frame(p.len())).sum::<usize>()
        + nq_blocks
            .iter()
            .map(|(nm, q)| frame(nm.len()) + frame(q.len()))
            .sum::<usize>()) as u64;
    Ok(Stats {
        reads: n as u64,
        blocks: ranges.len() as u64,
        out_bytes,
        group_size: g,
    })
}

/// The decoded whole-file reorder streams, before any un-permutation. `cl_reads`
/// is in clustered order; `names`/`lens`/`quals` are in clustered order without
/// `keep_order` and in original order with it (see [`encode_reordered`]).
struct ReorderStreams {
    n: usize,
    n_blocks: usize,
    flip: Vec<u8>,
    perm_c: Vec<u8>,
    cl_reads: Vec<Vec<u8>>,
    names: Vec<Vec<u8>>,
    lens: Vec<u32>,
    quals: Vec<u8>,
}

/// Read and entropy-decode the whole-file reorder layout. `r` is positioned just
/// past the header. Shared by [`decode_reordered_whole`] and
/// [`decode_reordered_split`].
fn read_reordered_streams<R: Read>(mut r: R, pool: &rayon::ThreadPool) -> Result<ReorderStreams> {
    let mut n_buf = [0u8; 8];
    r.read_exact(&mut n_buf)?;
    let n = u64::from_le_bytes(n_buf) as usize;
    let flip = read_framed(&mut r, "reorder flip bitmap")?;
    let perm_c = read_framed(&mut r, "reorder permutation")?;
    let mut nb = [0u8; 4];
    r.read_exact(&mut nb)?;
    let n_blocks = u32::from_le_bytes(nb) as usize;

    let mut seq_payloads: Vec<Vec<u8>> = Vec::with_capacity(n_blocks.min(1 << 20));
    for i in 0..n_blocks {
        seq_payloads.push(read_framed(&mut r, &format!("reorder sequence block {i}"))?);
    }
    let mut nq_payloads: Vec<(Vec<u8>, Vec<u8>)> = Vec::with_capacity(n_blocks.min(1 << 20));
    for i in 0..n_blocks {
        let names = read_framed(&mut r, &format!("reorder name block {i}"))?;
        let qual = read_framed(&mut r, &format!("reorder quality block {i}"))?;
        nq_payloads.push((names, qual));
    }

    // Decode both partitions in parallel.
    let seq_dec: Vec<Vec<Vec<u8>>> = pool.install(|| {
        seq_payloads
            .par_iter()
            .map(|p| -> Result<Vec<Vec<u8>>> { Ok(fqxv_reorder::decode_clustered_auto(p)?) })
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

    // Flatten the per-block vectors into whole-file streams.
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
    Ok(ReorderStreams {
        n,
        n_blocks,
        flip,
        perm_c,
        cl_reads,
        names,
        lens,
        quals,
    })
}

/// Un-permute a `keep_order` reorder archive: place each clustered sequence
/// (un-flipped) at its original position via the stored permutation, yielding the
/// sequences in original order. Consumes `s.cl_reads`.
fn unpermute_sequences(s: &mut ReorderStreams) -> Result<Vec<Vec<u8>>> {
    let n = s.n;
    let perm: Vec<u32> = {
        let pb = fqxv_rans::decode(&s.perm_c).map_err(|_| Error::Malformed("bad permutation"))?;
        if pb.len() != n * 4 {
            return Err(Error::Malformed("permutation length mismatch"));
        }
        (0..n)
            .map(|i| u32::from_le_bytes([pb[i], pb[n + i], pb[2 * n + i], pb[3 * n + i]]))
            .collect()
    };
    let mut seq_orig: Vec<Vec<u8>> = vec![Vec::new(); n];
    for (j, mut seq) in std::mem::take(&mut s.cl_reads).into_iter().enumerate() {
        if s.flip.get(j / 8).copied().unwrap_or(0) >> (j % 8) & 1 == 1 {
            seq = fqxv_reorder::revcomp(&seq);
        }
        let dest = perm[j] as usize;
        *seq_orig
            .get_mut(dest)
            .ok_or(Error::Malformed("permutation out of range"))? = seq;
    }
    Ok(seq_orig)
}

/// Decode a whole-file globally-clustered reorder archive to interleaved FASTQ on
/// a single writer (see [`compress_reordered_whole`]). `r` is positioned just past
/// the header. `keep_order` (from `FLAG_KEEP_ORDER`) selects the mode: un-permute
/// into original order, or emit in clustered order. `group_size` is recorded in
/// the returned [`Stats`]; grouped archives are always `keep_order`, so their
/// records emerge in original spot-interleaved order.
fn decode_reordered_whole<R: Read, W: Write>(
    r: R,
    writer: W,
    threads: usize,
    keep_order: bool,
    group_size: u8,
) -> Result<Stats> {
    let pool = build_pool(threads)?;
    let mut s = read_reordered_streams(r, &pool)?;
    let n = s.n;
    let n_blocks = s.n_blocks;
    let mut w = BufWriter::new(writer);
    if keep_order {
        // Un-permute, then emit in original order against the original-order
        // names/quality.
        let seq_orig = unpermute_sequences(&mut s)?;
        let mut qoff = 0usize;
        for i in 0..n {
            let l = s.lens[i] as usize;
            let qual = s
                .quals
                .get(qoff..qoff + l)
                .ok_or(Error::Malformed("quality underrun"))?;
            qoff += l;
            if seq_orig[i].len() != l {
                return Err(Error::Malformed("reordered sequence length mismatch"));
            }
            let mut rec = Vec::with_capacity(l * 2 + s.names[i].len() + 8);
            write_record(&mut rec, &s.names[i], &seq_orig[i], qual);
            w.write_all(&rec)?;
        }
    } else {
        // Reads emerge in clustered order; names/quality were coded clustered too.
        // Un-flip the reverse-complemented reads (sequence and quality) to restore
        // each record's original content, then emit in clustered order.
        let mut qoff = 0usize;
        for (j, mut seq) in std::mem::take(&mut s.cl_reads).into_iter().enumerate() {
            let l = s.lens[j] as usize;
            let mut qual = s
                .quals
                .get(qoff..qoff + l)
                .ok_or(Error::Malformed("quality underrun"))?
                .to_vec();
            qoff += l;
            if seq.len() != l {
                return Err(Error::Malformed("reordered sequence length mismatch"));
            }
            if s.flip.get(j / 8).copied().unwrap_or(0) >> (j % 8) & 1 == 1 {
                seq = fqxv_reorder::revcomp(&seq);
                qual.reverse();
            }
            let mut rec = Vec::with_capacity(l * 2 + s.names[j].len() + 8);
            write_record(&mut rec, &s.names[j], &seq, &qual);
            w.write_all(&rec)?;
        }
    }
    w.flush()?;
    Ok(Stats {
        reads: n as u64,
        blocks: n_blocks as u64,
        out_bytes: 0,
        group_size: group_size.max(1),
    })
}

/// Decode a grouped whole-file reorder archive, splitting the reads back into `g`
/// writers by their per-spot member. Only valid for `keep_order` archives (the
/// permutation reconstructs the mate interleaving); the caller guarantees this.
/// Record `i` in restored original order belongs to member `i % g`.
fn decode_reordered_split<R: Read, W: Write>(
    r: R,
    writers: &mut [W],
    threads: usize,
    g: usize,
) -> Result<Stats> {
    let pool = build_pool(threads)?;
    let mut s = read_reordered_streams(r, &pool)?;
    let n = s.n;
    let n_blocks = s.n_blocks;
    let seq_orig = unpermute_sequences(&mut s)?;
    let mut bufs: Vec<BufWriter<&mut W>> = writers.iter_mut().map(BufWriter::new).collect();
    let mut stats = Stats {
        reads: n as u64,
        blocks: n_blocks as u64,
        group_size: g as u8,
        ..Stats::default()
    };
    let mut qoff = 0usize;
    for i in 0..n {
        let l = s.lens[i] as usize;
        let qual = s
            .quals
            .get(qoff..qoff + l)
            .ok_or(Error::Malformed("quality underrun"))?;
        qoff += l;
        if seq_orig[i].len() != l {
            return Err(Error::Malformed("reordered sequence length mismatch"));
        }
        let mut rec = Vec::with_capacity(l * 2 + s.names[i].len() + 8);
        write_record(&mut rec, &s.names[i], &seq_orig[i], qual);
        bufs[i % g].write_all(&rec)?;
        stats.out_bytes += rec.len() as u64;
    }
    for b in &mut bufs {
        b.flush()?;
    }
    Ok(stats)
}

/// Code one non-reorder block: names (tokenizer), sequence (order-k), and quality
/// (fqzcomp), each length-prefixed. Reorder uses the whole-file path instead.
fn compress_block(b: &RawBlock, params: &Params) -> Result<Vec<u8>> {
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
        return decode_reordered_whole(r, writer, threads, keep_order, header.group_size);
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
            return decode_reordered_split(r, writers, threads, g);
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

/// Verify an archive's integrity without materializing the decoded FASTQ.
///
/// For the plain layout this checks the footer's own CRC, then re-hashes the
/// archive prefix and compares against the stored whole-file CRC-32C — a single
/// linear pass that catches any corruption in the header, any block payload,
/// framing, or the index. For the globally-clustered reorder layout (which has no
/// footer) it decodes into a sink, so every frame CRC and cross-stream length
/// check is exercised. Returns `Ok(())` iff the archive is intact.
pub fn verify<R: Read + Seek>(reader: R) -> Result<()> {
    let mut r = BufReader::new(reader);
    let header = read_header(&mut r)?;
    if header.flags & FLAG_GLOBAL_REORDER != 0 {
        // No footer/whole-file digest here; decoding drives every frame CRC.
        let keep_order = header.flags & FLAG_KEEP_ORDER != 0;
        decode_reordered_whole(r, io::sink(), 0, keep_order, header.group_size)?;
        return Ok(());
    }
    let footer = read_footer(&mut r)?;
    r.seek(SeekFrom::Start(0))?;
    let mut hasher = Hasher::new();
    let mut remaining = footer.covered_len;
    let mut buf = [0u8; 1 << 16];
    while remaining > 0 {
        let want = remaining.min(buf.len() as u64) as usize;
        r.read_exact(&mut buf[..want])?;
        hasher.update(&buf[..want]);
        remaining -= want as u64;
    }
    if hasher.finalize() != footer.whole_file_crc {
        return Err(Error::Corrupt {
            what: "archive (whole-file crc)".to_string(),
        });
    }
    Ok(())
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

/// Read blocks in batches of `batch`, invoking `f` on each batch.
fn for_each_block_batch<R: Read, F>(r: &mut R, batch: usize, mut f: F) -> Result<()>
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

fn decode_block(buf: &[u8]) -> Result<(u64, Vec<u8>)> {
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
fn read_slices<'a>(
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
fn decode_block_group(buf: &[u8], g: usize) -> Result<(u64, Vec<Vec<u8>>)> {
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
        keep_order: header.flags & FLAG_REORDERED == 0
            || header.flags & FLAG_KEEP_ORDER != 0,
        ..Info::default()
    })
}

/// Read a `[u32 len][u32 crc][bytes]` frame's length and skip the CRC + bytes
/// without allocating them (for metadata-only scans; the CRC is not verified
/// since the payload is discarded). Returns the payload length.
fn skip_framed<R: Read>(r: &mut R) -> Result<usize> {
    let mut lb = [0u8; 4];
    r.read_exact(&mut lb)?;
    let len = u32::from_le_bytes(lb) as usize;
    io::copy(
        &mut r.by_ref().take((CRC_LEN + len) as u64),
        &mut io::sink(),
    )?;
    Ok(len)
}

/// Read the header and per-stream sizes without decoding block payloads.
///
/// For the plain and per-block-reorder layouts this uses the v1 footer: it reads
/// the row-group index (O(1) for reads/blocks/total), then seeks to each row
/// group and reads only its small block header — the `n_reads`, reorder preamble,
/// and three stream length prefixes — skipping every coded payload. So the cost
/// is O(row groups) tiny reads rather than O(archive bytes). The globally
/// clustered layout has no footer and is scanned sequentially as before.
pub fn inspect<R: Read + Seek>(reader: R) -> Result<Info> {
    let mut r = BufReader::new(reader);
    let header = read_header(&mut r)?;
    let mut info = Info {
        seq_order: header.seq_order,
        quality_binning: header.quality_binning,
        plus_normalized: header.flags & FLAG_PLUS_NORMALIZED != 0,
        group_size: header.group_size,
        reordered: header.flags & FLAG_REORDERED != 0,
        keep_order: header.flags & FLAG_REORDERED == 0
            || header.flags & FLAG_KEEP_ORDER != 0,
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

    // Prefer the footer's O(row groups) index. If it's unreadable — a truncated
    // download loses the EOF trailer, corruption can fail its CRC — fall back to
    // scanning block frames forward from the header so a partial file still
    // reports what it contains.
    match read_footer(&mut r) {
        Ok(footer) => {
            info.reads = footer.total_reads;
            info.blocks = footer.groups.len() as u64;
            for &(off, _) in &footer.groups {
                // Position past the [8] payload length and [4] block CRC, at the
                // block header, then walk its stream frames.
                r.seek(SeekFrom::Start(off + 8 + CRC_LEN as u64))?;
                scan_block_header(&mut r, info.reordered, &mut info)?;
            }
        }
        Err(_) => scan_blocks_sequentially(&mut r, &mut info)?,
    }
    Ok(info)
}

/// Walk one block's header at the current position — `[4] n_reads`, an optional
/// reorder preamble, and the three `[4 len][bytes]` stream frames — accumulating
/// per-stream sizes into `info` and seeking past each payload. Returns the
/// block's read count. Leaves the cursor at the end of the block's payload.
fn scan_block_header<R: Read + Seek>(r: &mut R, reordered: bool, info: &mut Info) -> Result<u64> {
    let n = u64::from(read_u32(r)?);
    if reordered {
        r.seek(SeekFrom::Current(n.div_ceil(8) as i64))?; // flip bitmap
        let perm_len = read_u32(r)? as i64;
        r.seek(SeekFrom::Current(perm_len))?;
    }
    for bytes in [
        &mut info.names_bytes,
        &mut info.seq_bytes,
        &mut info.qual_bytes,
    ] {
        let len = read_u32(r)?;
        *bytes += u64::from(len);
        r.seek(SeekFrom::Current(i64::from(len)))?;
    }
    Ok(n)
}

/// Footer-less fallback for [`inspect`]: scan block frames forward from the
/// header until the terminator, a clean EOF, or the first structurally
/// implausible frame (a truncated download's boundary). Best-effort — the block
/// CRCs are not checked here, only the framing is followed.
fn scan_blocks_sequentially<R: Read + Seek>(r: &mut R, info: &mut Info) -> Result<()> {
    let mut off = HEADER_LEN as u64;
    loop {
        r.seek(SeekFrom::Start(off))?;
        let mut lenb = [0u8; 8];
        if r.read_exact(&mut lenb).is_err() {
            break; // clean EOF at a frame boundary
        }
        let plen = u64::from_le_bytes(lenb);
        if plen == 0 || plen > MAX_BLOCK_PAYLOAD {
            break; // terminator, or a length too large to be a real frame
        }
        r.seek(SeekFrom::Start(off + 8 + CRC_LEN as u64))?;
        let n = match scan_block_header(r, info.reordered, info) {
            Ok(n) => n,
            Err(_) => break, // ran off the end of a truncated block
        };
        info.reads += n;
        info.blocks += 1;
        off += 8 + CRC_LEN as u64 + plen;
    }
    Ok(())
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

/// Read one length-prefixed, CRC-checked block, or `None` at the terminator /
/// a clean EOF. `index` names the block in any corruption error.
///
/// A zero-length block is the terminator that separates the block region from
/// the footer, so a streaming (non-seekable) decoder stops here without reading
/// into the footer. A clean EOF (no length at all) is also treated as the end,
/// which keeps truncated pre-footer streams decoding what they can. The frame is
/// `[8 payload_len][4 crc32c(payload)][payload]`; the CRC is verified before the
/// payload is handed to the entropy decoders so corruption surfaces as a clean
/// [`Error::Corrupt`] rather than garbage output.
fn read_block<R: Read>(r: &mut R, index: u64) -> Result<Option<Vec<u8>>> {
    let mut len = [0u8; 8];
    match r.read_exact(&mut len) {
        Ok(()) => {}
        Err(e) if e.kind() == io::ErrorKind::UnexpectedEof => return Ok(None),
        Err(e) => return Err(e.into()),
    }
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

/// The footer index: per row group `(byte_offset, read_count)`, the total read
/// count, and the whole-archive CRC-32C, located via the EOF trailer's
/// back-pointer. Only the plain and per-block-reorder layouts carry a footer (the
/// globally-clustered layout is self-describing — see the module docs).
///
/// On-disk body (at `footer_offset`, i.e. just past the terminator):
/// `[4 n_groups] [8 off][4 read_count]* [8 total_reads] [4 whole_file_crc] [4 footer_crc]`.
/// `footer_crc` covers the body up to (not including) itself, so the index can be
/// trusted without rereading the whole archive; `whole_file_crc` covers every
/// byte from the header through `total_reads` and is checked only by the
/// whole-archive verify/recover path.
struct Footer {
    /// `(byte offset of the row group's [8] length field, its read count)`.
    groups: Vec<(u64, u32)>,
    total_reads: u64,
    /// CRC-32C of the archive from byte 0 through the `total_reads` field.
    whole_file_crc: u32,
    /// Number of leading bytes that `whole_file_crc` covers (the offset of the
    /// `whole_file_crc` field itself), so a verifier can re-hash exactly that
    /// prefix.
    covered_len: u64,
}

/// Bytes appended to the footer body after `total_reads`: `[4 whole_file_crc]
/// [4 footer_crc]`.
const FOOTER_CRC_TAIL: usize = 8;

/// Read the footer by seeking to the EOF trailer and following its back-pointer,
/// then verify the footer's own CRC before trusting any offset in it.
fn read_footer<R: Read + Seek>(r: &mut R) -> Result<Footer> {
    let end = r.seek(SeekFrom::End(0))?;
    if end < (HEADER_LEN + TRAILER_LEN) as u64 {
        return Err(Error::Truncated);
    }
    let mut trailer = [0u8; TRAILER_LEN];
    r.seek(SeekFrom::End(-(TRAILER_LEN as i64)))?;
    r.read_exact(&mut trailer)?;
    if trailer[8..] != FOOTER_MAGIC {
        return Err(Error::Malformed("missing footer trailer magic"));
    }
    let footer_offset = u64::from_le_bytes(trailer[..8].try_into().unwrap());
    let body_end = end - TRAILER_LEN as u64;
    // Body must hold at least n_groups(4) + total_reads(8) + the crc tail(8).
    if footer_offset < HEADER_LEN as u64
        || footer_offset + (4 + 8 + FOOTER_CRC_TAIL as u64) > body_end
    {
        return Err(Error::Malformed("footer offset out of range"));
    }
    // Read the whole footer body in one shot, then parse and CRC-check it from
    // memory — cheaper and simpler than incremental reads, and lets the CRC guard
    // the index before any offset is dereferenced.
    let body_len = (body_end - footer_offset) as usize;
    let mut body = Vec::new();
    body.try_reserve_exact(body_len)
        .map_err(|_| Error::Malformed("footer too large to allocate"))?;
    body.resize(body_len, 0);
    r.seek(SeekFrom::Start(footer_offset))?;
    r.read_exact(&mut body)?;

    let (covered, footer_crc_bytes) = body.split_at(body_len - CRC_LEN);
    let footer_crc = u32::from_le_bytes(footer_crc_bytes.try_into().unwrap());
    if crc32c(covered) != footer_crc {
        return Err(Error::Corrupt {
            what: "footer".to_string(),
        });
    }

    let mut c = Cursor::new(covered);
    let n_groups = c.u32()? as usize;
    // `covered` = [4 n_groups][12*n_groups][8 total_reads][4 whole_file_crc]; the
    // CRC just passed already implies a self-consistent length, but bound the
    // allocation independently in case of a hash collision.
    let max_groups = covered.len().saturating_sub(4 + 8 + CRC_LEN) / 12;
    if n_groups > max_groups {
        return Err(Error::Malformed("footer group count exceeds footer size"));
    }
    let mut groups = Vec::with_capacity(n_groups);
    for _ in 0..n_groups {
        let off = u64::from_le_bytes(c.take(8)?.try_into().unwrap());
        let rc = u32::from_le_bytes(c.take(4)?.try_into().unwrap());
        if off < HEADER_LEN as u64 || off >= footer_offset {
            return Err(Error::Malformed("footer row-group offset out of range"));
        }
        groups.push((off, rc));
    }
    let total_reads = u64::from_le_bytes(c.take(8)?.try_into().unwrap());
    let whole_file_crc = c.u32()?;
    // whole_file_crc covers everything up to its own field: the archive prefix
    // plus the footer body through total_reads.
    let covered_len = footer_offset + (body_len - FOOTER_CRC_TAIL) as u64;
    Ok(Footer {
        groups,
        total_reads,
        whole_file_crc,
        covered_len,
    })
}

/// Read a little-endian `u32` from `r`.
fn read_u32<R: Read>(r: &mut R) -> Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
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
        let info = inspect(io::Cursor::new(&archive[..])).expect("inspect");
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
                inspect(io::Cursor::new(&archive[..]))
                    .expect("inspect")
                    .quality_binning,
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
        assert_eq!(
            inspect(io::Cursor::new(&archive[..]))
                .unwrap()
                .quality_binning,
            0
        );
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
            rescue: false,
            ..Params::default()
        };
        let mut archive = Vec::new();
        compress(&input[..], &mut archive, params).unwrap();
        assert_eq!(
            inspect(io::Cursor::new(&archive[..])).unwrap().group_size,
            1
        );
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
            rescue: false,
            ..Params::default()
        };
        let mut archive = Vec::new();
        compress(&input[..], &mut archive, params).unwrap();
        let mut out = Vec::new();
        decompress(&archive[..], &mut out, 1).unwrap();
        assert_eq!(record_set(&out), record_set(&input));
    }

    #[test]
    fn reorder_rescue_preserves_records_as_a_set() {
        // The literal-rescue sequence codec must round-trip through the container
        // (decode auto-detects the version byte). Multi-thread to exercise the
        // parallel per-block encode path.
        let input = dup_rich_input('e');
        for threads in [1usize, 4] {
            let params = Params {
                reorder: true,
                keep_order: false,
                rescue: true,
                threads,
                ..Params::default()
            };
            let mut archive = Vec::new();
            compress(&input[..], &mut archive, params).unwrap();
            let mut out = Vec::new();
            decompress(&archive[..], &mut out, 1).unwrap();
            assert_eq!(record_set(&out), record_set(&input), "threads={threads}");
        }
    }

    /// FASTQ of `n` overlapping windows (length `win`) of a fixed pseudo-random
    /// reference, emitted in a SHUFFLED order with header text from `name(i)`.
    /// The windows share minimizers so clustering re-groups them; because file
    /// order is shuffled, clustered order differs from file order — so a
    /// positional counter in the name scrambles under clustering (the case where
    /// keep-order pays off). Bare `+`, so a keep-order archive round-trips
    /// byte-for-byte.
    fn windowed_input(name: impl Fn(usize) -> String, n: usize, win: usize) -> Vec<u8> {
        let bases = b"ACGT";
        let mut x = 0x1234_5678u32;
        let mut lcg = || {
            x = x.wrapping_mul(1_103_515_245).wrapping_add(12_345);
            x
        };
        let mut refseq = Vec::with_capacity(n + win);
        for _ in 0..n + win {
            refseq.push(bases[((lcg() >> 16) & 3) as usize]);
        }
        // Window starts, Fisher-Yates shuffled so file order != clustered order.
        let mut starts: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            starts.swap(i, lcg() as usize % (i + 1));
        }
        let mut v = Vec::new();
        for i in 0..n {
            v.extend_from_slice(name(i).as_bytes());
            v.push(b'\n');
            let s = starts[i];
            v.extend_from_slice(&refseq[s..s + win]);
            v.extend_from_slice(b"\n+\n");
            v.extend(std::iter::repeat_n(b'I', win));
            v.push(b'\n');
        }
        v
    }

    #[test]
    fn adaptive_keeps_order_for_counter_names() {
        // Counter-style names (the `.N N` pattern) delta-code to almost nothing in
        // original order, so the permutation is cheaper than the scrambled-counter
        // clustered-order stream: adaptive should keep order — and then restore it
        // byte-for-byte.
        let input = windowed_input(|i| format!("@read.{i} {i}"), 2000, 40);
        let params = Params {
            reorder: true,
            ..Params::default()
        };
        let mut archive = Vec::new();
        compress(&input[..], &mut archive, params).unwrap();
        assert!(
            peek(&archive[..]).unwrap().keep_order,
            "counter names should trigger keep_order"
        );
        let mut out = Vec::new();
        decompress(&archive[..], &mut out, 1).unwrap();
        assert_eq!(out, input, "keep_order must restore original order exactly");
    }

    #[test]
    fn adaptive_drops_order_for_random_names() {
        // Avalanched (splitmix64) names look i.i.d., so they carry no order
        // structure: original- and clustered-order coding cost the same and the
        // permutation is pure overhead — adaptive should NOT keep order.
        let splitmix = |i: usize| -> u64 {
            let mut z = (i as u64).wrapping_add(0x9E37_79B9_7F4A_7C15);
            z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            z ^ (z >> 31)
        };
        let input = windowed_input(|i| format!("@{:016x}", splitmix(i)), 2000, 40);
        let params = Params {
            reorder: true,
            ..Params::default()
        };
        let mut archive = Vec::new();
        compress(&input[..], &mut archive, params).unwrap();
        assert!(
            !peek(&archive[..]).unwrap().keep_order,
            "random names should not keep order (permutation is pure overhead)"
        );
        let mut out = Vec::new();
        decompress(&archive[..], &mut out, 1).unwrap();
        assert_eq!(record_set(&out), record_set(&input));
    }

    /// Duplicate-rich reads for one mate of a paired set (`n` spots), sharing
    /// sequences across spots so clustering has real work — including a
    /// reverse-complement so the flip path is exercised.
    fn dup_rich_mate(mate: u8, n: u32) -> Vec<u8> {
        let a = b"ACGTTTGACCGATTGCAACGT";
        let ra = fqxv_reorder::revcomp(a);
        let mut v = Vec::new();
        for i in 0..n {
            let s = match (i + mate as u32) % 3 {
                0 => a.to_vec(),
                1 => ra.clone(),
                _ => b"TTTTGGGGCCCCAAAATTTTG".to_vec(),
            };
            v.extend_from_slice(format!("@spot.{i}/{mate}\n").as_bytes());
            v.extend_from_slice(&s);
            v.extend_from_slice(format!("\n+\n{}\n", "I".repeat(s.len())).as_bytes());
        }
        v
    }

    fn paired_readers<'a>(r1: &'a [u8], r2: &'a [u8]) -> Vec<Box<dyn Read + Send + 'a>> {
        vec![Box::new(r1) as Box<dyn Read + Send>, Box::new(r2)]
    }

    #[test]
    fn reorder_paired_preserves_order_and_splits() {
        let r1 = dup_rich_mate(1, 30);
        let r2 = dup_rich_mate(2, 30);
        let params = Params {
            reorder: true,
            ..Params::default()
        };
        let mut archive = Vec::new();
        compress_multi(paired_readers(&r1, &r2), &mut archive, params).unwrap();

        // Grouped reorder records the real group size and stays a reorder archive.
        assert_eq!(peek(&archive[..]).unwrap().group_size, 2);
        assert!(inspect(io::Cursor::new(&archive[..])).unwrap().reordered);

        // Interleaved decode matches a plain (non-reorder) archive byte-for-byte,
        // i.e. the permutation fully restored the original spot interleaving.
        let mut plain = Vec::new();
        compress_multi(paired_readers(&r1, &r2), &mut plain, Params::default()).unwrap();
        let (mut expected, mut got) = (Vec::new(), Vec::new());
        decompress(&plain[..], &mut expected, 1).unwrap();
        decompress(&archive[..], &mut got, 1).unwrap();
        assert_eq!(got, expected, "grouped reorder must restore spot order");

        // Split decode reconstructs each mate file exactly.
        let (mut o1, mut o2) = (Vec::new(), Vec::new());
        {
            let mut outs: Vec<&mut Vec<u8>> = vec![&mut o1, &mut o2];
            decompress_split(&archive[..], &mut outs, 1).unwrap();
        }
        assert_eq!(o1, r1);
        assert_eq!(o2, r2);
    }

    #[test]
    fn reorder_paired_is_thread_count_deterministic() {
        let r1 = dup_rich_mate(1, 50);
        let r2 = dup_rich_mate(2, 50);
        let mut archives = Vec::new();
        for threads in [1usize, 4] {
            let params = Params {
                reorder: true,
                threads,
                ..Params::default()
            };
            let mut archive = Vec::new();
            compress_multi(paired_readers(&r1, &r2), &mut archive, params).unwrap();
            archives.push(archive);
        }
        assert_eq!(
            archives[0], archives[1],
            "reorder output must not vary by threads"
        );
    }

    #[test]
    fn reorder_interleaved_single_stream_splits() {
        // One already-interleaved stream (as `sracha get -Z` emits), reordered.
        let mut stream = Vec::new();
        for i in 0..20u32 {
            for mate in 1..=2u8 {
                let s = if (i + mate as u32).is_multiple_of(2) {
                    "ACGTACGTAC"
                } else {
                    "TTTTGGGGCC"
                };
                stream.extend_from_slice(
                    format!("@spot.{i}/{mate}\n{s}\n+\nIIIIIIIIII\n").as_bytes(),
                );
            }
        }
        let params = Params {
            reorder: true,
            ..Params::default()
        };
        let mut archive = Vec::new();
        compress_interleaved(&stream[..], &mut archive, params, 2).unwrap();
        assert_eq!(peek(&archive[..]).unwrap().group_size, 2);
        let mut out = Vec::new();
        decompress(&archive[..], &mut out, 1).unwrap();
        assert_eq!(out, stream, "interleaved reorder must be byte-exact");
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

    // --- v1 footer: index, determinism --------------------------------------

    #[test]
    fn block_ranges_cuts_on_reads_bytes_and_spots() {
        // 10 reads × 100 bp, one parse chunk.
        let mut fq = Vec::new();
        for i in 0..10 {
            fq.extend_from_slice(
                format!("@r{i}\n{}\n+\n{}\n", "A".repeat(100), "I".repeat(100)).as_bytes(),
            );
        }
        let chunks = vec![parse_chunk(&fq, 0, fq.len()).unwrap()];

        // Byte budget binds first: cut every 3 reads (300 B ≥ 250 B).
        assert_eq!(
            block_ranges(&chunks, 1000, 250, 1),
            vec![(0, 3), (3, 6), (6, 9), (9, 10)]
        );
        // Read budget binds first.
        assert_eq!(
            block_ranges(&chunks, 4, usize::MAX, 1),
            vec![(0, 4), (4, 8), (8, 10)]
        );
        // A tiny byte budget still only cuts on whole spots (g = 2).
        assert_eq!(
            block_ranges(&chunks, 1000, 1, 2),
            vec![(0, 2), (2, 4), (4, 6), (6, 8), (8, 10)]
        );
    }

    #[test]
    fn archive_is_deterministic_across_threads() {
        let input = make_reads("y", 500);
        let mk = |threads| {
            compress_bytes(
                &input,
                Params {
                    block_reads: 32,
                    threads,
                    ..Params::default()
                },
            )
        };
        // Byte-identical (header, blocks, and footer offsets) regardless of pool.
        assert_eq!(mk(1), mk(4));
    }

    #[test]
    fn inspect_falls_back_without_trailer() {
        let archive = compress_bytes(SAMPLE, Params::default());
        // Drop the trailing "FQXF" magic so the footer can't be located — a
        // partial download loses the EOF trailer this way. inspect must fall back
        // to a forward scan and still report the intact blocks rather than error.
        let truncated = &archive[..archive.len() - 4];
        let info = inspect(io::Cursor::new(truncated)).expect("fallback scan");
        assert_eq!(info.reads, 2);
        assert_eq!(info.blocks, 1);
    }

    #[test]
    fn ragged_lengths_roundtrip_multiblock() {
        // Mixed read lengths (10..=310 bp) exercise variable-length framing; a
        // small block target spreads them over several row groups.
        let mut input = Vec::new();
        for i in 0..30usize {
            let len = 10 + (i % 7) * 50;
            let seq: String = "ACGT".chars().cycle().take(len).collect();
            input.extend_from_slice(
                format!("@read.{i}\n{seq}\n+\n{}\n", "I".repeat(len)).as_bytes(),
            );
        }
        let params = Params {
            block_reads: 5,
            ..Params::default()
        };
        let archive = compress_bytes(&input, params);
        let mut out = Vec::new();
        decompress(&archive[..], &mut out, 1).unwrap();
        assert_eq!(out, input, "ragged variable-length reads must round-trip");
    }

    // --- integrity: CRC detection, recovery, truncation --------------------

    /// Build a multi-block archive of `n` uniform reads (small block target).
    fn multiblock_archive(n: usize, block_reads: usize) -> Vec<u8> {
        let input = make_reads("x", n);
        compress_bytes(
            &input,
            Params {
                block_reads,
                ..Params::default()
            },
        )
    }

    #[test]
    fn verify_accepts_intact_archive() {
        let archive = multiblock_archive(40, 8);
        verify(io::Cursor::new(&archive)).expect("intact archive verifies");
    }

    #[test]
    fn verify_rejects_payload_bit_flip() {
        let mut archive = multiblock_archive(40, 8);
        // Header(10) + [8 len][4 crc]; the first payload byte is at offset 22.
        archive[HEADER_LEN + 8 + CRC_LEN] ^= 0x01;
        let err = verify(io::Cursor::new(&archive)).unwrap_err();
        assert!(matches!(err, Error::Corrupt { .. }), "got {err:?}");
    }

    #[test]
    fn verify_rejects_footer_bit_flip() {
        let mut archive = multiblock_archive(40, 8);
        // Flip a byte inside the footer body (just before the EOF trailer).
        let i = archive.len() - TRAILER_LEN - 1;
        archive[i] ^= 0x01;
        let err = verify(io::Cursor::new(&archive)).unwrap_err();
        assert!(matches!(err, Error::Corrupt { .. }), "got {err:?}");
    }

    #[test]
    fn decompress_detects_block_corruption() {
        let mut archive = multiblock_archive(40, 8);
        archive[HEADER_LEN + 8 + CRC_LEN] ^= 0xFF;
        let mut out = Vec::new();
        let err = decompress(&archive[..], &mut out, 1).unwrap_err();
        assert!(matches!(err, Error::Corrupt { .. }), "got {err:?}");
    }

    #[test]
    fn oversized_block_length_is_rejected_not_allocated() {
        let mut archive = multiblock_archive(40, 8);
        // Overwrite the first block's [8] length with a hostile value; the reader
        // must reject it up front instead of trying to allocate exabytes.
        archive[HEADER_LEN..HEADER_LEN + 8].copy_from_slice(&u64::MAX.to_le_bytes());
        let mut out = Vec::new();
        let err = decompress(&archive[..], &mut out, 1).unwrap_err();
        assert!(matches!(err, Error::Malformed(_)), "got {err:?}");
    }

    #[test]
    fn recover_skips_corrupt_block_and_keeps_the_rest() {
        let mut archive = multiblock_archive(40, 8);
        let footer = read_footer(&mut io::Cursor::new(&archive)).unwrap();
        assert!(
            footer.groups.len() >= 3,
            "need several blocks for this test"
        );
        let (off1, rc1) = footer.groups[1];
        // Corrupt one byte in block 1's payload (past its [8 len][4 crc]).
        archive[off1 as usize + 8 + CRC_LEN] ^= 0xFF;

        let mut out = Vec::new();
        let rec = decompress_recover(io::Cursor::new(&archive), &mut out, 1).unwrap();
        assert_eq!(rec.blocks_skipped, 1);
        assert_eq!(rec.reads_lost, u64::from(rc1));
        assert_eq!(rec.stats.reads, footer.total_reads - u64::from(rc1));
        assert_eq!(
            rec.blocks_recovered,
            footer.groups.len() as u64 - 1,
            "every other block recovered"
        );
        // Output is valid FASTQ for the recovered reads (4 lines each).
        assert_eq!(
            out.iter().filter(|&&b| b == b'\n').count() as u64,
            rec.stats.reads * 4
        );
    }

    #[test]
    fn truncated_at_block_boundary_streams_prefix() {
        let full = multiblock_archive(40, 8);
        // The trailer's back-pointer gives footer_offset; the 8-byte terminator
        // sits just before it, so footer_offset - 8 is a clean block boundary.
        let n = full.len();
        let footer_offset =
            u64::from_le_bytes(full[n - TRAILER_LEN..n - 4].try_into().unwrap()) as usize;
        let truncated = &full[..footer_offset - 8];

        // Streaming decode reads every whole block, then stops at the clean EOF.
        let mut out_trunc = Vec::new();
        decompress(truncated, &mut out_trunc, 1).expect("prefix decodes");
        let mut out_full = Vec::new();
        decompress(&full[..], &mut out_full, 1).unwrap();
        assert_eq!(
            out_trunc, out_full,
            "boundary-truncated file yields all blocks"
        );
    }

    #[test]
    fn reorder_archive_detects_frame_corruption() {
        let input = make_reads("y", 200);
        let mut archive = Vec::new();
        compress(
            &input[..],
            &mut archive,
            Params {
                reorder: true,
                keep_order: true,
                rescue: false,
                block_reads: 64,
                ..Params::default()
            },
        )
        .unwrap();
        // Flip a byte well past the header, inside a framed payload. The frame
        // CRC (or a downstream consistency check) must catch it — never a silent
        // wrong decode.
        let mid = archive.len() / 2;
        archive[mid] ^= 0xFF;
        assert!(verify(io::Cursor::new(&archive)).is_err());
    }
}
