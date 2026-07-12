//! The `.fqxv` container: a header followed by independent, parallel-codable
//! blocks.
//!
//! ```text
//! [4] magic "FQXV"
//! [2] format version (LE)
//! [1] sequence context order (k)
//! [1] quality binning tag
//! [1] flags (bit0: '+' normalized)
//! [1] group size G (reads interleaved per spot: 1 single-end, 2 paired,
//!                   3-4 single-cell R1/R2/I1[/I2], ...)
//! repeated until EOF:
//!   [8] block payload length (LE)
//!   [ ] block payload
//! block payload:
//!   [4] n_reads (LE)
//!   [4] names_len (LE)  [ ] names   (fqxv-tokenizer)
//!   [4] seq_len   (LE)  [ ] seq     (fqxv-seq)
//!   [4] qual_len  (LE)  [ ] qual    (fqxv-fqzcomp)
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

use crate::{Error, Result, FORMAT_VERSION, MAGIC};
use fqxv_fqzcomp::QualityBinning;

/// Default reads per block. Larger blocks populate the sequence model's contexts
/// better (higher ratio) but reduce parallelism and raise memory.
const DEFAULT_BLOCK_READS: usize = 1 << 20;
const HEADER_LEN: usize = 10;
const FLAG_PLUS_NORMALIZED: u8 = 0x01;

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
    /// Worker threads (0 = all available cores).
    pub threads: usize,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            seq_order: 11,
            block_reads: DEFAULT_BLOCK_READS,
            quality_binning: QualityBinning::Lossless,
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

#[derive(Default)]
struct RawBlock {
    headers: Vec<Vec<u8>>,
    lens: Vec<u32>,
    seq: Vec<u8>,
    qual: Vec<u8>,
}

impl RawBlock {
    fn push(&mut self, name: &[u8], description: &[u8], seq: &[u8], qual: &[u8]) {
        let mut h = name.to_vec();
        if !description.is_empty() {
            h.push(b' ');
            h.extend_from_slice(description);
        }
        self.headers.push(h);
        self.lens.push(seq.len() as u32);
        self.seq.extend_from_slice(seq);
        self.qual.extend_from_slice(qual);
    }
}

/// Compress single-end FASTQ from `reader` into a `.fqxv` stream.
pub fn compress<R: Read, W: Write>(reader: R, writer: W, params: Params) -> Result<Stats> {
    let block_reads = params.block_reads.max(1);
    let mut fq = noodles_fastq::io::Reader::new(BufReader::new(reader));
    let mut rec = noodles_fastq::Record::default();
    drive(writer, params, 1, |b| {
        while b.headers.len() < block_reads {
            if fq.read_record(&mut rec)? == 0 {
                break;
            }
            b.push(
                rec.name(),
                rec.description(),
                rec.sequence(),
                rec.quality_scores(),
            );
        }
        Ok(b.headers.len())
    })
}

/// Compress `G >= 1` per-spot read files (paired mates, single-cell R1/R2/I1/I2,
/// …) into one `.fqxv` stream, interleaving them.
///
/// Readers are consumed in lockstep; unequal read counts are an error. Restore
/// with [`decompress_split`], or stream interleaved with [`decompress`].
pub fn compress_multi<'a, W: Write>(
    readers: Vec<Box<dyn Read + 'a>>,
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

    // Keep whole spots together: round the block target down to a multiple of g.
    let block_reads = (params.block_reads / g).max(1) * g;
    drive(writer, params, g as u8, |b| {
        while b.headers.len() < block_reads {
            // Read one record from each input; member 0 EOF ends cleanly.
            let mut got = 0;
            for j in 0..g {
                if fqs[j].read_record(&mut recs[j])? == 0 {
                    if j == 0 {
                        return Ok(b.headers.len());
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
        Ok(b.headers.len())
    })
}

/// Shared block driver: `fill` populates one [`RawBlock`] and returns the number
/// of reads it added (0 at EOF). Blocks are compressed in parallel, written in
/// order.
fn drive<W, F>(writer: W, params: Params, group_size: u8, mut fill: F) -> Result<Stats>
where
    W: Write,
    F: FnMut(&mut RawBlock) -> Result<usize>,
{
    let pool = build_pool(params.threads)?;
    let batch = pool.current_num_threads().max(1);
    let mut w = BufWriter::new(writer);

    w.write_all(&MAGIC)?;
    w.write_all(&FORMAT_VERSION.to_le_bytes())?;
    w.write_all(&[
        params.seq_order,
        binning_tag(params.quality_binning),
        FLAG_PLUS_NORMALIZED,
        group_size,
    ])?;

    let mut stats = Stats::default();
    let mut done = false;
    while !done {
        let mut blocks: Vec<RawBlock> = Vec::with_capacity(batch);
        for _ in 0..batch {
            let mut b = RawBlock::default();
            if fill(&mut b)? == 0 {
                done = true;
                break;
            }
            blocks.push(b);
        }
        if blocks.is_empty() {
            break;
        }
        let compressed: Vec<Result<Vec<u8>>> = pool.install(|| {
            blocks
                .par_iter()
                .map(|b| compress_block(b, &params))
                .collect()
        });
        for (b, payload) in blocks.iter().zip(compressed) {
            let payload = payload?;
            w.write_all(&(payload.len() as u64).to_le_bytes())?;
            w.write_all(&payload)?;
            stats.reads += b.headers.len() as u64;
            stats.blocks += 1;
            stats.out_bytes += 8 + payload.len() as u64;
        }
    }
    w.flush()?;
    stats.out_bytes += HEADER_LEN as u64;
    Ok(stats)
}

fn compress_block(b: &RawBlock, params: &Params) -> Result<Vec<u8>> {
    let header_refs: Vec<&[u8]> = b.headers.iter().map(Vec::as_slice).collect();
    let names_c = fqxv_tokenizer::encode(&header_refs)?;
    let seq_c = fqxv_seq::encode(&b.lens, &b.seq, params.seq_order as usize)?;
    let qual_c = fqxv_fqzcomp::encode(&b.lens, &b.qual, params.quality_binning)?;

    let mut out = Vec::with_capacity(16 + names_c.len() + seq_c.len() + qual_c.len());
    out.extend_from_slice(&(b.headers.len() as u32).to_le_bytes());
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
pub fn decompress<R: Read, W: Write>(reader: R, writer: W, threads: usize) -> Result<Stats> {
    let pool = build_pool(threads)?;
    let batch = pool.current_num_threads().max(1);
    let mut r = BufReader::new(reader);
    let _ = read_header(&mut r)?;
    let mut w = BufWriter::new(writer);

    let mut stats = Stats::default();
    for_each_block_batch(&mut r, batch, |raw_blocks| {
        let decoded: Vec<Result<(u64, Vec<u8>)>> =
            pool.install(|| raw_blocks.par_iter().map(|b| decode_block(b)).collect());
        for d in decoded {
            let (reads, fastq) = d?;
            w.write_all(&fastq)?;
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

    let mut stats = Stats::default();
    for_each_block_batch(&mut r, batch, |raw_blocks| {
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
    let names = fqxv_tokenizer::decode(c.slice_u32()?)?;
    let (seq_lens, seq) = fqxv_seq::decode(c.slice_u32()?)?;
    let (_qlens, qual) = fqxv_fqzcomp::decode(c.slice_u32()?)?;
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
        write_record(&mut out, &names[i], &seq[off..off + l], &qual[off..off + l]);
        off += l;
    }
    Ok((n_reads as u64, out))
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
        ..Info::default()
    })
}

/// Read the header and per-stream sizes without decoding block payloads.
pub fn inspect<R: Read>(reader: R) -> Result<Info> {
    let mut r = BufReader::new(reader);
    let mut info = peek(&mut r)?;
    while let Some(block) = read_block(&mut r)? {
        let mut c = Cursor::new(&block);
        info.reads += c.u32()? as u64;
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
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads) // 0 => rayon default (all cores)
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
        let readers: Vec<Box<dyn Read>> =
            vec![Box::new(&r1[..]) as Box<dyn Read>, Box::new(&r2[..])];
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
        let readers: Vec<Box<dyn Read>> = files
            .iter()
            .map(|f| Box::new(&f[..]) as Box<dyn Read>)
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
        let readers: Vec<Box<dyn Read>> =
            vec![Box::new(&r1[..]) as Box<dyn Read>, Box::new(&r2[..])];
        compress_multi(readers, &mut archive, Params::default()).unwrap();
        let mut out = Vec::new();
        decompress(&archive[..], &mut out, 1).unwrap();
        assert_eq!(out, b"@r.1 a\nACGT\n+\nIIII\n@r.1 b\nGGGG\n+\n####\n");
    }

    #[test]
    fn unequal_mate_counts_error() {
        let r1 = make_reads("a", 2);
        let r2 = make_reads("b", 1);
        let readers: Vec<Box<dyn Read>> =
            vec![Box::new(&r1[..]) as Box<dyn Read>, Box::new(&r2[..])];
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
