//! The `.fqxv` container: a header followed by independent, parallel-codable
//! blocks.
//!
//! ```text
//! [4] magic "FQXV"
//! [2] format version (LE)
//! [1] sequence context order (k)
//! [1] quality binning tag
//! [1] flags (bit0: '+' line normalized)
//! repeated until EOF:
//!   [8] block payload length (LE)
//!   [ ] block payload
//! block payload:
//!   [4] n_reads (LE)
//!   [4] names_len (LE)  [ ] names   (fqxv-tokenizer)
//!   [4] seq_len   (LE)  [ ] seq     (fqxv-seq)
//!   [4] qual_len  (LE)  [ ] qual    (fqxv-fqzcomp)
//! ```

use std::io::{BufReader, BufWriter, Read, Write};

use rayon::prelude::*;

use crate::{Error, Result, FORMAT_VERSION, MAGIC};
use fqxv_fqzcomp::QualityBinning;

/// Reads per block. Blocks are the unit of parallelism and random access.
const BLOCK_READS: usize = 256 * 1024;
const FLAG_PLUS_NORMALIZED: u8 = 0x01;

/// Compression parameters.
#[derive(Debug, Clone, Copy)]
pub struct Params {
    /// Sequence context-model order (higher = better ratio, more memory).
    pub seq_order: u8,
    /// Quality quantization (lossless by default).
    pub quality_binning: QualityBinning,
    /// Worker threads (0 = all available cores).
    pub threads: usize,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            seq_order: 11,
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

/// Container header + per-stream size summary, from [`inspect`].
#[derive(Debug, Default, Clone)]
pub struct Info {
    /// Sequence context order.
    pub seq_order: u8,
    /// Quality binning tag (0 = lossless).
    pub quality_binning: u8,
    /// Whether the `+` line was normalized.
    pub plus_normalized: bool,
    /// Number of blocks.
    pub blocks: u64,
    /// Total reads.
    pub reads: u64,
    /// Compressed bytes in the name / sequence / quality streams.
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

/// Compress FASTQ from `reader` into a `.fqxv` stream on `writer`.
pub fn compress<R: Read, W: Write>(reader: R, writer: W, params: Params) -> Result<Stats> {
    let pool = build_pool(params.threads)?;
    let batch = pool.current_num_threads().max(1);

    let mut fq = noodles_fastq::io::Reader::new(BufReader::new(reader));
    let mut w = BufWriter::new(writer);

    // Header.
    w.write_all(&MAGIC)?;
    w.write_all(&FORMAT_VERSION.to_le_bytes())?;
    w.write_all(&[
        params.seq_order,
        binning_tag(params.quality_binning),
        FLAG_PLUS_NORMALIZED,
    ])?;

    let mut stats = Stats::default();
    let mut rec = noodles_fastq::Record::default();
    let mut eof = false;
    while !eof {
        // Fill up to `batch` blocks for this parallel round.
        let mut blocks: Vec<RawBlock> = Vec::with_capacity(batch);
        for _ in 0..batch {
            let mut b = RawBlock::default();
            while b.headers.len() < BLOCK_READS {
                if fq.read_record(&mut rec)? == 0 {
                    eof = true;
                    break;
                }
                let mut h = rec.name().to_vec();
                if !rec.description().is_empty() {
                    h.push(b' ');
                    h.extend_from_slice(rec.description());
                }
                b.headers.push(h);
                b.lens.push(rec.sequence().len() as u32);
                b.seq.extend_from_slice(rec.sequence());
                b.qual.extend_from_slice(rec.quality_scores());
            }
            if b.headers.is_empty() {
                break;
            }
            blocks.push(b);
            if eof {
                break;
            }
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
    stats.out_bytes += 9; // header bytes
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

/// Decompress a `.fqxv` stream from `reader` into FASTQ on `writer`.
pub fn decompress<R: Read, W: Write>(reader: R, writer: W, threads: usize) -> Result<Stats> {
    let pool = build_pool(threads)?;
    let batch = pool.current_num_threads().max(1);

    let mut r = BufReader::new(reader);
    let _header = read_header(&mut r)?;
    let mut w = BufWriter::new(writer);

    let mut stats = Stats::default();
    let mut eof = false;
    while !eof {
        let mut raw_blocks: Vec<Vec<u8>> = Vec::with_capacity(batch);
        for _ in 0..batch {
            match read_block(&mut r)? {
                Some(block) => raw_blocks.push(block),
                None => {
                    eof = true;
                    break;
                }
            }
        }
        if raw_blocks.is_empty() {
            break;
        }

        let decoded: Vec<Result<(u64, Vec<u8>)>> =
            pool.install(|| raw_blocks.par_iter().map(|b| decode_block(b)).collect());

        for d in decoded {
            let (reads, fastq) = d?;
            w.write_all(&fastq)?;
            stats.reads += reads;
            stats.blocks += 1;
            stats.out_bytes += fastq.len() as u64;
        }
    }
    w.flush()?;
    Ok(stats)
}

fn decode_block(buf: &[u8]) -> Result<(u64, Vec<u8>)> {
    let mut c = Cursor::new(buf);
    let n_reads = c.u32()? as usize;
    let names = fqxv_tokenizer::decode(c.slice_u32()?)?;
    let (seq_lens, seq) = fqxv_seq::decode(c.slice_u32()?)?;
    let (_qlens, qual) = fqxv_fqzcomp::decode(c.slice_u32()?)?;

    if names.len() != n_reads || seq_lens.len() != n_reads {
        return Err(Error::Malformed("block stream length disagreement"));
    }

    // Reassemble FASTQ: @name / seq / + / qual.
    let mut out = Vec::with_capacity(seq.len() * 2 + qual.len());
    let mut off = 0usize;
    for i in 0..n_reads {
        let l = seq_lens[i] as usize;
        out.push(b'@');
        out.extend_from_slice(&names[i]);
        out.push(b'\n');
        out.extend_from_slice(&seq[off..off + l]);
        out.extend_from_slice(b"\n+\n");
        out.extend_from_slice(&qual[off..off + l]);
        out.push(b'\n');
        off += l;
    }
    Ok((n_reads as u64, out))
}

/// Read the header and per-stream sizes without decoding block payloads.
pub fn inspect<R: Read>(reader: R) -> Result<Info> {
    let mut r = BufReader::new(reader);
    let header = read_header(&mut r)?;
    let mut info = Info {
        seq_order: header.seq_order,
        quality_binning: header.quality_binning,
        plus_normalized: header.flags & FLAG_PLUS_NORMALIZED != 0,
        ..Info::default()
    };
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
}

fn read_header<R: Read>(r: &mut R) -> Result<Header> {
    let mut magic = [0u8; 4];
    r.read_exact(&mut magic)?;
    if magic != MAGIC {
        return Err(Error::BadMagic);
    }
    let mut ver = [0u8; 2];
    r.read_exact(&mut ver)?;
    let ver = u16::from_le_bytes(ver);
    if ver != FORMAT_VERSION {
        return Err(Error::UnsupportedVersion(ver));
    }
    let mut p = [0u8; 3];
    r.read_exact(&mut p)?;
    Ok(Header {
        seq_order: p[0],
        quality_binning: p[1],
        flags: p[2],
    })
}

/// Read one length-prefixed block, or `None` at a clean EOF.
fn read_block<R: Read>(r: &mut R) -> Result<Option<Vec<u8>>> {
    let mut len = [0u8; 8];
    match r.read_exact(&mut len) {
        Ok(()) => {}
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => return Ok(None),
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
        .map_err(|e| Error::Io(std::io::Error::other(e.to_string())))
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

        // Name+description, sequence, quality preserved; '+' normalized to bare.
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

    #[test]
    fn reparse_is_stable() {
        // Compressing our own normalized output must reproduce it exactly.
        let a1 = compress_bytes(SAMPLE, Params::default());
        let mut f1 = Vec::new();
        decompress(&a1[..], &mut f1, 1).unwrap();
        let a2 = compress_bytes(&f1, Params::default());
        let mut f2 = Vec::new();
        decompress(&a2[..], &mut f2, 1).unwrap();
        assert_eq!(f1, f2);
    }

    #[test]
    fn inspect_reports_streams() {
        let archive = compress_bytes(SAMPLE, Params::default());
        let info = inspect(&archive[..]).expect("inspect");
        assert_eq!(info.reads, 2);
        assert_eq!(info.blocks, 1);
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
