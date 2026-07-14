//! Compression-ratio estimation from a bounded input sample (`--estimate`).
//!
//! Predicts the archive size without a full compress by sampling the leading
//! records, coding them with the real per-block codecs (names → tokenizer,
//! sequence → order-k, quality → fqzcomp), and reporting the sample's per-stream
//! sizes. Blocks are coded independently and the context models are per-block
//! stationary, so the sample's compression ratio is a faithful proxy for the
//! whole file — the caller projects the full archive by scaling the sample by the
//! fraction of input it consumed.
//!
//! The reorder layout is deliberately NOT modeled: its cross-read redundancy
//! grows with read count, so a small sample would understate it. The estimate
//! always codes the sample with the non-reorder path, which makes it a
//! conservative lower bound for `--order any` (the real archive can only come out
//! equal or smaller).

use super::*;

/// Per-stream sizes from a bounded compression sample (see [`estimate`]).
#[derive(Debug, Clone, Copy)]
pub struct Estimate {
    /// Records coded in the sample.
    pub sample_reads: u64,
    /// Sequence bases in the sample (sum of read lengths).
    pub sample_bases: u64,
    /// Uncompressed, `+`-normalized FASTQ bytes the sample represents — the exact
    /// bytes `decompress` reconstructs for these reads (`@name\nseq\n+\nqual\n`).
    pub raw_bytes: u64,
    /// Compressed names stream (tokenizer).
    pub names_bytes: u64,
    /// Compressed sequence stream (order-k).
    pub seq_bytes: u64,
    /// Compressed quality stream (fqzcomp).
    pub qual_bytes: u64,
    /// Full sample archive bytes: the three streams plus the per-block payload
    /// header and frame (digest, counts, length prefixes, CRC) — everything that
    /// scales with the data. The fixed file header/footer is excluded (a constant
    /// few tens of bytes, amortized to nothing at scale).
    pub archive_bytes: u64,
    /// True when the whole input fit inside the sample, so these numbers are the
    /// actual full-file compression rather than an extrapolation base.
    pub exhausted: bool,
}

impl Estimate {
    /// Compression ratio the archive achieves on this data: uncompressed FASTQ
    /// divided by the coded archive. Scale-invariant, so it holds for the whole
    /// file even though only a sample was coded. `0.0` for an empty archive.
    pub fn ratio(&self) -> f64 {
        if self.archive_bytes == 0 {
            0.0
        } else {
            self.raw_bytes as f64 / self.archive_bytes as f64
        }
    }
}

/// Code a bounded sample of `reader`'s leading records with the real per-block
/// codecs and report the per-stream sizes, for `--estimate`.
///
/// Reads at most `sample_reads` records into one block and codes it. The
/// `reorder` flag in `params` is ignored — the sample is always coded with the
/// non-reorder path (see the module doc); every other parameter (sequence order,
/// hashed tier, quality binning) is honoured so the estimate tracks the settings
/// the real run would use. Errors if the input has no reads.
pub fn estimate<R: Read>(reader: R, params: Params, sample_reads: usize) -> Result<Estimate> {
    let target = sample_reads.max(1);
    let mut fq = noodles_fastq::io::Reader::new(BufReader::new(reader));
    let mut blk = RawBlock::default();
    let mut rec = noodles_fastq::Record::default();
    let mut raw_bytes = 0u64;
    let mut exhausted = false;
    while blk.n_reads() < target {
        if fq.read_record(&mut rec)? == 0 {
            exhausted = true;
            break;
        }
        let (name, desc) = (rec.name(), rec.description());
        let (seq, qual) = (rec.sequence(), rec.quality_scores());
        // Normalized FASTQ record size = the bytes `decompress` emits for it:
        // `@` + header + `\n` + seq + `\n+\n` + qual + `\n` (6 fixed bytes). The
        // header is `name`, plus a single space and the description when present.
        let header = name.len() + if desc.is_empty() { 0 } else { 1 + desc.len() };
        raw_bytes += (header + seq.len() + qual.len() + 6) as u64;
        blk.push(name, desc, seq, qual);
    }
    if blk.n_reads() == 0 {
        return Err(Error::Malformed("input has no reads to estimate from"));
    }

    // Always code the sample with the non-reorder path: it is the accurate,
    // scale-invariant baseline, and a lower bound for reorder (module doc).
    let mut p = params;
    p.reorder = false;
    let payload = compress_block(&blk, &p)?;

    // Recover the three stream sizes from the payload framing
    // (`[8 digest][4 n_reads][ (u32 len + bytes) × 3 ]`).
    let mut c = Cursor::new(&payload[..]);
    c.u64()?; // content digest
    c.u32()?; // n_reads
    let names_bytes = c.slice_u32()?.len() as u64;
    let seq_bytes = c.slice_u32()?.len() as u64;
    let qual_bytes = c.slice_u32()?.len() as u64;

    // On-disk cost of this block: the payload plus its frame (`[8 len][4 crc]`).
    let archive_bytes = payload.len() as u64 + 8 + CRC_LEN as u64;

    Ok(Estimate {
        sample_reads: blk.n_reads() as u64,
        sample_bases: blk.seq.len() as u64,
        raw_bytes,
        names_bytes,
        seq_bytes,
        qual_bytes,
        archive_bytes,
        exhausted,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn fastq(n: usize) -> Vec<u8> {
        // Deterministic, mildly compressible reads (fixed length, cycling bases).
        let mut buf = Vec::new();
        let bases = b"ACGT";
        for i in 0..n {
            let seq: Vec<u8> = (0..100).map(|j| bases[(i + j) % 4]).collect();
            let qual = vec![b'I'; 100];
            buf.extend_from_slice(format!("@read{i} extra\n").as_bytes());
            buf.extend_from_slice(&seq);
            buf.extend_from_slice(b"\n+\n");
            buf.extend_from_slice(&qual);
            buf.push(b'\n');
        }
        buf
    }

    #[test]
    fn estimate_reports_whole_small_input() {
        let data = fastq(1000);
        let est = estimate(&data[..], Params::default(), 100_000).unwrap();
        assert_eq!(est.sample_reads, 1000);
        assert_eq!(est.sample_bases, 100_000);
        assert!(
            est.exhausted,
            "sample cap exceeds input, so it is exhausted"
        );
        // Raw = per-record `@read{i} extra\n` header + 100 seq + `\n+\n` + 100 qual
        // + trailing `\n`. The header text is `read{i} extra`; the record adds 6
        // fixed bytes on top (see `write_record`).
        let expected: u64 = (0..1000)
            .map(|i| (format!("read{i} extra").len() + 100 + 100 + 6) as u64)
            .sum();
        assert_eq!(est.raw_bytes, expected);
        // The archive is smaller than the raw input and no smaller than its coded
        // streams.
        assert!(est.archive_bytes < est.raw_bytes);
        assert!(est.archive_bytes >= est.names_bytes + est.seq_bytes + est.qual_bytes);
        assert!(est.ratio() > 1.0);
    }

    #[test]
    fn estimate_caps_the_sample() {
        let est = estimate(&fastq(1000)[..], Params::default(), 250).unwrap();
        assert_eq!(est.sample_reads, 250);
        assert!(!est.exhausted, "input has more reads than the cap");
    }

    #[test]
    fn estimate_ignores_reorder_flag() {
        // Reorder is never modeled: the flag must not change the coded sample.
        let data = fastq(500);
        let plain = estimate(&data[..], Params::default(), 100_000).unwrap();
        let reorder = estimate(
            &data[..],
            Params {
                reorder: true,
                ..Params::default()
            },
            100_000,
        )
        .unwrap();
        assert_eq!(plain.archive_bytes, reorder.archive_bytes);
        assert_eq!(plain.seq_bytes, reorder.seq_bytes);
    }

    #[test]
    fn estimate_rejects_empty_input() {
        assert!(estimate(&b""[..], Params::default(), 1000).is_err());
    }
}
