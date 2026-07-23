//! Compression-ratio estimation from a bounded input sample (`--estimate`).
//!
//! Predicts the archive size *without coding it*. The real codecs are entropy
//! coders, so they asymptotically reach the data's entropy — which a histogram pass
//! over the sample measures directly, at a fraction of the cost (no arithmetic
//! coding). Per stream:
//!
//! - **names** — coded exactly with the real tokenizer (a small stream, and its
//!   per-column delta structure defeats any raw byte-entropy proxy).
//! - **sequence** — the sample's static order-`k` empirical entropy on the
//!   short-read path (where the coder is pure order-k); on the long-read path,
//!   blended with a k-mer *duplication sketch* that captures the cross-read
//!   redundancy the overlap/tiler codecs exploit (which a per-read entropy is blind
//!   to). See [`long_read_seq_bits`].
//! - **quality** — the sample's static order-1 empirical entropy.
//!
//! Each measured stream is scaled by a small per-platform calibration factor fit
//! against real `compress` runs (the coders pay a learning premium over the ideal
//! static code length, and use richer context than the proxy). The passes run in
//! parallel across read chunks, so a whole sample is measured in well under a
//! second. The caller projects the full archive by scaling the sample by the
//! fraction of input it consumed.
//!
//! The reorder layout is deliberately NOT modeled (its cross-read redundancy grows
//! with read count, which a bounded sample understates), so the estimate is a
//! conservative lower bound for `--order any` — the real archive can only come out
//! equal or smaller.

use super::*;
use rayon::prelude::*;

/// Cap on records measured for a short-read sample. Entropy converges well before
/// this, and short reads take the pure order-k path where cross-read coverage is
/// irrelevant, so a modest count keeps parsing fast.
const ESTIMATE_MAX_READS: usize = 300_000;

/// Cap on sequence bases in any sample: one block's worth ([`MAX_BLOCK_SEQ_BYTES`]).
/// This is the cap that binds on long reads — it makes the sample's cross-read
/// coverage (which drives the duplication sketch) match a real block's.
const ESTIMATE_MAX_BASES: usize = MAX_BLOCK_SEQ_BYTES;

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
    /// Platform resolved from the sample (names, or content when SRA-renamed).
    /// `Unknown` for an empty input or names matching no convention. A caller
    /// estimating grouped inputs (paired mates) can compare this across the group
    /// to reject an accidental cross-platform mix before summing.
    pub platform: Platform,
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

/// Read leading records from `reader` into one [`RawBlock`], stopping at whichever
/// cap binds first — `max_reads` records or `max_bases` sequence bases — and
/// returning the block, the uncompressed `+`-normalized FASTQ bytes it represents,
/// and whether the input was exhausted (the whole file fit in the sample). The dual
/// cap lets [`estimate`] bound short-read samples by count (fast, coverage-
/// irrelevant) and long-read samples by bases (one block's worth, so the cross-read
/// redundancy the sample sees matches a real block). An empty input yields an empty
/// block.
fn sample_block<R: Read>(
    reader: R,
    max_reads: usize,
    max_bases: usize,
) -> Result<(RawBlock, u64, bool)> {
    let target = max_reads.max(1);
    let mut fq = noodles_fastq::io::Reader::new(BufReader::new(reader));
    let mut blk = RawBlock::default();
    let mut rec = noodles_fastq::Record::default();
    let mut raw_bytes = 0u64;
    let mut exhausted = false;
    while blk.n_reads() < target && blk.seq.len() < max_bases {
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
    Ok((blk, raw_bytes, exhausted))
}

/// Order of the sequence entropy proxy. Capped well below the coder's `seq_order`
/// (which reaches 11) so the context table stays `4^K` = 64k rows; the calibration
/// factor absorbs the gap between this static order-`K` code length and the coder's
/// higher adaptive order.
const SEQ_PROXY_ORDER: usize = 8;

/// Calibration factors mapping measured entropy to the bytes the real adaptive
/// coders emit. Fit against real `compress` runs on the sample; each currently rests
/// on one dataset per platform (NovaSeq / PacBio HiFi / Nanopore), so widen them
/// against the corpus before treating the numbers as tight.
///
/// - `CAL_SEQ_SHORT` — short-read sequence: the coder runs a higher adaptive order
///   (11) than the order-8 proxy, so the proxy runs ~16% high; this pulls it back.
/// - Quality is calibrated per platform — see [`cal_qual`].
/// - Long-read sequence is not a single factor — see [`long_read_seq_bits`].
const CAL_SEQ_SHORT: f64 = 0.86;

/// Calibration factor for the order-1 quality entropy, by resolved platform. The
/// real fqzcomp coder conditions HiFi quality on sequence and homopolymer-run
/// context (its Phase-1 model), so an order-1 proxy runs high on PacBio and needs
/// the largest correction; Illumina and Nanopore quality is close to order-1.
fn cal_qual(platform: Platform) -> f64 {
    match platform {
        Platform::PacBio => 0.83,
        Platform::Nanopore => 0.95,
        _ => 0.94,
    }
}

/// K-mer length for the cross-read redundancy sketch (matches the overlap codec's
/// anchor k).
const SKETCH_K: usize = 15;

/// Subsample rate for the sketch: keep ~1 in `2^SKETCH_SUBSAMPLE_BITS` k-mers. The
/// duplication rate is a ratio over the kept k-mers, so a coarse subsample bounds
/// the set size and cost without disturbing the estimate.
const SKETCH_SUBSAMPLE_BITS: u32 = 6;

/// Residual per-base cost (bits) charged to sequence covered by an earlier read on
/// the long-read path: the overlap codec still pays for each redundant region's edit
/// script (mismatches/indels) and anchoring, so a duplicated base is cheap but not
/// free. This tracks the platform's *error rate* — a noisier read carries a fatter
/// edit script per redundant base — so it is set per resolved platform. Fit so the
/// blended estimate reproduces the overlap codec's bits/base (PacBio/HiFi ≈ 0.07,
/// Nanopore ≈ 0.28; Unknown long reads take a middle value).
fn lr_residual_bits(platform: Platform) -> f64 {
    match platform {
        Platform::PacBio => 0.07,
        Platform::Nanopore => 0.28,
        _ => 0.20,
    }
}

/// Estimate the archive size for `reader`'s FASTQ from a bounded sample, measuring
/// each stream's entropy rather than coding it (see the module doc). `sample_reads`
/// caps the record count; a long-read sample is additionally capped at one block of
/// bases so its cross-read coverage matches a real block. `params` is honoured for
/// sequence order and quality binning; `params.reorder` is ignored (the estimate is
/// a lower bound for reorder). An empty input yields a zero estimate (0 reads,
/// `exhausted`), mirroring `compress`.
pub fn estimate<R: Read>(reader: R, params: Params, sample_reads: usize) -> Result<Estimate> {
    // Short reads: bound by count (entropy converges fast; coverage is irrelevant on
    // the pure order-k path). Long reads: bound by one block's worth of bases, so the
    // sketch sees the same cross-read coverage a real block would.
    let (blk, raw_bytes, exhausted) = sample_block(
        reader,
        sample_reads.min(ESTIMATE_MAX_READS),
        ESTIMATE_MAX_BASES,
    )?;
    if blk.n_reads() == 0 {
        return Ok(Estimate {
            sample_reads: 0,
            sample_bases: 0,
            raw_bytes: 0,
            names_bytes: 0,
            seq_bytes: 0,
            qual_bytes: 0,
            archive_bytes: 0,
            exhausted: true,
            platform: Platform::Unknown,
        });
    }

    // Resolve the platform (from names, or content when SRA-renamed) so the
    // long-read sequence blend can pick the right residual cost. Cheap — it peeks
    // only the leading records.
    let platform = resolve_platform_block(params.platform, &blk);

    // The three per-stream analyses are independent; run them concurrently. Names
    // are coded exactly with the real tokenizer (small stream; its delta-column
    // structure makes any raw-byte entropy proxy wildly pessimistic). Sequence and
    // quality are predicted from measured entropy.
    let k = SEQ_PROXY_ORDER.min(params.seq_order as usize).max(1);
    let (names_r, (seq_bits, qual_bits)) = rayon::join(
        || fqxv_tokenizer::encode(&blk.header_refs()).map(|c| c.len() as u64),
        || {
            rayon::join(
                // Sequence: mirror the real codec's own branch. Short reads take pure
                // order-k (the coder does too, so the sample's static entropy is a
                // faithful proxy); long reads take the overlap/tiler codecs, whose
                // cross-read redundancy a per-read entropy cannot see, so blend in a
                // k-mer duplication sketch.
                || {
                    if super::reorder::is_long_read(&blk.lens) {
                        long_read_seq_bits(&blk.lens, &blk.seq, k, platform)
                    } else {
                        seq_orderk_bits(&blk.lens, &blk.seq, k) * CAL_SEQ_SHORT
                    }
                },
                // Quality: static order-1 empirical entropy (previous quality in the
                // read), post-binning so a lossy estimate measures what is stored.
                || qual_order1_bits(&blk.lens, &blk.qual, params.quality_binning),
            )
        },
    );
    let names_bytes = names_r?;
    let seq_bytes = (seq_bits / 8.0).ceil() as u64;
    let qual_bytes = (qual_bits / 8.0 * cal_qual(platform)).ceil() as u64;

    // Per-block frame the real payload carries and that scales with block count:
    // stream digests + n_reads + three u32 stream-length prefixes + the block's
    // own [8 len][4 crc] frame.
    let frame = (STREAM_DIGESTS_LEN + 4 + 12 + 8 + CRC_LEN) as u64;
    let archive_bytes = names_bytes + seq_bytes + qual_bytes + frame;

    Ok(Estimate {
        sample_reads: blk.n_reads() as u64,
        sample_bases: blk.seq.len() as u64,
        raw_bytes,
        names_bytes,
        seq_bytes,
        qual_bytes,
        archive_bytes,
        exhausted,
        platform,
    })
}

/// Split `lens` into contiguous read chunks for parallel measurement, one per worker
/// thread. Each entry is `(read_lo, read_hi, byte_lo, byte_hi)` — the read index
/// range and its span in the concatenated sequence/quality buffer. Context is reset
/// per read, so chunks are fully independent.
fn read_chunks(lens: &[u32]) -> Vec<(usize, usize, usize, usize)> {
    let r = lens.len();
    if r == 0 {
        return Vec::new();
    }
    let per = r.div_ceil(rayon::current_num_threads().max(1));
    let mut chunks = Vec::new();
    let (mut read_lo, mut byte_lo) = (0usize, 0usize);
    while read_lo < r {
        let read_hi = (read_lo + per).min(r);
        let bytes: usize = lens[read_lo..read_hi].iter().map(|&l| l as usize).sum();
        chunks.push((read_lo, read_hi, byte_lo, byte_lo + bytes));
        byte_lo += bytes;
        read_lo = read_hi;
    }
    chunks
}

/// Total bits to code `seq` under a static order-`k` model over 2-bit ACGT symbols,
/// context reset at each read boundary (`lens`). Non-ACGT bases are billed a flat
/// literal cost (the coder routes them to an exception list) and do not enter the
/// context. This is the ideal per-context code length `Σ_ctx Σ_sym c·log2(T/c)`,
/// the entropy an order-`k` coder converges to. Chunks are counted in parallel and
/// their tables summed (context resets per read, so chunk boundaries are exact).
fn seq_orderk_bits(lens: &[u32], seq: &[u8], k: usize) -> f64 {
    // Flat cost charged per non-ACGT base: routed to the exception list, roughly a
    // byte of position/base overhead each. A crude constant — these are rare on the
    // platforms where the fast path matters.
    const NON_ACGT_BITS: f64 = 8.0;
    let n_ctx = 1usize << (2 * k);
    let mask = n_ctx - 1;
    let (counts, non_acgt) = read_chunks(lens)
        .par_iter()
        .map(|&(rlo, rhi, blo, bhi)| {
            // counts[ctx*4 + sym]; a flat Vec keeps the hot loop index-only.
            let mut counts = vec![0u32; n_ctx * 4];
            let mut non_acgt = 0u64;
            let mut off = blo;
            for &l in &lens[rlo..rhi] {
                let end = off + l as usize;
                let mut ctx = 0usize;
                for &b in &seq[off..end] {
                    let sym = fqxv_dna::code_fold(b);
                    if sym == fqxv_dna::NON_ACGT {
                        non_acgt += 1;
                        continue;
                    }
                    counts[ctx * 4 + sym as usize] += 1;
                    ctx = ((ctx << 2) | sym as usize) & mask;
                }
                off = end;
            }
            debug_assert_eq!(off, bhi);
            (counts, non_acgt)
        })
        .reduce(
            || (vec![0u32; n_ctx * 4], 0u64),
            |mut a, b| {
                for (x, y) in a.0.iter_mut().zip(&b.0) {
                    *x += *y;
                }
                a.1 += b.1;
                a
            },
        );
    context_code_bits(&counts, 4) + non_acgt as f64 * NON_ACGT_BITS
}

/// Predict the long-read sequence stream, where the overlap/tiler codecs exploit
/// cross-read redundancy that a per-read order-`k` entropy is blind to. Blend the
/// two regimes by the sample's k-mer duplication rate `d`: a base whose context is
/// novel costs the measured order-`k` entropy, one already covered by an earlier
/// read costs only the platform residual (the overlap codec's residual edit/anchor
/// cost — see [`lr_residual_bits`]). `bits = bases · ((1-d)·H_k + d·R)`.
fn long_read_seq_bits(lens: &[u32], seq: &[u8], k: usize, platform: Platform) -> f64 {
    let bases = seq.len() as f64;
    if bases == 0.0 {
        return 0.0;
    }
    let hk_per_base = seq_orderk_bits(lens, seq, k) / bases;
    let d = kmer_duplication_rate(lens, seq);
    let r = lr_residual_bits(platform);
    if std::env::var_os("FQXV_ESTIMATE_DEBUG").is_some() {
        eprintln!(
            "[est dbg] long-read seq: platform={platform:?} bases={} dup={:.4} H_k={:.4} R={:.2} -> blended={:.4} b/base",
            bases as u64,
            d,
            hk_per_base,
            r,
            (1.0 - d) * hk_per_base + d * r,
        );
    }
    bases * ((1.0 - d) * hk_per_base + d * r)
}

/// Fraction of the sample's k-mers that repeat an earlier k-mer — a cheap proxy for
/// cross-read redundancy (coverage). A [`SKETCH_K`]-mer is rolled per read (reset at
/// read boundaries and non-ACGT bases), Fibonacci-hashed, and subsampled to ~1 in
/// `2^SKETCH_SUBSAMPLE_BITS` by its top bits. The duplicate fraction is
/// `1 - distinct/total` over the kept hashes — no ordering needed, so chunks emit
/// their kept hashes in parallel and the distinct count is a single sort+dedup.
/// `0.0` when the sample has no k-mers.
fn kmer_duplication_rate(lens: &[u32], seq: &[u8]) -> f64 {
    // Fibonacci hashing: multiply by 2^64/φ and read the well-mixed high bits.
    const PHI64: u64 = 0x9E37_79B9_7F4A_7C15;
    let kmask: u64 = if 2 * SKETCH_K >= 64 {
        u64::MAX
    } else {
        (1u64 << (2 * SKETCH_K)) - 1
    };
    let mut kept: Vec<u64> = read_chunks(lens)
        .par_iter()
        .map(|&(rlo, rhi, blo, _bhi)| {
            let mut local = Vec::new();
            let mut off = blo;
            for &l in &lens[rlo..rhi] {
                let end = off + l as usize;
                let mut kmer = 0u64;
                let mut valid = 0usize;
                for &b in &seq[off..end] {
                    let sym = fqxv_dna::code_fold(b);
                    if sym == fqxv_dna::NON_ACGT {
                        valid = 0;
                        kmer = 0;
                        continue;
                    }
                    kmer = ((kmer << 2) | u64::from(sym)) & kmask;
                    valid += 1;
                    if valid < SKETCH_K {
                        continue;
                    }
                    let h = kmer.wrapping_mul(PHI64);
                    if h >> (64 - SKETCH_SUBSAMPLE_BITS) == 0 {
                        local.push(h);
                    }
                }
                off = end;
            }
            local
        })
        .reduce(Vec::new, |mut a, mut b| {
            a.append(&mut b);
            a
        });
    let total = kept.len();
    if total == 0 {
        return 0.0;
    }
    kept.sort_unstable();
    kept.dedup();
    1.0 - kept.len() as f64 / total as f64
}

/// Total bits to code `qual` under a static order-1 model (previous quality byte in
/// the read as context), context reset at each read boundary (`lens`) and each byte
/// mapped through `binning` first so a lossy estimate measures the stored alphabet.
/// Chunks are counted in parallel and their tables summed.
fn qual_order1_bits(lens: &[u32], qual: &[u8], binning: QualityBinning) -> f64 {
    // Quality bytes span the printable Phred range; a full 256-wide context keeps
    // the mapping exact without assuming an alphabet.
    const A: usize = 256;
    let counts = read_chunks(lens)
        .par_iter()
        .map(|&(rlo, rhi, blo, _bhi)| {
            let mut counts = vec![0u32; A * A];
            let mut off = blo;
            for &l in &lens[rlo..rhi] {
                let end = off + l as usize;
                let mut prev = 0usize; // read-start context
                for &q in &qual[off..end] {
                    let cur = binning.apply(q) as usize;
                    counts[prev * A + cur] += 1;
                    prev = cur;
                }
                off = end;
            }
            counts
        })
        .reduce(
            || vec![0u32; A * A],
            |mut a, b| {
                for (x, y) in a.iter_mut().zip(&b) {
                    *x += *y;
                }
                a
            },
        );
    context_code_bits(&counts, A)
}

/// Ideal static code length in bits for a context→symbol count table laid out as
/// `row_width` consecutive symbol counts per context: `Σ_ctx Σ_sym c·log2(T_ctx/c)`,
/// where `T_ctx` is the row total. Empty rows contribute nothing.
fn context_code_bits(counts: &[u32], row_width: usize) -> f64 {
    let mut bits = 0.0f64;
    for row in counts.chunks_exact(row_width) {
        let total: u64 = row.iter().map(|&c| u64::from(c)).sum();
        if total == 0 {
            continue;
        }
        let lt = (total as f64).log2();
        for &c in row {
            if c != 0 {
                bits += f64::from(c) * (lt - (f64::from(c)).log2());
            }
        }
    }
    bits
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

    /// Long reads (mean length > the long-read threshold): `n` records of `len`
    /// high-entropy (LCG-random) bases. `distinct` seeds each read differently (low
    /// cross-read redundancy); otherwise every read is identical (maximal
    /// redundancy). High per-base entropy is the point — cyclic ACGT would collapse
    /// to nearly zero order-k entropy and hide the sketch's effect.
    fn long_fastq(n: usize, len: usize, distinct: bool) -> Vec<u8> {
        let bases = b"ACGT";
        let mut buf = Vec::new();
        for i in 0..n {
            let mut state: u64 = if distinct { i as u64 + 1 } else { 1 };
            let seq: Vec<u8> = (0..len)
                .map(|_| {
                    state = state
                        .wrapping_mul(6_364_136_223_846_793_005)
                        .wrapping_add(1_442_695_040_888_963_407);
                    bases[((state >> 33) & 3) as usize]
                })
                .collect();
            buf.extend_from_slice(format!("@read{i}\n").as_bytes());
            buf.extend_from_slice(&seq);
            buf.extend_from_slice(b"\n+\n");
            buf.extend_from_slice(&vec![b'5'; len]);
            buf.push(b'\n');
        }
        buf
    }

    #[test]
    fn estimate_is_deterministic() {
        // The measurement passes run in parallel (rayon); the result must not depend
        // on chunking or reduction order.
        let data = fastq(5000);
        let a = estimate(&data[..], Params::default(), 100_000).unwrap();
        let b = estimate(&data[..], Params::default(), 100_000).unwrap();
        assert_eq!(a.archive_bytes, b.archive_bytes);
        assert_eq!(a.seq_bytes, b.seq_bytes);
        assert_eq!(a.qual_bytes, b.qual_bytes);
    }

    #[test]
    fn long_read_duplication_shrinks_sequence() {
        // The k-mer sketch must credit cross-read redundancy: identical long reads
        // code their sequence far smaller than all-distinct long reads of the same
        // size, even though the per-read order-k entropy is similar.
        let redundant = estimate(
            &long_fastq(200, 1000, false)[..],
            Params::default(),
            100_000,
        )
        .unwrap();
        let distinct =
            estimate(&long_fastq(200, 1000, true)[..], Params::default(), 100_000).unwrap();
        assert!(
            redundant.seq_bytes * 2 < distinct.seq_bytes,
            "redundant seq {} should be well below distinct seq {}",
            redundant.seq_bytes,
            distinct.seq_bytes
        );
    }

    #[test]
    fn estimate_empty_input_is_zero() {
        // Empty input is valid (compresses to an empty archive), so it yields a
        // zero estimate rather than an error.
        let est = estimate(&b""[..], Params::default(), 1000).expect("empty is Ok");
        assert_eq!(est.sample_reads, 0);
        assert_eq!(est.sample_bases, 0);
        assert_eq!(est.raw_bytes, 0);
        assert_eq!(est.archive_bytes, 0);
        assert!(est.exhausted);
        assert_eq!(est.ratio(), 0.0);
    }
}
