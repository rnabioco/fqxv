//! Public compression entry points and the block-driving machinery.

use super::*;
use rayon::prelude::*;
use tracing::{debug, info, instrument, warn};

/// Default reads per block. Larger blocks populate the sequence model's contexts
/// better (higher ratio) but reduce parallelism and raise memory.
pub(crate) const DEFAULT_BLOCK_READS: usize = 1 << 20;

/// Per-block raw-sequence byte budget for Nanopore blocks when
/// [`Params::block_seq_bytes`] is 0 (auto). Long-read files hold few enough reads
/// that the read-count budget never binds, so this byte budget alone sets the
/// block count — and with the short-read 256 MiB cap a typical ONT file collapses
/// into ~2 blocks, leaving block-level decode parallelism nothing to work with
/// (issue #273).
///
/// 64 MiB is the measured knee of the ratio-vs-parallelism curve on a 576 MB ONT
/// MinION file (DRR205413, 21 k × ~14 kb reads): 2 → 5 blocks for **+1.08%**
/// archive size, full decode 56.4 s → 17.6 s at 8 threads (3.2×) and seq-only
/// (`--fasta`) 4.6 s → 1.9 s — where 256 MiB was flat at every thread count.
/// Compress drops 65 s → 25 s too (block-parallel encode was equally starved).
/// The cost is almost entirely the sequence stream (the per-block overlap codec
/// sees less coverage; quality pays only ~+0.5%), and it is what halving again
/// doubles: 32 MiB measured +2.7% for 5.3× at 16 threads. Block count scales
/// with file size at a fixed budget, so larger files gain parallelism while the
/// per-block ratio cost stays put.
pub(crate) const LONGREAD_BLOCK_SEQ_BYTES: usize = 64 << 20;

/// Resolve [`Params::block_seq_bytes`] against the resolved platform: 0 means
/// auto (the platform default — [`LONGREAD_BLOCK_SEQ_BYTES`] for Nanopore,
/// [`MAX_BLOCK_SEQ_BYTES`] otherwise); an explicit value is clamped to the hard
/// cap. Platform is a pure function of the input and the caller's `Params`, so
/// the budget — and with it the block boundaries — stays thread-count invariant.
pub(crate) fn resolve_block_seq_bytes(params: &Params, platform: Platform) -> usize {
    match params.block_seq_bytes {
        0 => match platform {
            Platform::Nanopore => LONGREAD_BLOCK_SEQ_BYTES,
            _ => MAX_BLOCK_SEQ_BYTES,
        },
        n => n.min(MAX_BLOCK_SEQ_BYTES),
    }
}

/// One interleaved spot's records — `(raw_header, sequence, quality)` per member —
/// owned so the platform can be detected before the streaming header is written
/// (see `compress_multi`). The header is the byte-exact definition line.
pub(crate) type PrimedSpot = Vec<(Vec<u8>, Vec<u8>, Vec<u8>)>;

/// Compression parameters.
#[derive(Debug, Clone)]
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
    /// Per-block raw-sequence byte budget; 0 (the default) means auto.
    ///
    /// Blocks are cut at whichever comes first — `block_reads` reads or this many
    /// raw sequence bytes. Long-read files hold few enough reads that the byte
    /// budget alone sets the block count, and blocks are the unit of decode
    /// parallelism and random access, so this is the granularity lever for
    /// long-read archives (issue #273). Auto resolves per platform: Nanopore gets
    /// [`LONGREAD_BLOCK_SEQ_BYTES`], everything else the 256 MiB cap. An explicit
    /// value is clamped to that cap; the CLI keeps `--block-reads`'s historical
    /// read-count-only semantics by pinning this to the cap when that flag is
    /// given. Ignored by the reorder path, which sizes its own blocks.
    pub block_seq_bytes: usize,
    /// Per-member slot labels recorded in the archive, one per interleaved member
    /// (`"R1"`, `"I1"`, `"2"`, …). Purely descriptive: nothing in decode depends on
    /// them, they exist so `decompress_split` can restore the *original* per-slot
    /// file names rather than renaming positionally. A run whose read slots are
    /// empty over part of the run yields files numbered by original slot with gaps
    /// (`_2` and `_4`, no `_1`/`_3`), and without these that identity is lost —
    /// positional naming turns `_2`/`_4` back into `_1`/`_2`.
    ///
    /// Empty (the default) writes no labels, which is exactly the 1.0 header. When
    /// non-empty the length must equal the group size, and the labels ride in a
    /// *non-critical* header extension record, so a reader that predates them skips
    /// it and decodes the archive normally.
    pub member_labels: Vec<String>,
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
            tile_band: 256,
            tile_max_refs: 1,
            block_seq_bytes: 0,
            member_labels: Vec::new(),
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

/// Widest interleaving [`detect_group_size`] will infer from read names: 4, the
/// single-cell R1/R2/I1/I2 layout. Wider spots need an explicit group size.
pub(crate) const MAX_AUTO_GROUP: usize = 4;

/// How many leading records [`compress_auto`] reads to decide whether a single
/// stream is interleaved per-spot data. Four spots at the widest layout we infer
/// is plenty to be confident while staying cheap for the common single-end case.
pub(crate) const AUTODETECT_PEEK: usize = 4 * MAX_AUTO_GROUP;

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

/// Split a read name into the base that identifies its *spot* and an optional
/// member marker. Handles the two conventions that actually mark a member: a
/// `/1`..`/4` name suffix, and the Casava mate field at the head of the
/// description (`@id 1:N:0:…` / `@id 2:N:0:…`).
///
/// The `:` is what makes a leading description digit a Casava mate field rather
/// than a numeric spot name. SRA deflines are `@RUN.5 5 length=150`, where that
/// second `5` is the spot name repeated on *every* member of the spot — reading
/// it as a member number made both mates of an interleaved SRA dump look like
/// member 5, and since [`detect_group_size`] wants members to be *distinct*,
/// equal markers vetoed the spot outright and the pairing was missed.
pub(crate) fn mate_key(rec: &noodles_fastq::Record) -> (&[u8], Option<u8>) {
    let name: &[u8] = rec.name().as_ref();
    if let [base @ .., b'/', m @ b'1'..=b'4'] = name {
        return (base, Some(*m));
    }
    let desc: &[u8] = rec.description().as_ref();
    let marker = match desc {
        [m @ b'1'..=b'4', b':', ..] => Some(*m),
        _ => None,
    };
    (name, marker)
}

/// True when a run of same-base records can be one spot: either no member is
/// marked (a bare repeated name, as interleaved SRA dumps emit) or every member
/// is marked and no two markers agree.
fn members_distinct(spot: &[(&[u8], Option<u8>)]) -> bool {
    if spot.iter().all(|(_, m)| m.is_none()) {
        return true;
    }
    let mut seen = 0u16;
    for (_, m) in spot {
        let Some(m) = m else { return false };
        let bit = 1u16 << (m - b'0');
        if seen & bit != 0 {
            return false;
        }
        seen |= bit;
    }
    true
}

/// Infer the interleaving of a single stream from its leading records.
///
/// Records of one spot share a name — bare-repeated, or plus a distinct member
/// marker — so the peek is cut into maximal runs of consecutive same-base names,
/// each run being one spot, and the layout is the common run length (up to
/// [`MAX_AUTO_GROUP`]). Anything ambiguous falls back to 1 (single-end), which is
/// always safe to archive.
///
/// Only *complete* runs count. The peek's final run is still open unless the
/// whole stream fit inside the peek (`saw_stream_end`), and counting a truncated
/// run would infer a narrower layout than the file really has — a 3-member
/// layout peeked 16 records deep ends on a 1-record fragment, which would
/// otherwise disagree with the five complete runs before it and fall back to 1.
///
/// The previous rule only ever compared records pairwise, so a 4-member spot
/// sharing one name — 10x R1/R2/I1/I2 as `sracha`/`fasterq-dump` emit it —
/// satisfied it at every pair and was archived as *paired*. `--split` then wrote
/// two files that each interleaved two real slots, at two different read lengths,
/// with no warning: the same shape rnabioco/sracha-rs#84 reports.
pub(crate) fn detect_group_size(peeked: &[noodles_fastq::Record], saw_stream_end: bool) -> u8 {
    let keys: Vec<(&[u8], Option<u8>)> = peeked.iter().map(mate_key).collect();
    let mut runs: Vec<usize> = Vec::new();
    let mut start = 0usize;
    for i in 1..=keys.len() {
        if i < keys.len() && keys[i].0 == keys[start].0 {
            continue;
        }
        if i < keys.len() || saw_stream_end {
            if !members_distinct(&keys[start..i]) {
                return 1;
            }
            runs.push(i - start);
        }
        start = i;
    }
    let g = match runs.split_first() {
        Some((&first, rest)) if rest.iter().all(|&r| r == first) => first,
        _ => return 1,
    };
    if (2..=MAX_AUTO_GROUP).contains(&g) {
        g as u8
    } else {
        1
    }
}

/// Compress a single FASTQ stream, auto-detecting its per-spot interleaving from
/// the leading read names (see `detect_group_size`) — paired, or single-cell up to
/// four members (R1/R2/I1/I2). This is what the CLI uses for a lone input, so
/// `sracha get -Z … | fqxv compress -` archives a download with the right spot
/// grouping and no flag. Detection only promotes on unambiguous spot names;
/// otherwise it behaves exactly like [`compress`]. `reorder` mode honours the
/// detected grouping too: grouped input is globally clustered and a permutation
/// restores the interleaving (see `encode_reordered`).
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
    // The peek's last record is genuinely the stream's last only if the reader
    // ran dry *and* the peek stopped short of its own cap — at the cap there may
    // be more records sitting in `prefix` beyond the ones we decoded.
    let g = detect_group_size(&peeked, eof && peeked.len() < AUTODETECT_PEEK);
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
    if g > 1 {
        check_spot_coherence(&prefix, g as usize);
    }
    if params.reorder {
        return encode_reordered(buffer_records(&prefix)?, writer, params, g);
    }
    compress_buffered(&prefix, writer, params, g)
}

/// The base name shared by every member of one spot: the definition line up to
/// its first separator, minus a `/1`..`/4` member suffix. The raw-bytes analogue
/// of [`mate_key`]'s first element, for the paths that never build a `noodles`
/// record.
pub(crate) fn spot_base(def: &[u8]) -> &[u8] {
    let end = def
        .iter()
        .position(|&b| b == b' ' || b == b'\t')
        .unwrap_or(def.len());
    match &def[..end] {
        [base @ .., b'/', b'1'..=b'4'] => base,
        name => name,
    }
}

/// Running tally of spots whose members disagree on their base name.
///
/// With `G > 1` the layout is purely positional — member `i` of every spot is
/// written to output `i` — so a stream that loses one member somewhere in the
/// middle shifts every later read into the wrong mate file. The read count stays
/// a clean multiple of `G`, so neither the spot-multiple check nor any CRC can
/// see it; only the names can. That is exactly what
/// `fasterq-dump --split-files` produces for a spot whose `READ_LEN` slot is 0 —
/// it writes no record at all for the empty slot (rnabioco/sracha-rs#76).
///
/// This warns rather than rejects: mates whose names genuinely disagree are
/// unusual but legal, and the archive is still byte-lossless either way. What is
/// wrong is only the spot grouping, which does not bite until `--split`.
#[derive(Default)]
pub(crate) struct SpotCoherence {
    spots: u64,
    bad: u64,
    first_bad: Option<u64>,
}

impl SpotCoherence {
    /// Record one spot, given its members' raw definition lines.
    fn observe<T: AsRef<[u8]>>(&mut self, defs: &[T]) {
        let base = spot_base(defs[0].as_ref());
        if defs[1..].iter().any(|d| spot_base(d.as_ref()) != base) {
            self.bad += 1;
            self.first_bad.get_or_insert(self.spots);
        }
        self.spots += 1;
    }

    /// Emit one warning covering the whole run, if any spot disagreed.
    fn finish(&self, g: usize) {
        if let Some(first) = self.first_bad {
            warn!(
                group_size = g,
                spots = self.spots,
                mismatched_spots = self.bad,
                first_mismatched_spot = first,
                "members of a spot do not share a read name — the inputs look \
                 out of step, so reads may be assigned to the wrong mate on \
                 --split (a converter that drops zero-length reads instead of \
                 emitting them does this)"
            );
        }
    }
}

/// [`SpotCoherence`] over an in-memory interleaved stream, without a second full
/// parse: FASTQ is strictly four lines per record, so every fourth line is a
/// definition. Bails out silently the moment that assumption does not hold —
/// the real parser runs next and reports malformed input far better than a
/// heuristic scan could.
fn check_spot_coherence(buf: &[u8], g: usize) {
    let mut cohere = SpotCoherence::default();
    let mut defs: Vec<&[u8]> = Vec::with_capacity(g);
    let mut line_no = 0usize;
    let mut pos = 0usize;
    while pos < buf.len() {
        let end = memchr::memchr(b'\n', &buf[pos..]).map_or(buf.len(), |k| pos + k);
        if line_no.is_multiple_of(4) {
            let mut def = &buf[pos..end];
            if def.last() == Some(&b'\r') {
                def = &def[..def.len() - 1];
            }
            match def.split_first() {
                Some((b'@', rest)) => defs.push(rest),
                _ => return,
            }
            if defs.len() == g {
                cohere.observe(&defs);
                defs.clear();
            }
        }
        line_no += 1;
        pos = end + 1;
    }
    cohere.finish(g);
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
    check_spot_coherence(&buf, g as usize);
    if params.reorder {
        return encode_reordered(buffer_records(&buf)?, writer, params, g);
    }
    compress_buffered(&buf, writer, params, g)
}

/// Confirm members 1..G are also spent once member 0 hits EOF. The interleaving
/// loops only ever touch a member after member 0, so a *short* member 0 ends the
/// archive early and drops the surplus reads from the longer files into a
/// well-formed, plausible-looking archive — the one unequal-count direction that
/// fails silently rather than loudly.
fn members_at_eof<R: BufRead>(fqs: &mut [R]) -> Result<()> {
    let (mut def, mut seq, mut qual) = (Vec::new(), Vec::new(), Vec::new());
    for fq in fqs.iter_mut().skip(1) {
        if read_raw_record(fq, &mut def, &mut seq, &mut qual)? {
            return Err(Error::Malformed("inputs have unequal read counts"));
        }
    }
    Ok(())
}

/// Compress `G >= 1` per-spot read files (paired mates, single-cell R1/R2/I1/I2,
/// …) into one `.fqxv` stream, interleaving them.
///
/// Readers are consumed in lockstep; unequal read counts are an error, whichever
/// member runs out first. Restore with [`decompress_split`], or stream
/// interleaved with [`decompress`].
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

    // Separate inputs can fall out of step too: a converter that drops a spot's
    // zero-length read from one file and a different spot's from another leaves
    // the counts equal but the spots misaligned (see `SpotCoherence`).
    let mut cohere = SpotCoherence::default();

    if params.reorder {
        // Buffer every spot in interleaved order (m0₀, m1₀, …), then globally
        // cluster; the stored permutation restores this spot order on decode, so
        // grouped reorder is order-preserving and de-interleaves cleanly.
        let mut all = RawBlock::default();
        loop {
            for j in 0..g {
                if !read_raw_record(&mut fqs[j], &mut defs[j], &mut seqs[j], &mut quals[j])? {
                    if j == 0 {
                        members_at_eof(&mut fqs)?;
                        cohere.finish(g);
                        return encode_reordered(all, writer, params, g as u8);
                    }
                    return Err(Error::Malformed("inputs have unequal read counts"));
                }
            }
            if g > 1 {
                cohere.observe(&defs);
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
                members_at_eof(&mut fqs)?;
                break; // empty input
            }
            return Err(Error::Malformed("inputs have unequal read counts"));
        }
        primed.push((defs[j].clone(), seqs[j].clone(), quals[j].clone()));
    }
    if g > 1 && primed.len() == g {
        cohere.observe(&defs);
    }
    let refs: Vec<&[u8]> = primed.iter().map(|(d, _, _)| d.as_slice()).collect();
    let platform = params.platform.unwrap_or_else(|| detect_platform(&refs));
    // Platform-resolved byte budget: must match `compress_buffered_plain`'s so an
    // explicit-Nanopore stream cuts the same blocks as the buffered path (#211's
    // byte-identical claim).
    let seq_budget = resolve_block_seq_bytes(&params, platform);
    let mut primed = Some(primed).filter(|p| !p.is_empty());
    let stats = drive(writer, params, g as u8, platform, |b| {
        // Emit the primed first spot into the first block before reading on.
        if let Some(spot) = primed.take() {
            for (header, seq, qual) in &spot {
                b.push_raw(header, seq, qual);
            }
        }
        // Cut on reads OR the raw-sequence byte budget, whichever comes first;
        // the loop reads whole spots, so a byte cut still lands on a spot
        // boundary. Matches the byte budgeting in `block_ranges`.
        while b.n_reads() < block_reads && b.seq.len() < seq_budget {
            // Read one record from each input; member 0 EOF ends cleanly.
            for j in 0..g {
                if !read_raw_record(&mut fqs[j], &mut defs[j], &mut seqs[j], &mut quals[j])? {
                    if j == 0 {
                        members_at_eof(&mut fqs)?;
                        return Ok(b.n_reads());
                    }
                    return Err(Error::Malformed("inputs have unequal read counts"));
                }
            }
            if g > 1 {
                cohere.observe(&defs);
            }
            for j in 0..g {
                b.push_raw(&defs[j], &seqs[j], &quals[j]);
            }
        }
        Ok(b.n_reads())
    })?;
    cohere.finish(g);
    Ok(stats)
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
    // cap, on whole-spot boundaries) — a pure function of the read lengths and the
    // resolved platform, so determinism holds regardless of thread count.
    let ranges = block_ranges(
        &chunks,
        block_reads,
        resolve_block_seq_bytes(&params, platform),
        g,
    );
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
    let header_len = write_header(&mut w, params, group_size, platform)?;
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
    let mut index = FooterIndex::new_at(header_len);
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
    stats.out_bytes += header_len + footer_bytes;
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
    // Auto keeps the full 256 MiB here (this path never runs for Nanopore), but an
    // explicit `block_seq_bytes` is honored: the shared reference is stored once
    // regardless of the cut, so granularity stays a per-block-stream decision.
    let ranges = block_ranges(
        &chunks,
        block_reads,
        resolve_block_seq_bytes(&params, platform),
        g,
    );
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
    let header_len = write_header_prefix(
        &mut w,
        &HeaderPrefix {
            seq_order: params.seq_order,
            binning: binning_tag(params.quality_binning),
            flags,
            group_size,
            platform,
            required_features: crate::feature::GLOBAL_REFERENCE,
        },
        &encode_member_labels(&params.member_labels, group_size),
    )?;
    write_framed(&mut w, &ref_frame)?;
    // Framed slice on disk is [4 len][4 crc][bytes]; blocks begin past it.
    let ref_frame_bytes = (4 + CRC_LEN + ref_frame.len()) as u64;

    let mut stats = Stats {
        group_size,
        ..Stats::default()
    };
    let mut index = FooterIndex::new_at(header_len + ref_frame_bytes);
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
    stats.out_bytes += header_len + ref_frame_bytes + footer_bytes;
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
    let header_len = write_header(&mut w, &params, group_size, platform)?;

    let mut stats = Stats {
        group_size,
        ..Stats::default()
    };
    let mut index = FooterIndex::new_at(header_len);

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
    stats.out_bytes += header_len + footer_bytes;
    Ok(stats)
}
