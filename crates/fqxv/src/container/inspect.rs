//! Container inspection, platform detection, and metadata reporting.

use super::*;

/// Leading reads sampled to guess the platform (see [`detect_platform`]).
pub(crate) const PLATFORM_PEEK: usize = 16;

/// Sequencing platform recorded in the archive metadata. Detected from read-name
/// grammar at compress time (see [`detect_platform`]) and overridable via
/// [`Params::platform`]; stored once in the container header, since it is a
/// per-archive fact, not a per-read one.
#[derive(Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum Platform {
    /// Not recorded, or the read names matched no known convention.
    #[default]
    Unknown,
    /// Illumina — `instrument:run:flowcell:lane:tile:x:y` colon-delimited names.
    Illumina,
    /// Oxford Nanopore — UUID read names, `runid=`/`ch=` description tags.
    Nanopore,
    /// PacBio — `movie/zmw[/ccs]` read names.
    PacBio,
    /// MGI / BGI — `V`/`E`/`DP`-prefixed `…L<lane>C…R…` read names.
    MgiBgi,
}

impl Platform {
    /// Numeric tag stored in the header's dedicated platform byte.
    pub(crate) fn to_code(self) -> u8 {
        match self {
            Platform::Unknown => 0,
            Platform::Illumina => 1,
            Platform::Nanopore => 2,
            Platform::PacBio => 3,
            Platform::MgiBgi => 4,
        }
    }

    pub(crate) fn from_code(code: u8) -> Self {
        match code {
            1 => Platform::Illumina,
            2 => Platform::Nanopore,
            3 => Platform::PacBio,
            4 => Platform::MgiBgi,
            _ => Platform::Unknown,
        }
    }

    /// Human-facing label for the `info` report.
    pub fn label(self) -> &'static str {
        match self {
            Platform::Unknown => "unknown",
            Platform::Illumina => "Illumina",
            Platform::Nanopore => "Oxford Nanopore",
            Platform::PacBio => "PacBio",
            Platform::MgiBgi => "MGI/BGI",
        }
    }

    /// Stable lowercase token for TSV/JSON output and the `--platform` flag.
    pub fn token(self) -> &'static str {
        match self {
            Platform::Unknown => "unknown",
            Platform::Illumina => "illumina",
            Platform::Nanopore => "nanopore",
            Platform::PacBio => "pacbio",
            Platform::MgiBgi => "mgi",
        }
    }
}

/// Guess the sequencing platform from a sample of read headers (each the name
/// plus any description, as stored in a [`RawBlock`]). Read-name grammar is a
/// dead giveaway per platform; a per-header vote is taken and the most common
/// non-`Unknown` verdict wins. Returns [`Platform::Unknown`] when nothing
/// matches, so a wrong platform is never recorded.
pub(crate) fn detect_platform(headers: &[&[u8]]) -> Platform {
    let mut votes = [0u32; 5];
    for &h in headers {
        votes[classify_header(h).to_code() as usize] += 1;
    }
    // Most-voted platform among the real ones; ties resolve to the lower code.
    let mut best = Platform::Unknown;
    let mut best_votes = 0;
    for code in 1..5u8 {
        if votes[code as usize] > best_votes {
            best_votes = votes[code as usize];
            best = Platform::from_code(code);
        }
    }
    best
}

/// Classify one read header by its name grammar. The header is the first
/// whitespace-delimited token (the name) plus an optional description tail; both
/// carry platform signal (Illumina packs everything in the name; Nanopore's
/// `runid=`/`ch=` tags live in the description).
pub(crate) fn classify_header(header: &[u8]) -> Platform {
    let (name, desc) = match header.iter().position(|&b| b == b' ') {
        Some(i) => (&header[..i], &header[i + 1..]),
        None => (header, &[][..]),
    };
    if is_uuid(name)
        || [b"runid=".as_slice(), b"read=", b"ch=", b"start_time="]
            .iter()
            .any(|&needle| contains_sub(desc, needle))
    {
        return Platform::Nanopore;
    }
    if is_pacbio_name(name) {
        return Platform::PacBio;
    }
    if is_mgi_name(name) {
        return Platform::MgiBgi;
    }
    if is_illumina_name(name) {
        return Platform::Illumina;
    }
    Platform::Unknown
}

/// Strip a trailing `/1` or `/2` mate marker from a read name.
pub(crate) fn strip_mate(name: &[u8]) -> &[u8] {
    match name {
        [base @ .., b'/', b'1' | b'2'] => base,
        _ => name,
    }
}

/// True if `needle` occurs in `hay` (small-slice substring search).
pub(crate) fn contains_sub(hay: &[u8], needle: &[u8]) -> bool {
    !needle.is_empty() && hay.windows(needle.len()).any(|w| w == needle)
}

/// True if `s` is a canonical 8-4-4-4-12 hex UUID (a Nanopore read id).
pub(crate) fn is_uuid(s: &[u8]) -> bool {
    if s.len() != 36 {
        return false;
    }
    s.iter().enumerate().all(|(i, &b)| {
        if matches!(i, 8 | 13 | 18 | 23) {
            b == b'-'
        } else {
            b.is_ascii_hexdigit()
        }
    })
}

/// True for PacBio `movie/zmw[/…]` names: a movie id (`m…` with digits), a numeric
/// ZMW hole number, and at least one more `/`-part (subread coords or `ccs`).
pub(crate) fn is_pacbio_name(name: &[u8]) -> bool {
    let mut parts = name.split(|&b| b == b'/');
    let movie = parts.next().unwrap_or_default();
    let zmw = match parts.next() {
        Some(z) => z,
        None => return false,
    };
    if parts.next().is_none() {
        return false;
    }
    let movie_ok = movie.first() == Some(&b'm') && movie.iter().any(u8::is_ascii_digit);
    let zmw_ok = !zmw.is_empty() && zmw.iter().all(u8::is_ascii_digit);
    movie_ok && zmw_ok
}

/// True for MGI/BGI names: a 1-3 letter prefix (`V`, `E`, `CL`, `DP`), a digit
/// flowcell id, then the `L<lane>…C…R…` tile layout.
pub(crate) fn is_mgi_name(name: &[u8]) -> bool {
    let base = strip_mate(name);
    let letters = base.iter().take_while(|b| b.is_ascii_uppercase()).count();
    if !(1..=3).contains(&letters) {
        return false;
    }
    let after_letters = &base[letters..];
    let digits = after_letters
        .iter()
        .take_while(|b| b.is_ascii_digit())
        .count();
    if digits == 0 {
        return false;
    }
    let rest = &after_letters[digits..];
    // `L<digit>` lane marker, with column/row markers following.
    matches!(rest, [b'L', d, ..] if d.is_ascii_digit())
        && rest.contains(&b'C')
        && rest.contains(&b'R')
}

/// True for Illumina names: 5+ colon-delimited fields (7 in Casava 1.8+, 5 in
/// older machines) ending in numeric x/y coordinates.
pub(crate) fn is_illumina_name(name: &[u8]) -> bool {
    let base = strip_mate(name);
    // Older `…#index/mate` names carry an index after a `#`; drop it.
    let base = match base.iter().position(|&b| b == b'#') {
        Some(i) => &base[..i],
        None => base,
    };
    let fields: Vec<&[u8]> = base.split(|&b| b == b':').collect();
    if fields.len() < 5 {
        return false;
    }
    let numeric = |f: &[u8]| !f.is_empty() && f.iter().all(u8::is_ascii_digit);
    numeric(fields[fields.len() - 1]) && numeric(fields[fields.len() - 2])
}

/// Resolve the platform to record: the caller's override if set, else a guess
/// from the leading read headers of an in-memory FASTQ buffer.
pub(crate) fn resolve_platform_buf(forced: Option<Platform>, buf: &[u8]) -> Platform {
    if let Some(p) = forced {
        return p;
    }
    let mut fq = noodles_fastq::io::Reader::new(buf);
    let mut rec = noodles_fastq::Record::default();
    let mut headers: Vec<Vec<u8>> = Vec::with_capacity(PLATFORM_PEEK);
    for _ in 0..PLATFORM_PEEK {
        match fq.read_record(&mut rec) {
            Ok(0) | Err(_) => break,
            Ok(_) => headers.push(join_header(rec.name(), rec.description())),
        }
    }
    let refs: Vec<&[u8]> = headers.iter().map(Vec::as_slice).collect();
    detect_platform(&refs)
}

/// Resolve the platform from the leading headers of an already-buffered block
/// (the reorder path buffers every read before writing the header).
pub(crate) fn resolve_platform_block(forced: Option<Platform>, all: &RawBlock) -> Platform {
    if let Some(p) = forced {
        return p;
    }
    let n = all.n_reads().min(PLATFORM_PEEK);
    let refs: Vec<&[u8]> = (0..n).map(|i| all.header(i)).collect();
    detect_platform(&refs)
}

/// Join a read `name` and optional `description` the way a [`RawBlock`] stores a
/// header, so platform detection sees the same bytes on every code path.
pub(crate) fn join_header(name: &[u8], description: &[u8]) -> Vec<u8> {
    let mut h = name.to_vec();
    if !description.is_empty() {
        h.push(b' ');
        h.extend_from_slice(description);
    }
    h
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
    /// Whether names were regenerated from a counter template (discard-order,
    /// reorder-lossy — reads were renumbered) rather than coded per read.
    pub regenerated_names: bool,
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
    /// Sequencing platform recorded at compress time (from read-name grammar or
    /// an explicit override); [`Platform::Unknown`] if none was recorded.
    pub platform: Platform,
    /// On-disk container format version. Always equals [`crate::FORMAT_VERSION`]
    /// for a readable archive (`read_header` rejects any other), but surfaced so
    /// tooling can report it without re-parsing the header.
    pub format_version: u16,
    /// The archive's stored whole-file CRC-32C (footer field), a stable
    /// fingerprint of the on-disk bytes through the `total_reads` field. `None`
    /// for the footer-less whole-file-reorder layout and for truncated archives
    /// whose footer could not be read (metadata then comes from a forward scan).
    /// This is the value `verify` recomputes and checks; reporting it lets a user
    /// record the expected checksum without a full pass.
    pub whole_file_crc: Option<u32>,
}

/// Highest Phred quality value tracked in [`ContentStats::qual_hist`]. Raw
/// quality bytes are printable ASCII (33..=126), so the Phred value (`byte - 33`)
/// tops out at 93; values are clamped into `0..QUAL_MAX` defensively.
pub const QUAL_MAX: usize = 94;

/// Content-level statistics over an archive's *decoded* reads — the data the
/// `--stats` pass reports (read-length spread, base composition, quality
/// distribution). Distinct from [`Info`], which is container metadata read from
/// the header and footer without decoding; computing these requires a full
/// decode, so [`content_stats`] runs the normal decompressor and folds its
/// output rather than re-implementing any codec.
///
/// Base counts are over the stored (post-binning for quality; sequence is never
/// lossy) content and are order-independent, so they match regardless of any
/// reordering. `a`/`c`/`g`/`t`/`n` count uppercase ACGT/N; every other byte
/// (IUPAC ambiguity codes, lowercase) falls in `other`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ContentStats {
    /// Number of reads decoded.
    pub reads: u64,
    /// Total sequence bases across all reads.
    pub bases: u64,
    /// Shortest / longest read length seen (both 0 when there are no reads).
    pub min_len: u32,
    /// Longest read length seen (0 when there are no reads).
    pub max_len: u32,
    /// Uppercase `A` base count.
    pub a: u64,
    /// Uppercase `C` base count.
    pub c: u64,
    /// Uppercase `G` base count.
    pub g: u64,
    /// Uppercase `T` base count.
    pub t: u64,
    /// Uppercase `N` base count.
    pub n: u64,
    /// Bases that are not uppercase A/C/G/T/N (IUPAC codes, lowercase, etc.).
    pub other: u64,
    /// Sum of every quality byte's Phred value (raw byte − 33), for the mean.
    pub qual_sum: u64,
    /// Count of quality bytes at each Phred value `0..QUAL_MAX`.
    pub qual_hist: [u64; QUAL_MAX],
}

impl Default for ContentStats {
    fn default() -> Self {
        ContentStats {
            reads: 0,
            bases: 0,
            min_len: 0,
            max_len: 0,
            a: 0,
            c: 0,
            g: 0,
            t: 0,
            n: 0,
            other: 0,
            qual_sum: 0,
            qual_hist: [0; QUAL_MAX],
        }
    }
}

impl ContentStats {
    /// GC fraction over unambiguous bases: `(G + C) / (A + C + G + T)`. `None`
    /// when there are no A/C/G/T bases (an all-`N` or empty archive).
    pub fn gc_fraction(&self) -> Option<f64> {
        let acgt = self.a + self.c + self.g + self.t;
        (acgt > 0).then(|| (self.g + self.c) as f64 / acgt as f64)
    }

    /// Mean read length, or `None` when there are no reads.
    pub fn mean_len(&self) -> Option<f64> {
        (self.reads > 0).then(|| self.bases as f64 / self.reads as f64)
    }

    /// Mean Phred quality over every base, or `None` when there are no bases.
    pub fn mean_quality(&self) -> Option<f64> {
        (self.bases > 0).then(|| self.qual_sum as f64 / self.bases as f64)
    }

    /// Whether every read has the same length (fixed-length run). True for the
    /// empty archive vacuously; check [`reads`](Self::reads) first if that
    /// matters.
    pub fn fixed_length(&self) -> bool {
        self.min_len == self.max_len
    }
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
        keep_order: header.flags & FLAG_REORDERED == 0 || header.flags & FLAG_KEEP_ORDER != 0,
        regenerated_names: header.flags & FLAG_REGEN_NAMES != 0,
        platform: Platform::from_code(header.platform),
        format_version: FORMAT_VERSION,
        ..Info::default()
    })
}

/// Read a `[u32 len][u32 crc][bytes]` frame's length and skip the CRC + bytes
/// without allocating them (for metadata-only scans; the CRC is not verified
/// since the payload is discarded). Returns the payload length.
pub(crate) fn skip_framed<R: Read>(r: &mut R) -> Result<usize> {
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
        keep_order: header.flags & FLAG_REORDERED == 0 || header.flags & FLAG_KEEP_ORDER != 0,
        regenerated_names: header.flags & FLAG_REGEN_NAMES != 0,
        platform: Platform::from_code(header.platform),
        format_version: FORMAT_VERSION,
        ..Info::default()
    };
    // Whole-file global-cluster layout: [u64 n][flip][perm][name template]
    // [global reference?][u32 n_blocks][seq blocks][name+qual blocks]. Permutation
    // is charged to seq; the name template (non-empty only in discard-order mode)
    // to names; the shared global reference (present iff FLAG_GLOBAL_REFERENCE) to
    // seq, since every version-4 block decodes against it.
    if header.flags & FLAG_GLOBAL_REORDER != 0 {
        let mut n8 = [0u8; 8];
        r.read_exact(&mut n8)?;
        info.reads = u64::from_le_bytes(n8);
        skip_framed(&mut r)?; // flip bitmap
        info.seq_bytes += skip_framed(&mut r)? as u64; // permutation
        info.names_bytes += skip_framed(&mut r)? as u64; // name template (regen mode)
        if header.flags & FLAG_GLOBAL_REFERENCE != 0 {
            info.seq_bytes += skip_framed(&mut r)? as u64; // shared global reference
        }
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
            info.whole_file_crc = Some(footer.whole_file_crc);
            // Per-stream sizes are recorded in the footer index itself (v3), so
            // summing them needs no per-block seeks — one footer read is the whole
            // metadata cost. The three streams are names, sequence, quality.
            for streams in &footer.stream_locs {
                info.names_bytes += u64::from(streams[0].len);
                info.seq_bytes += u64::from(streams[1].len);
                info.qual_bytes += u64::from(streams[2].len);
            }
        }
        Err(_) => scan_blocks_sequentially(&mut r, &mut info)?,
    }
    Ok(info)
}

/// Walk one block's header at the current position — the `[24] per-stream content
/// digests` prefix, `[4] n_reads`, an optional reorder preamble, and the three
/// `[4 len][bytes]` stream frames — accumulating per-stream sizes into `info` and
/// seeking past each payload. Returns the block's read count. Leaves the cursor at
/// the end of the block's payload.
pub(crate) fn scan_block_header<R: Read + Seek>(
    r: &mut R,
    reordered: bool,
    info: &mut Info,
) -> Result<u64> {
    // Skip the payload's leading content digests (see the block-payload layout).
    r.seek(SeekFrom::Current(STREAM_DIGESTS_LEN as i64))?;
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
pub(crate) fn scan_blocks_sequentially<R: Read + Seek>(r: &mut R, info: &mut Info) -> Result<()> {
    let mut off = HEADER_LEN as u64;
    loop {
        r.seek(SeekFrom::Start(off))?;
        // Frame head: [4] marker, [8] length, [4] CRC.
        let mut head = [0u8; FRAME_HEAD_LEN];
        if r.read_exact(&mut head).is_err() {
            break; // clean EOF (or a partial marker) at a frame boundary
        }
        if head[..BLOCK_MAGIC.len()] != BLOCK_MAGIC {
            break; // not a block boundary — corrupt, or past the block region
        }
        let plen = u64::from_le_bytes(head[4..12].try_into().unwrap());
        if plen == 0 || plen > MAX_BLOCK_PAYLOAD {
            break; // terminator, or a length too large to be a real frame
        }
        r.seek(SeekFrom::Start(off + FRAME_HEAD_LEN as u64))?;
        let n = match scan_block_header(r, info.reordered, info) {
            Ok(n) => n,
            Err(_) => break, // ran off the end of a truncated block
        };
        info.reads += n;
        info.blocks += 1;
        off += FRAME_HEAD_LEN as u64 + plen;
    }
    Ok(())
}
