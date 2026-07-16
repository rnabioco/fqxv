//! On-disk framing: header, footer/row-group index, and low-level frame IO.

use super::*;

/// Raw-sequence byte budget per row group. A group is cut at whichever comes
/// first — `block_reads` reads or this many raw sequence bytes — so long-read
/// (nanopore-style) data does not collapse into one enormous row group that
/// destroys parallelism and random-access granularity and could overflow the
/// `u32` per-stream compressed length. For fixed short reads the read count is
/// the binding limit and this never triggers.
pub(crate) const MAX_BLOCK_SEQ_BYTES: usize = 256 << 20;
/// Header fields covered by the header CRC: magic(4) + version(2) + seq_order(1)
/// + quality-binning tag(1) + flags(1) + group size(1) + platform(1).
///
/// The platform tag has its own byte rather than sharing the flags byte. It used
/// to occupy flags bits 5-7, which collided with [`FLAG_GLOBAL_REFERENCE`]
/// (0x20 == bit 5): `Platform::Illumina`'s code is 1, so its bits were exactly
/// the reference flag. An Illumina archive in the global-reorder layout that did
/// not adopt a reference therefore still advertised one and failed to decode.
/// Six flags plus a 3-bit platform tag do not fit in one byte, so the tag moved
/// out. (Pre-existing archives are not a concern — the format is not yet stable
/// and the widened CRC coverage rejects a stale header rather than misreading
/// it.)
pub(crate) const HEADER_FIELDS_LEN: usize = 11;
/// Full on-disk header prefix: the [`HEADER_FIELDS_LEN`] fields followed by a
/// CRC-32C over them, so a flipped header byte (version, flags, the lossy binning
/// tag, group size) is caught on read rather than silently changing decode — in
/// both the plain and reorder layouts. Also the byte offset at which the first
/// block / reorder frame begins.
pub(crate) const HEADER_LEN: usize = HEADER_FIELDS_LEN + CRC_LEN;
/// Bytes of CRC-32C appended after a frame's length field (plain block frames)
/// or after a `[u32 len]` framed slice (reorder layout).
pub(crate) const CRC_LEN: usize = 4;
/// Bytes of xxh3-64 content digest prepended to each plain block payload — an
/// end-to-end round-trip check over the block's *decoded* content, distinct from
/// the frame CRC (which only covers stored/compressed bytes). See
/// [`content_digest`].
pub(crate) const DIGEST_LEN: usize = 8;
/// Upper bound on a single block payload's declared length. A block holds at most
/// `block_reads` reads and `MAX_BLOCK_SEQ_BYTES` of raw sequence, and the three
/// compressed streams are each smaller than their raw input in the common case,
/// so a real payload is comfortably under this. It exists only so a corrupted
/// length field can't drive a multi-exabyte allocation before the CRC is even
/// checked; anything larger is rejected as malformed.
pub(crate) const MAX_BLOCK_PAYLOAD: u64 = (MAX_BLOCK_SEQ_BYTES as u64) * 8;
/// Magic at the very end of a v1 archive, just after the `[8] footer_offset`
/// back-pointer, so a reader can confirm it found a real footer.
pub(crate) const FOOTER_MAGIC: [u8; 4] = *b"FQXF";
/// Bytes in the fixed EOF trailer: `[8] footer_offset` + `[4] FOOTER_MAGIC`.
pub(crate) const TRAILER_LEN: usize = 12;
/// Per-block sync marker at the head of every block frame (`[4] BLOCK_MAGIC
/// [8] len [4] crc [payload]`). Recovery scans for it to resynchronize to a block
/// boundary, so a corrupt length prefix or a lost footer no longer strands the
/// blocks that follow it (see [`decompress_recover`](super::decompress_recover)).
pub(crate) const BLOCK_MAGIC: [u8; 4] = *b"FQXB";
/// Frame overhead before a block payload: `BLOCK_MAGIC + [8] len + [4] crc`.
pub(crate) const FRAME_HEAD_LEN: usize = BLOCK_MAGIC.len() + 8 + CRC_LEN;
pub(crate) const FLAG_PLUS_NORMALIZED: u8 = 0x01;
pub(crate) const FLAG_REORDERED: u8 = 0x02;
pub(crate) const FLAG_KEEP_ORDER: u8 = 0x04;
/// Whole-file, globally-clustered reorder layout (see `compress_reordered_whole`)
/// as opposed to the older per-block reorder blocks.
pub(crate) const FLAG_GLOBAL_REORDER: u8 = 0x08;
/// Names are regenerated from a stored counter template, not coded per read
/// (reorder-lossy: reads are renumbered). Set only in the discard-order layout.
pub(crate) const FLAG_REGEN_NAMES: u8 = 0x10;
/// The whole-file reorder layout carries a shared frozen global reference frame
/// (SPRING-style), and one or more sequence blocks are version-4 positions on it
/// (see `fqxv_reorder::GlobalReference`). Only set within a `FLAG_GLOBAL_REORDER`
/// archive, and only when the reference nets a whole-file byte win over the
/// block-local (v2/v3) codecs; otherwise no reference frame is written.
pub(crate) const FLAG_GLOBAL_REFERENCE: u8 = 0x20;
// Bits 0x40 and 0x80 are free: the platform tag that used to occupy bits 5-7
// now has its own header byte (see `HEADER_FIELDS_LEN`).

/// Write the `HEADER_FIELDS_LEN` header fields followed by a CRC-32C over them,
/// so a flipped header byte is caught on read instead of silently altering decode.
/// Shared by the plain ([`write_header`]) and reorder ([`encode_reordered`])
/// layouts, which differ only in their `flags`.
pub(crate) fn write_header_prefix<W: Write>(
    w: &mut W,
    seq_order: u8,
    binning: u8,
    flags: u8,
    group_size: u8,
    platform: Platform,
) -> Result<()> {
    let mut hdr = [0u8; HEADER_FIELDS_LEN];
    hdr[..4].copy_from_slice(&MAGIC);
    hdr[4..6].copy_from_slice(&FORMAT_VERSION.to_le_bytes());
    hdr[6] = seq_order;
    hdr[7] = binning;
    hdr[8] = flags;
    hdr[9] = group_size;
    hdr[10] = platform.to_code();
    w.write_all(&hdr)?;
    w.write_all(&crc32c(&hdr).to_le_bytes())?;
    Ok(())
}

/// Write the container header (plain layout).
pub(crate) fn write_header<W: Write>(
    w: &mut W,
    params: &Params,
    group_size: u8,
    platform: Platform,
) -> Result<()> {
    // The block layout is always non-reorder — reorder (both keep-order modes)
    // uses the whole-file path, which writes its own header.
    debug_assert!(!params.reorder);
    write_header_prefix(
        w,
        params.seq_order,
        binning_tag(params.quality_binning),
        FLAG_PLUS_NORMALIZED,
        group_size,
        platform,
    )
}

/// Row-group index accumulated as blocks are written, serialized into the v1
/// footer. `offset` tracks the current byte position from the start of the file,
/// so each entry records where its row group's `[8] length` field begins.
pub(crate) struct FooterIndex {
    pub(crate) entries: Vec<(u64, u32)>,
    pub(crate) offset: u64,
}

impl FooterIndex {
    pub(crate) fn new() -> Self {
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
pub(crate) fn write_footer<W: Write>(
    w: &mut CrcWriter<W>,
    index: &FooterIndex,
    total_reads: u64,
) -> Result<u64> {
    // Zero-length terminator block — a full frame (marker + `len == 0`) so the
    // scan-based recovery treats it like any other boundary. Fed through the tee.
    w.write_all(&BLOCK_MAGIC)?;
    w.write_all(&0u64.to_le_bytes())?;
    let footer_offset = index.offset + BLOCK_MAGIC.len() as u64 + 8;

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
    // body = n_groups..whole_file_crc; +CRC_LEN footer_crc, + marker+len terminator.
    Ok((BLOCK_MAGIC.len() + 8) as u64 + body.len() as u64 + CRC_LEN as u64 + TRAILER_LEN as u64)
}

/// Write `[u32 len][u32 crc32c][bytes]`.
pub(crate) fn write_framed<W: Write>(w: &mut W, bytes: &[u8]) -> Result<()> {
    // The reorder layout frames have no raw-byte budget, so guard the u32 length
    // cast explicitly: a >4 GiB frame must error, not silently truncate.
    let len = u32::try_from(bytes.len())
        .map_err(|_| Error::Malformed("framed slice exceeds u32 length"))?;
    w.write_all(&len.to_le_bytes())?;
    w.write_all(&crc32c(bytes).to_le_bytes())?;
    w.write_all(bytes)?;
    Ok(())
}

/// Read a `[u32 len][u32 crc32c][bytes]` frame, guarding the length allocation and
/// verifying the CRC. `what` names the frame in any corruption error.
pub(crate) fn read_framed<R: Read>(r: &mut R, what: &str) -> Result<Vec<u8>> {
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

// --- header / block framing --------------------------------------------------

pub(crate) struct Header {
    pub(crate) seq_order: u8,
    pub(crate) quality_binning: u8,
    pub(crate) flags: u8,
    pub(crate) group_size: u8,
    /// Sequencing platform tag ([`Platform::to_code`]) — its own byte, never
    /// packed into `flags` (see [`HEADER_FIELDS_LEN`]).
    pub(crate) platform: u8,
}

pub(crate) fn read_header<R: Read>(r: &mut R) -> Result<Header> {
    let mut buf = [0u8; HEADER_LEN];
    r.read_exact(&mut buf)?;
    let (fields, crc) = buf.split_at(HEADER_FIELDS_LEN);
    if fields[..4] != MAGIC {
        return Err(Error::BadMagic);
    }
    let ver = u16::from_le_bytes([fields[4], fields[5]]);
    if ver != FORMAT_VERSION {
        return Err(Error::UnsupportedVersion(ver));
    }
    // Verify the header CRC only after magic/version, so a genuinely foreign or
    // wrong-version file reports that (more useful) error rather than a CRC
    // mismatch. Within this version, a flipped field byte is caught here before
    // it can silently change decode (group size, flags, the lossy binning tag).
    if u32::from_le_bytes(crc.try_into().unwrap()) != crc32c(fields) {
        return Err(Error::Corrupt {
            what: "header".to_string(),
        });
    }
    let group_size = fields[9].max(1);
    Ok(Header {
        seq_order: fields[6],
        quality_binning: fields[7],
        flags: fields[8],
        group_size,
        platform: fields[10],
    })
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
pub(crate) struct Footer {
    /// `(byte offset of the row group's [8] length field, its read count)`.
    pub(crate) groups: Vec<(u64, u32)>,
    pub(crate) total_reads: u64,
    /// CRC-32C of the archive from byte 0 through the `total_reads` field.
    pub(crate) whole_file_crc: u32,
    /// Number of leading bytes that `whole_file_crc` covers (the offset of the
    /// `whole_file_crc` field itself), so a verifier can re-hash exactly that
    /// prefix.
    pub(crate) covered_len: u64,
}

/// Bytes appended to the footer body after `total_reads`: `[4 whole_file_crc]
/// [4 footer_crc]`.
pub(crate) const FOOTER_CRC_TAIL: usize = 8;

/// Read the footer by seeking to the EOF trailer and following its back-pointer,
/// then verify the footer's own CRC before trusting any offset in it.
pub(crate) fn read_footer<R: Read + Seek>(r: &mut R) -> Result<Footer> {
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
pub(crate) fn read_u32<R: Read>(r: &mut R) -> Result<u32> {
    let mut b = [0u8; 4];
    r.read_exact(&mut b)?;
    Ok(u32::from_le_bytes(b))
}

pub(crate) fn build_pool(threads: usize) -> Result<rayon::ThreadPool> {
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

pub(crate) fn binning_tag(b: QualityBinning) -> u8 {
    match b {
        QualityBinning::Lossless => 0,
        QualityBinning::Bin8 => 1,
        QualityBinning::Bin4 => 2,
        QualityBinning::Bin2 => 3,
    }
}

/// Shared byte cursor specialized to this crate's [`Error`].
///
/// The container reads fixed-width little-endian fields and length-prefixed
/// slices; every out-of-bounds read maps to [`Error::Truncated`].
pub(crate) type Cursor<'a> = fqxv_bytes::Reader<'a, Error>;

impl fqxv_bytes::ReaderError for Error {
    fn truncated() -> Self {
        Error::Truncated
    }
    fn bad_varint() -> Self {
        Error::Truncated
    }
    fn oversized() -> Self {
        Error::Truncated
    }
}
