//! On-disk framing: header, footer/row-group index, and low-level frame IO.

use super::*;

/// Raw-sequence byte budget per row group. A group is cut at whichever comes
/// first — `block_reads` reads or this many raw sequence bytes — so long-read
/// (nanopore-style) data does not collapse into one enormous row group that
/// destroys parallelism and random-access granularity and could overflow the
/// `u32` per-stream compressed length. For fixed short reads the read count is
/// the binding limit and this never triggers.
pub(crate) const MAX_BLOCK_SEQ_BYTES: usize = 256 << 20;
/// Fixed header prefix, byte offsets within it, up to and including the 2-byte
/// extension-length field. Layout:
///
/// ```text
/// [0..4]   magic "FQXV"
/// [4]      version major        -- reader refuses a differing major
/// [5]      version minor        -- tolerated within a major; informational
/// [6..14]  required_features    -- u64 LE; a set bit outside KNOWN_FEATURES is
///                                  refused ([`Error::UnsupportedFeature`])
/// [14]     sequence context order (k)
/// [15]     quality-binning tag
/// [16]     flags
/// [17]     group size G
/// [18]     platform tag
/// [19..21] ext_len              -- u16 LE; bytes of the TLV extension region
/// ```
///
/// The platform tag has its own byte rather than sharing the flags byte. It used
/// to occupy flags bits 5-7, which collided with [`FLAG_GLOBAL_REFERENCE`]
/// (0x20 == bit 5): `Platform::Illumina`'s code is 1, so its bits were exactly
/// the reference flag. An Illumina archive in the global-reorder layout that did
/// not adopt a reference therefore still advertised one and failed to decode.
/// Six flags plus a 3-bit platform tag do not fit in one byte, so the tag moved
/// out.
pub(crate) const HEADER_PREFIX_LEN: usize = 21;
/// Byte offset of the quality-binning tag within the header (see the layout
/// above). Exposed so tests that poke at a raw header field track the layout by
/// name instead of a bare literal that silently rots when the field shifts.
#[cfg(test)]
pub(crate) const HDR_OFF_BINNING: usize = 15;
/// Byte offset of the flags byte within the header.
#[cfg(test)]
pub(crate) const HDR_OFF_FLAGS: usize = 16;
/// Header length for an archive with an empty extension region — the size this
/// build writes (1.0 emits no TLV records). A later minor that appends TLVs makes
/// the on-disk header longer; readers must use [`Header::header_len`] (not this
/// constant) for any seek/scan start offset so they stay correct across minors.
/// Blocks begin right after the header, so this is also the ext-empty block-region
/// start offset.
pub(crate) const HEADER_LEN: usize = HEADER_PREFIX_LEN + CRC_LEN;
/// A header extension record is `[1 tag][2 len LE][len bytes]`. The tag's high bit
/// marks a *critical* record: a reader that doesn't recognize a critical tag must
/// refuse the archive ([`Error::UnsupportedExtension`]) rather than skip it; an
/// unknown non-critical tag is skipped. No tags are defined at 1.0 — the region
/// exists so a later minor can add skippable header fields without a major bump.
pub(crate) const EXT_CRITICAL_BIT: u8 = 0x80;
/// Bytes of CRC-32C appended after a frame's length field (plain block frames)
/// or after a `[u32 len]` framed slice (reorder layout).
pub(crate) const CRC_LEN: usize = 4;
/// Bytes of one xxh3-64 digest. Also the width of the reorder layout's single
/// `output_digest`.
pub(crate) const DIGEST_LEN: usize = 8;
/// Bytes of the three per-stream content digests (names, sequence, quality)
/// prepended to each plain block payload — an end-to-end round-trip check over the
/// block's *decoded* content, one digest per stream so a mismatch localizes which
/// stream a codec round-tripped wrong. Distinct from the frame CRC (which only
/// covers stored/compressed bytes). See [`stream_digests`].
pub(crate) const STREAM_DIGESTS_LEN: usize = 3 * DIGEST_LEN;
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
// now has its own header byte (see `HEADER_PREFIX_LEN`).

/// Write the fixed header prefix (see [`HEADER_PREFIX_LEN`]), then the extension
/// region (empty at 1.0), then a CRC-32C over prefix+extension — so a flipped
/// header byte is caught on read instead of silently altering decode. Shared by
/// the plain ([`write_header`]) and reorder (`encode_reordered`) layouts, which
/// differ only in their `flags`.
///
/// `required_features` is the coarse, intent-level capability word: set a bit
/// here only for a capability a reader must possess to decode the archive *at
/// all* and that is knowable before the blocks are written (e.g. a whole-file
/// global reference frame). Fine-grained "this one block used codec X" is carried
/// by the per-stream method byte and surfaces as an [`Error::UnsupportedMethod`]
/// on decode, not here. Pass 0 when the archive requires nothing beyond the base
/// format for its major.
pub(crate) fn write_header_prefix<W: Write>(
    w: &mut W,
    seq_order: u8,
    binning: u8,
    flags: u8,
    group_size: u8,
    platform: Platform,
    required_features: u64,
) -> Result<()> {
    // 1.0 emits an empty extension region; a later minor appends TLV records here
    // and must update HEADER_LEN-based start offsets to use `Header::header_len`.
    const EXT: &[u8] = &[];
    let mut hdr = Vec::with_capacity(HEADER_PREFIX_LEN + EXT.len());
    hdr.extend_from_slice(&MAGIC);
    hdr.push(FORMAT_MAJOR);
    hdr.push(FORMAT_MINOR);
    hdr.extend_from_slice(&required_features.to_le_bytes());
    hdr.push(seq_order);
    hdr.push(binning);
    hdr.push(flags);
    hdr.push(group_size);
    hdr.push(platform.to_code());
    hdr.extend_from_slice(&(EXT.len() as u16).to_le_bytes());
    hdr.extend_from_slice(EXT);
    debug_assert_eq!(hdr.len(), HEADER_PREFIX_LEN + EXT.len());
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
    // The plain layout needs nothing beyond the base format for its major, with one
    // exception: sequence-only (`no_quality`) archives set the coarse
    // `NO_QUALITY` feature bit so an older reader refuses them rather than
    // reconstructing quality-less records into mis-framed FASTQ. Per-block codec
    // choices (e.g. the long-read overlap codec) are still recorded by the sequence
    // stream's method byte and rejected per block on decode, not gated here.
    let required_features = if params.no_quality {
        crate::feature::NO_QUALITY
    } else {
        0
    };
    write_header_prefix(
        w,
        params.seq_order,
        binning_tag(params.quality_binning),
        FLAG_PLUS_NORMALIZED,
        group_size,
        platform,
        required_features,
    )
}

/// Absolute on-disk location of one coded stream (names, sequence, or quality)
/// within a block, recorded in the footer index so a remote client can fetch and
/// verify a single stream without reading the whole block.
///
/// `offset` is the absolute byte offset of the stream's *coded* bytes (past the
/// block frame head, the three content digests, `n_reads`, and this stream's length
/// prefix); `len` is that coded length; `crc` is CRC-32C of exactly those `len`
/// bytes. The block's content digests cover the *decoded* streams and so cannot be
/// checked on a projected fetch of *coded* bytes — this per-stream CRC is the
/// substitute that keeps a single-stream read verifiable.
#[derive(Debug, Clone, Copy, Default)]
pub(crate) struct StreamLoc {
    pub(crate) offset: u64,
    pub(crate) len: u32,
    pub(crate) crc: u32,
}

/// Per-group footer bytes: `[8 off][4 read_count]` plus a `[8 off][4 len][4 crc]`
/// triple for each of the three streams (names, sequence, quality).
pub(crate) const FOOTER_GROUP_BYTES: usize = 12 + 3 * 16;

/// Row-group index accumulated as blocks are written, serialized into the v1
/// footer. `offset` tracks the current byte position from the start of the file,
/// so each entry records where its row group's frame marker begins. `streams` is
/// parallel to `entries`: `streams[i]` locates group `i`'s three coded streams.
pub(crate) struct FooterIndex {
    pub(crate) entries: Vec<(u64, u32)>,
    pub(crate) streams: Vec<[StreamLoc; 3]>,
    pub(crate) offset: u64,
}

impl FooterIndex {
    pub(crate) fn new() -> Self {
        // Blocks begin right after the fixed header.
        Self::new_at(HEADER_LEN as u64)
    }

    /// Like [`FooterIndex::new`] but with the first block at an explicit offset —
    /// used by the shared-reference plain layout (issue #168), where a whole-file
    /// reference frame sits between the header and the first block, so block
    /// offsets recorded in the footer must start past it.
    pub(crate) fn new_at(offset: u64) -> Self {
        FooterIndex {
            entries: Vec::new(),
            streams: Vec::new(),
            offset,
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

    let mut body =
        Vec::with_capacity(4 + index.entries.len() * FOOTER_GROUP_BYTES + 8 + FOOTER_CRC_TAIL);
    body.extend_from_slice(&(index.entries.len() as u32).to_le_bytes());
    for (&(off, read_count), streams) in index.entries.iter().zip(&index.streams) {
        body.extend_from_slice(&off.to_le_bytes());
        body.extend_from_slice(&read_count.to_le_bytes());
        for s in streams {
            body.extend_from_slice(&s.offset.to_le_bytes());
            body.extend_from_slice(&s.len.to_le_bytes());
            body.extend_from_slice(&s.crc.to_le_bytes());
        }
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

/// Read the plain layout's whole-file reference frame if the header's
/// `FLAG_GLOBAL_REFERENCE` bit is set (issue #168). The frame is a single
/// [`write_framed`] slice immediately after the header, before the first block; a
/// reader positioned just past the header consumes it here and threads the decoded
/// [`fqxv_lroverlap::Reference`] into every block's sequence decode. Returns `None`
/// when the bit is clear (no frame is present). A corrupt frame fails closed.
pub(crate) fn read_reference_frame<R: Read>(
    r: &mut R,
    flags: u8,
) -> Result<Option<fqxv_lroverlap::Reference>> {
    if flags & FLAG_GLOBAL_REFERENCE == 0 {
        return Ok(None);
    }
    let bytes = read_framed(r, "plain-layout global reference")?;
    Ok(Some(fqxv_lroverlap::Reference::decode(&bytes)?))
}

/// Read a `[u32 len][u32 crc32c][bytes]` frame, guarding the length allocation and
/// verifying the CRC. `what` names the frame in any corruption error.
pub(crate) fn read_framed<R: Read>(r: &mut R, what: &str) -> Result<Vec<u8>> {
    let mut lb = [0u8; 4];
    r.read_exact(&mut lb)?;
    let len = u32::from_le_bytes(lb) as u64;
    // Cap the declared length like `read_block` does. `read_framed` drives every
    // section of the reorder decode path (flip bitmap, permutation, template,
    // global reference, per-block frames), all reached before any content check,
    // so an uncapped u32 `len` let a ~30-byte archive claim ~4 GB (#142).
    if len > MAX_BLOCK_PAYLOAD {
        return Err(Error::Malformed("framed slice length exceeds the maximum"));
    }
    let len = len as usize;
    let mut cb = [0u8; CRC_LEN];
    r.read_exact(&mut cb)?;
    let expected = u32::from_le_bytes(cb);
    // Read incrementally rather than `resize(len, 0)` + `read_exact`: a truncated
    // stream that claims the maximum length then allocates only what is actually
    // present, not the full claim zero-filled up front. `take` bounds the read so
    // the CRC still guards the payload, and a short read surfaces as `Truncated`.
    let mut buf = Vec::new();
    let got = r.by_ref().take(len as u64).read_to_end(&mut buf)?;
    if got != len {
        return Err(Error::Truncated);
    }
    if crc32c(&buf) != expected {
        return Err(Error::Corrupt {
            what: what.to_string(),
        });
    }
    Ok(buf)
}

// --- header / block framing --------------------------------------------------

pub(crate) struct Header {
    /// On-disk format major/minor. `major` always equals [`FORMAT_MAJOR`] for a
    /// readable archive (a differing major is refused); `minor` may be newer than
    /// this build's ([`FORMAT_MINOR`]) — additions within a major are tolerated.
    pub(crate) major: u8,
    pub(crate) minor: u8,
    /// Coarse capability word ([`crate::feature`]); every set bit is within
    /// [`crate::KNOWN_FEATURES`] for a readable archive (unknown bits are refused).
    pub(crate) required_features: u64,
    pub(crate) seq_order: u8,
    pub(crate) quality_binning: u8,
    pub(crate) flags: u8,
    pub(crate) group_size: u8,
    /// Sequencing platform tag ([`Platform::to_code`]) — its own byte, never
    /// packed into `flags` (see [`HEADER_PREFIX_LEN`]).
    pub(crate) platform: u8,
    /// Actual on-disk header length in bytes: prefix + extension region + CRC.
    /// Equals [`HEADER_LEN`] for a 1.0-written archive (empty extension region);
    /// larger if a later minor appended TLV records. Any seek/scan that starts at
    /// the block region must use this, not the [`HEADER_LEN`] constant.
    pub(crate) header_len: u64,
}

pub(crate) fn read_header<R: Read>(r: &mut R) -> Result<Header> {
    let mut prefix = [0u8; HEADER_PREFIX_LEN];
    // A file too short to even hold the fixed header prefix (empty, a stray byte, a
    // badly truncated download) surfaces `read_exact`'s raw "failed to fill whole
    // buffer" otherwise — map it to the clear `Truncated` message the rest of the
    // decode path uses, keeping the corruption diagnostics consistent.
    r.read_exact(&mut prefix).map_err(|e| {
        if e.kind() == std::io::ErrorKind::UnexpectedEof {
            Error::Truncated
        } else {
            Error::Io(e)
        }
    })?;
    if prefix[..4] != MAGIC {
        return Err(Error::BadMagic);
    }
    let major = prefix[4];
    let minor = prefix[5];
    // Refuse a wire-incompatible major before anything else, so a genuinely
    // foreign or future-major file reports that (most useful) error rather than a
    // CRC mismatch. A newer minor within our major is tolerated — its additions
    // ride in the skippable extension region and per-block method bytes.
    if major != FORMAT_MAJOR {
        return Err(Error::UnsupportedVersion { major, minor });
    }
    let required_features = u64::from_le_bytes(prefix[6..14].try_into().unwrap());
    // A required capability this build doesn't have means the archive cannot be
    // decoded — refuse cleanly, naming the unknown bits, rather than mis-reading a
    // stream produced by a codec/mode we lack.
    let unknown = required_features & !KNOWN_FEATURES;
    if unknown != 0 {
        return Err(Error::UnsupportedFeature(unknown));
    }
    let ext_len = u16::from_le_bytes([prefix[19], prefix[20]]) as usize;

    // The extension region is bounded by a u16 so this can't over-allocate.
    let mut ext = vec![0u8; ext_len];
    r.read_exact(&mut ext)?;
    let mut crcb = [0u8; CRC_LEN];
    r.read_exact(&mut crcb)?;

    // Verify the header CRC (over prefix + extension) only after magic/version/
    // features, so those errors win over a CRC mismatch. Within a compatible
    // header a flipped field byte is caught here before it can silently change
    // decode (group size, flags, the lossy binning tag, a feature bit).
    let mut covered = Vec::with_capacity(HEADER_PREFIX_LEN + ext.len());
    covered.extend_from_slice(&prefix);
    covered.extend_from_slice(&ext);
    if u32::from_le_bytes(crcb) != crc32c(&covered) {
        return Err(Error::Corrupt {
            what: "header".to_string(),
        });
    }
    // Walk the TLV records: an unknown *critical* tag is fatal; unknown
    // non-critical tags are skipped. (No tags are defined at 1.0.)
    check_header_extensions(&ext)?;

    let group_size = prefix[17].max(1);
    Ok(Header {
        major,
        minor,
        required_features,
        seq_order: prefix[14],
        quality_binning: prefix[15],
        flags: prefix[16],
        group_size,
        platform: prefix[18],
        header_len: (HEADER_PREFIX_LEN + ext_len + CRC_LEN) as u64,
    })
}

/// Walk the header extension region's `[1 tag][2 len][len bytes]` records. No
/// tags are known at 1.0, so a known tag is never consumed here yet; an unknown
/// *critical* tag (high bit set, [`EXT_CRITICAL_BIT`]) is refused, and an unknown
/// non-critical tag is skipped. A record that overruns the region is malformed.
fn check_header_extensions(mut ext: &[u8]) -> Result<()> {
    while !ext.is_empty() {
        if ext.len() < 3 {
            return Err(Error::Malformed("truncated header extension record"));
        }
        let tag = ext[0];
        let len = u16::from_le_bytes([ext[1], ext[2]]) as usize;
        let end = 3 + len;
        if end > ext.len() {
            return Err(Error::Malformed("header extension record overruns region"));
        }
        if tag & EXT_CRITICAL_BIT != 0 {
            return Err(Error::UnsupportedExtension(tag));
        }
        ext = &ext[end..];
    }
    Ok(())
}

/// The footer index: per row group `(byte_offset, read_count)`, the total read
/// count, and the whole-archive CRC-32C, located via the EOF trailer's
/// back-pointer. Only the plain and per-block-reorder layouts carry a footer (the
/// globally-clustered layout is self-describing — see the module docs).
///
/// On-disk body (at `footer_offset`, i.e. just past the terminator):
/// `[4 n_groups] ([8 off][4 read_count] ([8 s_off][4 s_len][4 s_crc])×3)*
///  [8 total_reads] [4 whole_file_crc] [4 footer_crc]`.
/// Each row group records its block offset and read count, then a
/// `(offset, len, crc32c)` triple locating each of its three coded streams
/// (names, sequence, quality) for remote column projection. `footer_crc` covers
/// the body up to (not including) itself, so the index can be trusted without
/// rereading the whole archive; `whole_file_crc` covers every byte from the header
/// through `total_reads` and is checked only by the whole-archive verify/recover
/// path.
pub(crate) struct Footer {
    /// `(byte offset of the row group's frame marker, its read count)`.
    pub(crate) groups: Vec<(u64, u32)>,
    /// Parallel to `groups`: `stream_locs[i]` locates group `i`'s three coded
    /// streams (names, sequence, quality) for single-stream projection.
    pub(crate) stream_locs: Vec<[StreamLoc; 3]>,
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
    parse_footer_body(&body, footer_offset)
}

/// Parse and CRC-check a footer body already in memory. `body` is the bytes from
/// `footer_offset` up to (not including) the EOF trailer — the `[4 n_groups] …
/// [8 total_reads] [4 whole_file_crc] [4 footer_crc]` region. Shared by
/// `read_footer` (seek-based) and the IO-free suffix parser in `random_access`,
/// so both apply exactly the same validation before any offset is trusted.
pub(crate) fn parse_footer_body(body: &[u8], footer_offset: u64) -> Result<Footer> {
    let body_len = body.len();
    if body_len < 4 + 8 + FOOTER_CRC_TAIL {
        return Err(Error::Malformed("footer body too short"));
    }
    let (covered, footer_crc_bytes) = body.split_at(body_len - CRC_LEN);
    let footer_crc = u32::from_le_bytes(footer_crc_bytes.try_into().unwrap());
    if crc32c(covered) != footer_crc {
        return Err(Error::Corrupt {
            what: "footer".to_string(),
        });
    }

    let mut c = Cursor::new(covered);
    let n_groups = c.u32()? as usize;
    // `covered` = [4 n_groups][FOOTER_GROUP_BYTES*n_groups][8 total_reads][4
    // whole_file_crc]; the CRC just passed already implies a self-consistent
    // length, but bound the allocation independently in case of a hash collision.
    let max_groups = covered.len().saturating_sub(4 + 8 + CRC_LEN) / FOOTER_GROUP_BYTES;
    if n_groups > max_groups {
        return Err(Error::Malformed("footer group count exceeds footer size"));
    }
    let mut groups = Vec::with_capacity(n_groups);
    let mut stream_locs = Vec::with_capacity(n_groups);
    for _ in 0..n_groups {
        let off = u64::from_le_bytes(c.take(8)?.try_into().unwrap());
        let rc = u32::from_le_bytes(c.take(4)?.try_into().unwrap());
        if off < HEADER_LEN as u64 || off >= footer_offset {
            return Err(Error::Malformed("footer row-group offset out of range"));
        }
        let mut streams = [StreamLoc::default(); 3];
        for s in &mut streams {
            let s_off = u64::from_le_bytes(c.take(8)?.try_into().unwrap());
            let s_len = u32::from_le_bytes(c.take(4)?.try_into().unwrap());
            let s_crc = u32::from_le_bytes(c.take(4)?.try_into().unwrap());
            // A stream lives inside its block, which starts at `off` and ends
            // before the footer; reject an index that points a stream outside the
            // block region before any caller dereferences it over the wire.
            let end = s_off
                .checked_add(u64::from(s_len))
                .ok_or(Error::Malformed("footer stream offset overflow"))?;
            if s_off < off || end > footer_offset {
                return Err(Error::Malformed("footer stream location out of range"));
            }
            *s = StreamLoc {
                offset: s_off,
                len: s_len,
                crc: s_crc,
            };
        }
        groups.push((off, rc));
        stream_locs.push(streams);
    }
    let total_reads = u64::from_le_bytes(c.take(8)?.try_into().unwrap());
    let whole_file_crc = c.u32()?;
    // whole_file_crc covers everything up to its own field: the archive prefix
    // plus the footer body through total_reads.
    let covered_len = footer_offset + (body_len - FOOTER_CRC_TAIL) as u64;
    Ok(Footer {
        groups,
        stream_locs,
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

/// Resolve the effective worker count: 0 means "all available cores", and any
/// explicit request is clamped to what physically exists so we never
/// over-subscribe.
pub(crate) fn resolve_threads(threads: usize) -> usize {
    let available = std::thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    if threads == 0 {
        available
    } else {
        threads.min(available)
    }
    .max(1)
}

pub(crate) fn build_pool(threads: usize) -> Result<rayon::ThreadPool> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(resolve_threads(threads))
        .build()
        .map_err(|e| Error::Io(io::Error::other(e.to_string())))
}

pub(crate) fn binning_tag(b: QualityBinning) -> u8 {
    match b {
        QualityBinning::Lossless => 0,
        QualityBinning::Bin8 => 1,
        QualityBinning::Bin4 => 2,
        QualityBinning::Bin2 => 3,
        QualityBinning::BinOnt => 4,
        QualityBinning::BinHifi => 5,
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

#[cfg(test)]
mod header_tests {
    use super::*;
    use std::io::Cursor;

    /// Build a header with arbitrary major/minor/features/extension bytes and a
    /// correct trailing CRC — the low-level mirror of [`write_header_prefix`] used
    /// to forge headers a current encoder would never emit.
    fn forge(major: u8, minor: u8, features: u64, ext: &[u8]) -> Vec<u8> {
        let mut hdr = Vec::new();
        hdr.extend_from_slice(&MAGIC);
        hdr.push(major);
        hdr.push(minor);
        hdr.extend_from_slice(&features.to_le_bytes());
        hdr.push(3); // seq_order
        hdr.push(0); // binning
        hdr.push(FLAG_PLUS_NORMALIZED);
        hdr.push(2); // group_size
        hdr.push(1); // platform (Illumina)
        hdr.extend_from_slice(&(ext.len() as u16).to_le_bytes());
        hdr.extend_from_slice(ext);
        hdr.extend_from_slice(&crc32c(&hdr).to_le_bytes());
        hdr
    }

    #[test]
    fn roundtrips_current_version() {
        let mut buf = Vec::new();
        write_header_prefix(
            &mut buf,
            3,
            0,
            FLAG_PLUS_NORMALIZED,
            2,
            Platform::Illumina,
            0,
        )
        .unwrap();
        assert_eq!(buf.len(), HEADER_LEN, "1.0 writes an ext-empty header");
        let h = read_header(&mut Cursor::new(&buf)).unwrap();
        assert_eq!((h.major, h.minor), (FORMAT_MAJOR, FORMAT_MINOR));
        assert_eq!(h.seq_order, 3);
        assert_eq!(h.group_size, 2);
        assert_eq!(h.platform, 1);
        assert_eq!(h.required_features, 0);
        assert_eq!(h.header_len, HEADER_LEN as u64);
    }

    #[test]
    fn refuses_foreign_major() {
        let buf = forge(FORMAT_MAJOR + 1, 0, 0, &[]);
        assert!(matches!(
            read_header(&mut Cursor::new(&buf)),
            Err(Error::UnsupportedVersion { major, .. }) if major == FORMAT_MAJOR + 1
        ));
    }

    #[test]
    fn tolerates_newer_minor() {
        // A newer minor within our major must read: additions are additive.
        let buf = forge(FORMAT_MAJOR, FORMAT_MINOR + 7, 0, &[]);
        let h = read_header(&mut Cursor::new(&buf)).unwrap();
        assert_eq!(h.minor, FORMAT_MINOR + 7);
    }

    #[test]
    fn refuses_unknown_required_feature() {
        let unknown = 1u64 << 40; // outside KNOWN_FEATURES
        let buf = forge(FORMAT_MAJOR, FORMAT_MINOR, unknown, &[]);
        assert!(matches!(
            read_header(&mut Cursor::new(&buf)),
            Err(Error::UnsupportedFeature(bits)) if bits == unknown
        ));
    }

    #[test]
    fn accepts_no_quality_feature() {
        // A current reader knows NO_QUALITY and must accept it, surfacing the bit;
        // a reader without it in KNOWN_FEATURES refuses via the mechanism above.
        let buf = forge(FORMAT_MAJOR, FORMAT_MINOR, crate::feature::NO_QUALITY, &[]);
        let h = read_header(&mut Cursor::new(&buf)).unwrap();
        assert_ne!(h.required_features & crate::feature::NO_QUALITY, 0);
    }

    #[test]
    fn skips_unknown_noncritical_extension() {
        // tag 0x01 (critical bit clear), 2-byte payload — a later minor's field.
        let ext = [0x01, 0x02, 0x00, 0xaa, 0xbb];
        let buf = forge(FORMAT_MAJOR, FORMAT_MINOR, 0, &ext);
        let h = read_header(&mut Cursor::new(&buf)).unwrap();
        // Block region starts past the extension region, not at the const.
        assert_eq!(
            h.header_len,
            (HEADER_PREFIX_LEN + ext.len() + CRC_LEN) as u64
        );
    }

    #[test]
    fn refuses_unknown_critical_extension() {
        let ext = [EXT_CRITICAL_BIT | 0x01, 0x00, 0x00];
        let buf = forge(FORMAT_MAJOR, FORMAT_MINOR, 0, &ext);
        assert!(matches!(
            read_header(&mut Cursor::new(&buf)),
            Err(Error::UnsupportedExtension(tag)) if tag == EXT_CRITICAL_BIT | 0x01
        ));
    }

    #[test]
    fn catches_flipped_header_byte() {
        let mut buf = forge(FORMAT_MAJOR, FORMAT_MINOR, 0, &[]);
        buf[16] ^= 0xff; // flip the flags byte; CRC no longer matches
        assert!(matches!(
            read_header(&mut Cursor::new(&buf)),
            Err(Error::Corrupt { .. })
        ));
    }

    // #142: a framed length is attacker-controlled and drives the whole reorder
    // decode path. It must be capped and read incrementally, so a tiny truncated
    // frame claiming a huge length errors without allocating the claim.
    #[test]
    fn read_framed_rejects_over_cap_length() {
        let mut buf = Vec::new();
        buf.extend_from_slice(&u32::MAX.to_le_bytes()); // len ~4 GB > MAX_BLOCK_PAYLOAD
        buf.extend_from_slice(&0u32.to_le_bytes()); // crc placeholder
        // No payload. Must reject on the length cap, not resize(4 GB).
        assert!(matches!(
            read_framed(&mut Cursor::new(&buf), "frame"),
            Err(Error::Malformed(_))
        ));
    }

    #[test]
    fn read_framed_truncated_body_is_truncated_not_alloc() {
        // Largest allowed length, but the body is short: `Truncated`, and only the
        // bytes actually present are read (no full-claim zero-fill).
        let mut buf = Vec::new();
        buf.extend_from_slice(&(MAX_BLOCK_PAYLOAD as u32).to_le_bytes());
        buf.extend_from_slice(&0u32.to_le_bytes());
        buf.extend_from_slice(b"only a few bytes");
        assert!(matches!(
            read_framed(&mut Cursor::new(&buf), "frame"),
            Err(Error::Truncated)
        ));
    }

    #[test]
    fn read_framed_roundtrips_a_valid_frame() {
        let payload = b"a legitimate framed payload";
        let mut buf = Vec::new();
        write_framed(&mut buf, payload).unwrap();
        let got = read_framed(&mut Cursor::new(&buf), "frame").unwrap();
        assert_eq!(got, payload);
    }
}
