//! Remote / parallel random access: an IO-free row-group index and per-stream
//! decoders for column projection over an `.fqxv` archive.
//!
//! The container is already Parquet-shaped for this: blocks are independently
//! decodable and individually CRC'd, and the footer row-group index plus the EOF
//! trailer back-pointer let a client fetch the tail, parse the index, and issue
//! parallel range requests for exactly the blocks — or, since v3, the individual
//! streams — it wants. This module exposes that surface without any IO of its own,
//! so a caller drives fetching however it likes (a local `File`, `object_store`,
//! async HTTP, …):
//!
//! - [`Index`] parses the footer, either by seeking a reader ([`Index::read`]) or
//!   from a fetched suffix buffer ([`Index::from_suffix`]).
//! - [`Index::byte_ranges`] turns `(groups, stream)` into the exact byte ranges to
//!   `GET`, and [`Index::verify_stream`] checks a fetched stream against its CRC.
//! - [`decode_names`] / [`decode_sequence`] / [`decode_quality`] decode a single
//!   fetched stream; [`decode_block_contents`] decodes a whole fetched block
//!   payload.
//!
//! Only the plain layout carries this index. The globally-clustered reorder
//! layout is self-describing and has no footer, so [`Index::read`] rejects it
//! (its streams are mutually dependent and cannot be projected).

use super::*;
use std::ops::Range;

/// One of the three per-read streams every block splits into. Their on-disk order
/// within a block payload is names, then sequence, then quality.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Stream {
    /// Read names and descriptions (fqxv-tokenizer). Typically <1% of the archive,
    /// so a client that only needs read IDs can fetch ~100× less than the whole.
    Names,
    /// Sequence bases (fqxv-seq) — the stream a k-mer screen or classifier wants.
    Sequence,
    /// Quality scores (fqxv-fqzcomp).
    Quality,
}

impl Stream {
    /// The three streams in on-disk order.
    pub const ALL: [Stream; 3] = [Stream::Names, Stream::Sequence, Stream::Quality];

    fn idx(self) -> usize {
        match self {
            Stream::Names => 0,
            Stream::Sequence => 1,
            Stream::Quality => 2,
        }
    }
}

/// On-disk location of one row group and its three coded streams.
#[derive(Debug, Clone, Copy)]
pub struct GroupLoc {
    /// Absolute byte offset of the block's frame marker (the start of
    /// `[4 BLOCK_MAGIC][8 len][4 crc][payload]`).
    pub block_offset: u64,
    /// Reads in this group.
    pub read_count: u32,
    streams: [StreamLoc; 3],
}

impl GroupLoc {
    /// Byte range of one coded stream's bytes within the archive — what to `GET`.
    pub fn stream_range(&self, stream: Stream) -> Range<u64> {
        let s = self.streams[stream.idx()];
        s.offset..s.offset + u64::from(s.len)
    }

    /// CRC-32C of one coded stream, for verifying a fetched projection.
    pub fn stream_crc(&self, stream: Stream) -> u32 {
        self.streams[stream.idx()].crc
    }
}

/// A parsed footer row-group index: enough to issue range requests for specific
/// blocks or streams and verify what comes back, with no decoding.
#[derive(Debug, Clone)]
pub struct Index {
    groups: Vec<GroupLoc>,
    total_reads: u64,
    whole_file_crc: u32,
}

/// Result of [`Index::from_suffix`]: either the parsed index, or a request for a
/// longer tail. A client fetches a modest suffix (say `Range: bytes=-65536`),
/// calls `from_suffix`, and on [`SuffixParse::NeedAtLeast`] refetches that many
/// trailing bytes — at most one extra round trip, since the second length is
/// exact.
#[derive(Debug, Clone)]
pub enum SuffixParse {
    /// The suffix reached the footer; here is the parsed index.
    Parsed(Index),
    /// The suffix was too short. Refetch at least this many trailing bytes
    /// (`Range: bytes=-N`) and call [`Index::from_suffix`] again.
    NeedAtLeast(u64),
}

impl Index {
    /// Minimum trailing bytes that could contain the footer trailer. A first
    /// suffix fetch should be at least this large; bigger (e.g. 64 KiB) usually
    /// captures the whole footer in one round trip.
    pub const MIN_SUFFIX: usize = TRAILER_LEN;

    /// Parse the index by seeking the reader to its footer (via the EOF trailer's
    /// back-pointer) and validating the footer CRC. Works on any seekable source —
    /// a local `File`, or an in-memory `Cursor` over a fully-fetched archive.
    ///
    /// Returns [`Error::Malformed`] for the footer-less globally-clustered reorder
    /// layout, whose streams cannot be projected.
    pub fn read<R: Read + Seek>(mut reader: R) -> Result<Index> {
        let header = read_header(&mut reader)?;
        if header.flags & FLAG_GLOBAL_REORDER != 0 {
            return Err(Error::Malformed(
                "reorder layout has no row-group index; random access is unsupported",
            ));
        }
        Index::from_footer(read_footer(&mut reader)?)
    }

    /// Parse the index from a buffer holding the archive's *tail*, IO-free.
    ///
    /// `suffix` must be the final `suffix.len()` bytes of the archive (i.e.
    /// `archive[file_len - suffix.len() .. file_len]`) — the body of a
    /// `Range: bytes=-N` request — and `file_len` the archive's total size. If the
    /// suffix does not reach back to the footer start, returns
    /// [`SuffixParse::NeedAtLeast`] with the exact tail length to refetch.
    pub fn from_suffix(suffix: &[u8], file_len: u64) -> Result<SuffixParse> {
        if suffix.len() as u64 > file_len {
            return Err(Error::Malformed("suffix longer than the archive"));
        }
        if suffix.len() < TRAILER_LEN {
            return Ok(SuffixParse::NeedAtLeast(TRAILER_LEN as u64));
        }
        let trailer = &suffix[suffix.len() - TRAILER_LEN..];
        if trailer[8..] != FOOTER_MAGIC {
            return Err(Error::Malformed("missing footer trailer magic"));
        }
        let footer_offset = u64::from_le_bytes(trailer[..8].try_into().unwrap());
        let body_end = file_len - TRAILER_LEN as u64;
        if footer_offset < HEADER_LEN as u64
            || footer_offset + (4 + 8 + FOOTER_CRC_TAIL as u64) > body_end
        {
            return Err(Error::Malformed("footer offset out of range"));
        }
        // Bytes of the archive this suffix covers begin at this absolute offset.
        let suffix_start = file_len - suffix.len() as u64;
        if footer_offset < suffix_start {
            // Need the tail to reach back to (and include) the footer body.
            return Ok(SuffixParse::NeedAtLeast(file_len - footer_offset));
        }
        let start = (footer_offset - suffix_start) as usize;
        let end = (body_end - suffix_start) as usize;
        let footer = parse_footer_body(&suffix[start..end], footer_offset)?;
        Ok(SuffixParse::Parsed(Index::from_footer(footer)?))
    }

    fn from_footer(footer: Footer) -> Result<Index> {
        let groups = footer
            .groups
            .iter()
            .zip(&footer.stream_locs)
            .map(|(&(block_offset, read_count), &streams)| GroupLoc {
                block_offset,
                read_count,
                streams,
            })
            .collect();
        Ok(Index {
            groups,
            total_reads: footer.total_reads,
            whole_file_crc: footer.whole_file_crc,
        })
    }

    /// The row groups, in archive order.
    pub fn groups(&self) -> &[GroupLoc] {
        &self.groups
    }

    /// Total reads across every group.
    pub fn total_reads(&self) -> u64 {
        self.total_reads
    }

    /// The archive's stored whole-file CRC-32C (a stable on-disk fingerprint).
    pub fn whole_file_crc(&self) -> u32 {
        self.whole_file_crc
    }

    /// Byte ranges to fetch for `stream` across the given row groups, in the order
    /// requested — one `GET` per range, then [`Index::verify_stream`] and a
    /// `decode_*` on each. Errors if any index is out of range.
    pub fn byte_ranges(&self, groups: &[usize], stream: Stream) -> Result<Vec<Range<u64>>> {
        groups
            .iter()
            .map(|&g| {
                self.groups
                    .get(g)
                    .map(|loc| loc.stream_range(stream))
                    .ok_or(Error::Malformed("row-group index out of range"))
            })
            .collect()
    }

    /// Verify a fetched coded stream against the index's CRC-32C, erroring if it
    /// was corrupted in transit or storage. Call before decoding a projection —
    /// the block's joint content digest cannot cover a single fetched stream, so
    /// this is the integrity check for projected reads.
    pub fn verify_stream(&self, group: usize, stream: Stream, coded: &[u8]) -> Result<()> {
        let loc = self
            .groups
            .get(group)
            .ok_or(Error::Malformed("row-group index out of range"))?;
        if crc32c(coded) != loc.stream_crc(stream) {
            return Err(Error::Corrupt {
                what: format!("row group {group} {stream:?} stream"),
            });
        }
        Ok(())
    }
}

/// Decode a projected **names** stream into per-read name/description bytes. Verify
/// it against [`Index::verify_stream`] first.
pub fn decode_names(coded: &[u8]) -> Result<Vec<Vec<u8>>> {
    Ok(fqxv_tokenizer::decode(coded)?)
}

/// Decode a projected **sequence** stream into `(per-read lengths, concatenated
/// bases)`. Slice read `i` out with the lengths' running sum.
pub fn decode_sequence(coded: &[u8]) -> Result<(Vec<u32>, Vec<u8>)> {
    Ok(fqxv_seq::decode(coded)?)
}

/// Decode a projected **quality** stream into `(per-read lengths, concatenated
/// quality)`. The lengths match the sequence stream's.
pub fn decode_quality(coded: &[u8]) -> Result<(Vec<u32>, Vec<u8>)> {
    Ok(fqxv_fqzcomp::decode(coded)?)
}

/// Decoded contents of one block: parallel per-read names and lengths, plus the
/// concatenated sequence and quality (slice each read out with the lengths).
#[derive(Debug, Clone)]
pub struct BlockContents {
    /// Per-read name/description bytes.
    pub names: Vec<Vec<u8>>,
    /// Per-read sequence (and quality) lengths.
    pub lengths: Vec<u32>,
    /// Concatenated sequence bases across the block.
    pub sequence: Vec<u8>,
    /// Concatenated quality scores across the block.
    pub quality: Vec<u8>,
}

/// Decode a whole fetched block **payload** — the bytes after the frame head
/// (`[4 BLOCK_MAGIC][8 len][4 crc]`), i.e. what remains once a fetched block frame
/// has been CRC-checked. Runs the same IO-free decoder and end-to-end content
/// digest check as streaming decompression, for a client that fetched an entire
/// block rather than one stream.
pub fn decode_block_contents(payload: &[u8]) -> Result<BlockContents> {
    let (_n, names, lengths, sequence, quality) = decode_block_parts(payload)?;
    Ok(BlockContents {
        names,
        lengths,
        sequence,
        quality,
    })
}
