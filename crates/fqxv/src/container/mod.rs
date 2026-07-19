//! The `.fqxv` container: a header followed by independent, parallel-codable
//! blocks.
//!
//! ```text
//! [4] magic "FQXV"
//! [1] format version major (LE) -- a reader refuses any archive whose major
//!     differs from its own FORMAT_MAJOR; the archive is wire-incompatible.
//! [1] format version minor (LE) -- backward-compatible additions; a reader
//!     tolerates any minor within its major (informational).
//! [8] required_features (LE u64) -- coarse capability bits (fqxv::feature). A set
//!     bit outside the reader's KNOWN_FEATURES is refused (UnsupportedFeature)
//!     rather than mis-decoded. Per-block codec choices are NOT gated here — they
//!     ride the sequence stream's method byte (UnsupportedMethod on decode).
//! [1] sequence context order (k)
//! [1] quality binning tag
//! [1] flags (bit0: '+' normalized; bit1: reordered; bit2: order preserved;
//!            bit3: global-cluster reorder; bit4: names regenerated;
//!            bit5: global reference frame present; bits6-7: free)
//! [1] group size G (reads interleaved per spot: 1 single-end, 2 paired,
//!                   3-4 single-cell R1/R2/I1[/I2], ...)
//! [1] platform tag (0 unknown, 1 Illumina, 2 Nanopore, 3 PacBio, 4 MGI/BGI).
//!     Its own byte: it previously shared flags bits5-7, where Illumina's code
//!     collided with the global-reference bit and made Illumina reorder archives
//!     undecodable.
//! [2] ext_len (LE u16) -- bytes of the extension region that follows.
//! [ext_len] extension records -- each [1 tag][2 len LE][len bytes]. Tag high bit
//!     marks a *critical* record: a reader that doesn't know a critical tag refuses
//!     the archive (UnsupportedExtension); an unknown non-critical tag is skipped.
//!     Empty at 1.0 -- the region lets a later minor add skippable header fields
//!     without a major bump. Covered by the header CRC.
//! [4] header_crc (LE) -- CRC-32C over every header byte above (prefix + extension
//!     region), verified on read so a flipped version/features/flags/binning-tag/
//!     group-size/platform byte is caught rather than silently changing decode.
//!     Present in both layouts. The block region (or the reference frame below, if
//!     present) begins right after this; seek/scan start offsets use the actual
//!     header length (prefix + ext_len + 4), which equals the ext-empty HEADER_LEN
//!     for a 1.0-written archive.
//! optional whole-file reference frame (plain layout only, present iff the flags
//!   bit5 FLAG_GLOBAL_REFERENCE and the GLOBAL_REFERENCE feature bit are set):
//!   [4] len (LE)  [4] CRC-32C  [len] reference bytes -- a single framed
//!       `fqxv_lroverlap::Reference` (consensus contigs, entropy-coded), assembled
//!       once over the whole file. Long-read sequence blocks then code against it
//!       (sequence method byte 2) instead of re-storing a reference per block, so
//!       the genome is stored once, not once per 256 MiB block. The footer's block
//!       offsets start past this frame. A reader without the GLOBAL_REFERENCE
//!       feature refuses the archive (upgrade signal); a block referencing a
//!       missing frame fails closed.
//! repeated until the terminator:
//!   [4] BLOCK_MAGIC "FQXB" -- per-block sync marker; recovery scans for it to
//!       resynchronize to a block boundary when the footer is lost or a length
//!       prefix is corrupt (see `decompress_recover`)
//!   [8] block payload length (LE, nonzero)
//!   [4] CRC-32C of the payload (LE) -- verified before decode, so corruption is
//!       caught and localized to one block instead of decoded into garbage
//!   [ ] block payload
//! [4] BLOCK_MAGIC  [8] 0  (zero-length terminator frame: a streaming, non-seekable
//!         decoder stops here; seekable readers jump to the footer via the trailer)
//! footer (row-group index — lets `inspect`/random access seek, not scan):
//!   [4] n_row_groups (LE)
//!   per row group: [8] byte_offset (LE, points at the group's frame marker)
//!                  [4] read_count  (LE)
//!                  per stream (names, sequence, quality, in that order):
//!                    [8] stream_offset (LE, absolute offset of the coded bytes)
//!                    [4] stream_len    (LE, coded length)
//!                    [4] stream_crc32c (LE, over exactly those coded bytes)
//!     The per-stream triple lets a remote client project one column — fetch just
//!     the names (~1% of the archive) or just the sequence — with a single range
//!     request and verify it against `stream_crc32c` (the block content digests
//!     cover the DECODED streams, so they can't check a projected fetch of coded
//!     bytes in isolation).
//!   [8] total_reads (LE)
//!   [4] whole_file_crc (LE)  -- CRC-32C of the archive from byte 0 through the
//!       total_reads field; a one-pass end-to-end integrity check (`verify`)
//!   [4] footer_crc (LE)      -- CRC-32C of the footer body above, checked before
//!       any offset in the index is trusted
//! trailer (fixed, at EOF):
//!   [8] footer_offset (LE)   -- seek straight to the footer
//!   [4] magic "FQXF"
//! block payload:
//!   [8] names_digest (LE) -- xxh3-64 of this block's DECODED names
//!   [8] seq_digest   (LE) -- xxh3-64 of this block's DECODED sequence
//!   [8] qual_digest  (LE) -- xxh3-64 of this block's DECODED post-binning quality
//!       One digest per stream, each verified after decode so a codec bug that
//!       decodes CRC-valid bytes into wrong-but-in-bounds output is caught at
//!       runtime AND localized to the offending stream. Distinct from the frame
//!       CRC, which only covers stored bytes. Each folds in n_reads + the stream's
//!       per-read lengths so no byte can slide across a read or stream boundary.
//!       Sit inside the payload, so the frame CRC covers them too.
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
//! `[8] n  [ ] flip  [ ] perm  [ ] template  [4] n_blocks  [seq block]*n
//!  [ [names][qual] ]*n  [ ] output_digest`
//! (each `[ ]` is a `[u32 len][u32 crc32c][bytes]` frame, CRC-verified on decode;
//! `perm` is empty without keep-order, `template` is empty unless regenerating
//! names). The trailing `output_digest` frame holds an xxh3-64 over the reads in
//! output order (the reorder analog of the per-block content digest), verified
//! after decode so a codec bug that reconstructs wrong reads is caught.
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

pub(crate) use crate::crc::{CrcWriter, crc32c, crc32c_combine};
pub(crate) use crate::{Error, FORMAT_MAJOR, FORMAT_MINOR, KNOWN_FEATURES, MAGIC, Result};
pub(crate) use fqxv_fqzcomp::QualityBinning;
pub(crate) use std::borrow::Cow;
pub(crate) use std::fs::File;
pub(crate) use std::io::{self, BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
pub(crate) use xxhash_rust::xxh3::Xxh3;

mod block;
mod compress;
mod decompress;
mod estimate;
mod format;
mod inspect;
mod parse;
mod random_access;
mod records;
mod reorder;
mod verify;

// Internal flat namespace: every submodule item is `pub(crate)`, re-globbed here
// so each submodule's `use super::*` sees the whole set.
pub(crate) use block::*;
pub(crate) use compress::*;
pub(crate) use format::*;
pub(crate) use inspect::*;
pub(crate) use parse::*;
pub(crate) use reorder::*;
// `decompress` and `verify` expose only their public surface (below) to the rest
// of the crate; their internals stay module-private, so no `*` glob here.

// Public surface, re-exported unchanged from `lib.rs` at the crate root. These
// explicit `pub use`s elevate the names above the `pub(crate)` globs above so the
// external API is preserved.
pub use compress::{Params, Stats, compress, compress_auto, compress_interleaved, compress_multi};
pub use decompress::{Recovery, content_stats, decompress, decompress_recover, decompress_split};
pub use estimate::{Estimate, estimate};
pub use inspect::{ContentStats, Info, Platform, QUAL_MAX, inspect, peek};
pub use random_access::{
    BlockContents, GroupLoc, Index, Stream, SuffixParse, decode_block_contents, decode_names,
    decode_quality, decode_quality_with_seq, decode_sequence, quality_needs_sequence,
};
pub use records::{Record, RecordReader, decompress_records};
pub use verify::{
    VerifyCheck, VerifyReport, expected_reads, verify, verify_quick, verify_report,
    verify_roundtrip,
};

#[cfg(test)]
mod tests;
