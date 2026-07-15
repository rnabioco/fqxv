//! The `.fqxv` container: a header followed by independent, parallel-codable
//! blocks.
//!
//! ```text
//! [4] magic "FQXV"
//! [2] format version (LE)
//! [1] sequence context order (k)
//! [1] quality binning tag
//! [1] flags (bit0: '+' normalized; bit1: reordered; bit2: order preserved;
//!            bit3: global-cluster reorder; bit4: names regenerated;
//!            bits5-7: platform tag)
//! [1] group size G (reads interleaved per spot: 1 single-end, 2 paired,
//!                   3-4 single-cell R1/R2/I1[/I2], ...)
//! [4] header_crc (LE) -- CRC-32C over the 10 header-field bytes above, verified
//!     on read so a flipped version/flags/binning-tag/group-size byte is caught
//!     rather than silently changing decode. Present in both layouts.
//! repeated until the terminator:
//!   [8] block payload length (LE, nonzero)
//!   [4] CRC-32C of the payload (LE) -- verified before decode, so corruption is
//!       caught and localized to one block instead of decoded into garbage
//!   [ ] block payload
//! [8] 0  (zero-length terminator block: a streaming, non-seekable decoder
//!         stops here; seekable readers jump to the footer via the trailer)
//! footer (row-group index — lets `inspect`/random access seek, not scan):
//!   [4] n_row_groups (LE)
//!   per row group: [8] byte_offset (LE, points at the group's length field)
//!                  [4] read_count  (LE)
//!   [8] total_reads (LE)
//!   [4] whole_file_crc (LE)  -- CRC-32C of the archive from byte 0 through the
//!       total_reads field; a one-pass end-to-end integrity check (`verify`)
//!   [4] footer_crc (LE)      -- CRC-32C of the footer body above, checked before
//!       any offset in the index is trusted
//! trailer (fixed, at EOF):
//!   [8] footer_offset (LE)   -- seek straight to the footer
//!   [4] magic "FQXF"
//! block payload:
//!   [8] content_digest (LE) -- xxh3-64 of this block's DECODED content (names,
//!       sequence, post-binning quality), verified after decode so a codec bug
//!       that decodes CRC-valid bytes into wrong-but-in-bounds output is caught
//!       at runtime. Distinct from the frame CRC, which only covers stored bytes.
//!       Sits inside the payload, so the frame CRC covers it too.
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

pub(crate) use crate::crc::{crc32c, crc32c_combine, CrcWriter};
pub(crate) use crate::{Error, Result, FORMAT_VERSION, MAGIC};
pub(crate) use fqxv_fqzcomp::QualityBinning;
pub(crate) use std::borrow::Cow;
pub(crate) use std::fs::File;
pub(crate) use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
pub(crate) use xxhash_rust::xxh3::Xxh3;

mod block;
mod compress;
mod decompress;
mod estimate;
mod format;
mod inspect;
mod parse;
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
pub use compress::{compress, compress_auto, compress_interleaved, compress_multi, Params, Stats};
pub use decompress::{content_stats, decompress, decompress_recover, decompress_split, Recovery};
pub use estimate::{estimate, Estimate};
pub use inspect::{inspect, peek, ContentStats, Info, Platform, QUAL_MAX};
pub use verify::{
    expected_reads, verify, verify_quick, verify_report, verify_roundtrip, VerifyCheck,
    VerifyReport,
};

#[cfg(test)]
mod tests;
