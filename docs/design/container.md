# Container Format

A `.fqxv` file is a small header followed by a sequence of independently-coded
**row groups** (blocks), then a footer index and an EOF trailer. Row groups are
the unit of parallelism and of coarse random access; the footer lets a reader
seek to any of them without scanning the file.

## Layout

```text
[4]  magic "FQXV"
[1]  format version major (LE)  (reader refuses a differing major)
[1]  format version minor (LE)  (tolerated within a major; informational)
[8]  required_features (LE u64) (coarse capability bits; an unknown set bit is
                                 refused rather than mis-decoded)
[1]  sequence context order (k)
[1]  quality binning tag        (0 lossless, 1 bin8, 2 bin4, 3 bin2,
                                 4 ont, 5 hifi)
[1]  flags                      (bit0 '+' normalized; bit1 reordered;
                                 bit2 keep-order; bit3 global-reorder;
                                 bit4 regen-names; bit5 global-reference;
                                 bits6-7 free. An unknown set bit is refused,
                                 not ignored — every flag changes how the
                                 archive is read)
[1]  group size G               (1 single-end, 2 paired, 3-4 single-cell)
[1]  platform tag               (0 unknown, 1 Illumina, 2 Nanopore,
                                 3 PacBio, 4 MGI/BGI)
[2]  ext_len (LE u16)           (bytes of the extension region below)
[ ]  extension region           (TLV records [1 tag][2 len][len bytes]; the tag's
                                 high bit marks a record as critical — an unknown
                                 critical tag is refused, an unknown non-critical
                                 tag is skipped. Tag 0x01 (non-critical) carries
                                 the per-member slot labels)
[4]  header CRC-32C (LE)        (over the header prefix + extension region)
optional whole-file reference frame (plain layout only; present iff flag bit5 and
the GLOBAL_REFERENCE feature bit are set):
  [4]  len (LE)  [4] CRC-32C  [len] reference bytes
                                (one framed long-read consensus reference,
                                 assembled once over the whole file; the footer's
                                 block offsets start past it)
repeated until the terminator:
  [4]  block sync marker "FQXB"          (recovery scans for it to resync)
  [8]  block payload length (LE, nonzero)
  [4]  CRC-32C of the payload (LE)       (verified before the block is decoded)
  [ ]  block payload
[4]  "FQXB" [8] 0               (zero-length terminator frame)
footer:
  [4]  n_row_groups (LE)
  per row group: [8] byte_offset (LE)   (points at the group's frame marker)
                 [4] read_count  (LE)
                 per stream (names, sequence, quality, in that order):
                   [8] stream_offset (LE)   (absolute offset of the coded bytes)
                   [4] stream_len    (LE)   (coded length)
                   [4] stream_crc32c (LE)   (over exactly those coded bytes)
  [8]  total_reads (LE)
  [4]  whole_file_crc (LE)              (CRC-32C of byte 0 .. total_reads)
  [4]  footer_crc (LE)                  (CRC-32C of the footer body above)
trailer (at EOF):
  [8]  footer_offset (LE)
  [4]  magic "FQXF"
```

The on-disk format is versioned `major.minor` (currently **1.0**, independent of
the crate version). A reader refuses an archive whose **major** differs from its
own — that archive is wire-incompatible — but tolerates a newer **minor**, since a
minor bump only adds backward-compatible things (a skippable header extension
record, or a new optional codec gated by a `required_features` bit). See
[Versioning and evolution policy](#versioning-and-evolution-policy) for what
warrants which.

**The header is variable-length.** The fixed prefix is 21 bytes and the trailing
CRC is 4, but the extension region between them is *not* always empty — a 1.0
archive records per-member slot labels there (tag `0x01`), so its header is longer
than the 25-byte extension-empty minimum. A reader MUST parse the extension region
(regardless of the minor version — 1.0 archives already use it) and MUST derive the
header length as `prefix + ext_len + 4` from the archive's own `ext_len`. Anything
positioned after the header — the first block frame, a recovery scan's start, the
block offsets recorded in the footer — is relative to *that* length, never a
hardcoded 25.

Each block payload:

```text
[8]  names_digest (LE)                 (xxh3-64 of the block's decoded names)
[8]  seq_digest   (LE)                 (xxh3-64 of the block's decoded sequence)
[8]  qual_digest  (LE)                 (xxh3-64 of the block's decoded quality)
[4]  n_reads (LE)
[4]  names_len (LE)  [ ] names   (fqxv-tokenizer)
[4]  seq_len   (LE)  [ ] seq     (fqxv-seq)
[4]  qual_len  (LE)  [ ] qual    (fqxv-fqzcomp)
```

The three streams are contiguous within the payload, so each one's absolute
`(offset, len, crc32c)` is recorded in the footer index (above). That is what
lets a reader fetch a single **column** — just the names, or just the
sequence — without reading the whole block; see
[Column projection](#footer-index-info-and-random-access). The three
`*_digest` fields are an end-to-end round-trip check over the block's *decoded*
content — one digest per stream (names, sequence, post-binning quality), verified
after decode, so a mismatch localizes which stream a codec round-tripped wrong —
see [Integrity](#integrity).

The stream count is a constant, not a field on disk, so a payload is **exactly**
three streams: a reader that finds bytes left over after the quality stream must
reject the block rather than decode the three it recognizes. Otherwise a payload
carrying a fourth stream would decode and pass every CRC and every content digest
while silently discarding data. Adding a stream within format major 1 is not
permitted; see [the gating rule](#versioning-and-evolution-policy).

### Sequence method bytes

The `seq` stream leads with one byte naming its codec, so a block is
self-describing and one file may mix codecs block to block:

| byte | codec |
| --- | --- |
| `0` | order-k context model (`fqxv-seq`) — always coded, the floor |
| `1` | long-read overlap-consensus, block-local (`fqxv-lroverlap`) |
| `2` | long-read overlap against the archive's shared reference frame |
| `3` | raw large-window LZMA over the block's bases (`fqxv-seq`'s clean-room LZMA) |
| `4` | multi-reference tiling against earlier raw reads (Nanopore) |

The encoder codes several candidates and keeps the smallest, so the method byte
records a measured outcome rather than a mode the user selected; order-k is always
among them, which is what makes a long-read codec unable to enlarge a block. An
unrecognized value is refused (`UnsupportedMethod`) rather than guessed at, and
method `2` fails closed when the archive carries no reference frame. See
[long reads](longread.md#wiring-and-the-per-block-coverage-cap).

## Streaming vs. seeking

The block region carries an inline `[8]` length before every payload, so a
**streaming** decoder reads row groups one after another and stops at the
zero-length terminator — it never needs to seek, so `fqxv decompress` stays
bounded-memory and pipe-friendly (`fqxv decompress x.fqxv -Z | bwa mem -p`). A
**seeking** reader instead jumps to the EOF trailer, follows `footer_offset` back
to the footer, and reads the row-group index directly. The terminator is what
lets the same file serve both: the streaming reader stops before the footer, the
seeking reader skips straight past it.

## Footer index, `info`, and random access

The footer records, per row group, its frame-marker byte offset, its read count
(read-start is the running sum), and an `(offset, len, crc32c)` triple for each
of the three coded streams (names, sequence, quality). It covers
the plain layout; reordered archives use a distinct self-describing layout
(below). This buys three things over a scan-everything approach:

- **`info` is O(row groups), not O(bytes).** The per-stream sizes live in the
  footer itself, so `info` reads the footer once and sums them — no per-block
  seeks at all.
- **Coarse random access.** The per-group `read_count` (read-start is the running
  sum) locates the row groups overlapping any read range, so a reader can seek to
  and decode just those. Granularity is the row group, not the read — every codec
  carries model state across the reads within a group, so a group is the smallest
  independently decodable unit.
- **Column projection.** Because the footer pins each stream's absolute
  `(offset, len)`, a reader can fetch a **single column** across some or all row
  groups — just the read names, or just the sequence — with one range request per
  group and never touch the others. Read names are typically <1% of the archive,
  so a client that only wants IDs (counting reads, building an external index)
  fetches ~100× less than the whole file; a sequence-only client (k-mer
  screening, taxonomic classification, alignment-free QC) skips the quality
  stream entirely. The per-stream `crc32c` makes the projection verifiable: the
  block's content digests cover the *decoded* streams and so cannot check a
  projected fetch of *coded* bytes, so each coded stream also carries its own
  CRC-32C in the footer.

The per-group record size is a compile-time constant recorded nowhere on disk, so
a reader checks the footer body's length **exactly** against its `n_row_groups` —
`4 + 60 × n_row_groups + 8 + 4` bytes — and rejects any other size. Merely
requiring "large enough" would let a footer written at a different stride pass its
own CRC and then be walked at the wrong one, producing offsets that survive the
range checks and get reported as fact. The consequence is deliberate: **the footer
cannot grow within format major 1**, and any change to its shape must set a
`required_features` bit so an older reader refuses at the header, long before it
seeks here.

### Random-access API

This is exposed IO-free so a caller can drive fetching however it likes (a local
`File`, `object_store`, async HTTP over `Range` requests). The `fqxv` library
provides an `Index` that parses from a seekable reader (`Index::read`) or from a
fetched archive **tail** (`Index::from_suffix`, which reports the exact longer
tail to refetch when the first suffix fell short of the footer). From it,
`Index::byte_ranges(groups, stream)` yields the byte ranges to fetch,
`Index::verify_stream` checks a fetched column against its CRC, and
`decode_names` / `decode_sequence` / `decode_quality` (or `decode_block_contents`
for a whole fetched block) turn the fetched bytes back into reads. The reordered
layout has no footer index, so `Index::read` rejects it — its streams are
mutually dependent and cannot be projected.

The Python `fqxv.remote` module surfaces this projection over HTTP `Range`
requests: `RemoteArchive` fetches the footer tail with `parse_index_suffix`, then
reads one column with `Index.stream_range` / `verify_stream` and `decode_*_bytes`
(a custom async client can drive those primitives directly for concurrent fetches).
Whole-archive **streaming** decode needs no index at all — it reads the blocks
front-to-back and stops at the terminator frame — so it works over any forward byte
stream: `fqxv decompress -` (pipe in `aws s3 cp s3://… -`, `curl`, …) on the CLI,
and `fqxv.open(file_like)` / `fqxv.remote.stream(url)` in Python.

## Blocks and parallelism

The compressor reads FASTQ into row groups and hands a batch to `rayon` to
compress concurrently, writing them back in order and recording each in the
footer index. Decompression reads length-prefixed blocks and decodes a batch in
parallel. Because each block is self-contained (its own name/sequence/quality
models), no block depends on another — so the archive, **including the footer
offsets**, is byte-identical regardless of thread count, and a block can be
decoded without touching the rest of the file.

## Row-group sizing (short and long reads)

A row group is cut at whichever comes first: `block_reads` reads (set by
`--level`) or a raw-sequence **byte budget**. For fixed short reads
(Illumina) the read count binds and the byte budget never triggers. For long,
ragged reads (nanopore) the byte budget binds first — otherwise a read-count
block of, say, 1M × 10 kb reads would be a ~10 GB row group that destroys
parallelism and random-access granularity and could overflow the `u32` per-stream
compressed length. Byte cuts still land on whole-spot boundaries. Boundaries
depend only on the read lengths and the two limits, never on thread scheduling,
so determinism holds.

The byte budget resolves per platform (`Params::block_seq_bytes`, 0 = auto):
**64 MiB for Nanopore, 256 MiB otherwise**. Because the byte budget is the sole
cut for long reads, it alone sets the archive's block count — the unit of decode
parallelism — and 256 MiB left a typical single-flowcell ONT file at ~2 blocks
with a decode curve flat in the thread count (issue #273). 64 MiB is the
measured knee: +1.08% archive for 3.2× full decode at 8 threads on a 576 MB
MinION file, the cost concentrated in the overlap codec's per-block coverage.
The CLI's `--max` pins the cap back (smallest-archive contract), and an explicit
`--block-reads` also pins it so the read count alone governs. The budget is a
pure function of input and parameters — platform comes from the caller or
content/name detection, never the pool — so thread-count determinism survives.

Capping a block at the byte budget also caps how much coverage a long-read block's
overlap codec sees, and each block otherwise self-assembles (and re-stores) its own
consensus reference. For long-read input the container hoists that reference out
(explicitly-declared Nanopore excepted — at ~10% error the whole-file consensus is
too noisy to pay, so it stays on the per-block layout and streams):
it assembles one consensus over the whole file and stores it **once** in a framed
region between the header and the first block (gated by the `GLOBAL_REFERENCE`
feature bit and flag bit5), then codes every block's reads against that frozen
frame (sequence method byte 2). Placement is per-read against an immutable frame,
so blocks stay byte-budgeted, parallel, and independently decodable given the
frame, while the genome is stored once rather than once per block. A whole-file
never-worse gate adopts the layout only when the frame plus the reference-coded
sequence beats the plain layout it would otherwise fall back to (which is itself
the per-block best of the other methods, not order-k alone). See
[long reads](longread.md#wiring-and-the-per-block-coverage-cap).

## Grouping (paired-end and single-cell)

When `G > 1`, per-spot reads are interleaved (`m0₀, m1₀, …, m0₁, m1₁, …`). Blocks
always contain whole spots and start on member 0, so a block splits back into
its `G` files by local read index mod `G`. Interleaving lets the name tokenizer
collapse the near-identical mate names and keeps a spot's reads adjacent for the
sequence model.

`fqxv decompress` emits interleaved FASTQ (ideal for aligners); `--split`
restores the `G` separate files.

## Reordered archives

`--order any` uses a distinct whole-file, globally-clustered layout (flag bit3,
global-reorder) that is self-describing and carries no footer or terminator —
decode dispatches on the flag before reading any block, so the footer index
applies only to the plain layout. Both reorder modes share this one path: with
keep-order (flag bit2) names/quality are coded in original order and a
permutation restores it; without it they follow the clustered order and no
permutation is written. Grouped (paired / single-cell) input reorders too: the
reads are clustered ignoring mate structure, but the group size is recorded in
the header and the permutation reconstructs the original spot interleaving, so
keep-order is forced on and the archive de-interleaves cleanly on `--split`. See
[Read Reordering](reordering.md).

After the shared header, the reorder layout is a run of length-prefixed,
CRC-guarded frames (each `[4] len · [4] CRC-32C · payload`):

```text
[8]  n_reads (LE)
     flip bitmap        (one bit per clustered read: stored reverse-complemented)
     permutation        (byte-plane split → rANS; empty when order is not kept)
     name template      (counter template for regenerated names; empty otherwise)
     global reference   (present only when flag bit5, global-reference, is set)
[4]  n_blocks (LE)
     per block: sequence payload
     per block: names payload, quality payload
     output digest      (xxh3-64 over the reads in emit order)
```

**Sequence codecs (chosen per block, smaller wins).** Each sequence block is
differential-coded against the clustered order and tagged with a leading version
byte, so blocks may mix versions and decode dispatches on the byte:

- **v2** — single-contig clustered coding (block-local).
- **v3** — adds literal-rescue (block-local): reads the single-contig codec would
  strand as literals are recovered against a block-local reference.
- **v4** — codes reads as `(contig, offset, mismatches)` positions on **one frozen
  whole-file global reference** assembled (SPRING-style) over every clustered read,
  so cross-block overlaps v3 strands as literals collapse to a cheap
  back-reference. The reference is assembled once, then **overlap-merged** (contigs
  whose suffix overlaps another's prefix are chained into fewer, longer
  super-contigs) and stored once in the `global reference` frame.

Exactly one block-local codec is coded per archive: v3 on the adaptive rescue path
(the default under `--order any`) and v2 under `--no-rescue`, since v3 generalizes
v2 and degenerates to the same coding, making it the block-local floor on its own.
v4 is a candidate only on the rescue path, coded per block alongside the
block-local one; the smaller wins, and the whole-file layout is adopted only when
`reference frame + Σ min(block-local, v4)` is strictly smaller than `Σ block-local`.
So the flag bit5 reference is written only when it nets a whole-file win — v4 can
never enlarge the archive. Ties keep the lower version, for determinism.

**Reference-frame coding.** The `global reference` frame begins with its own method
byte. Every method is clean-room and pure-Rust — there is no external compressor
here (the earlier xz/`liblzma` reference path was removed):

- `0` — the whole reference in a single order-k `fqxv-seq` pass. Decode-only: no
  current writer emits it, but archives that carry it still read.
- `1` — the reference split into a fixed 64 contig blocks, each coded with order-k
  `fqxv-seq` **in parallel**. The block count is fixed (never derived from
  `--threads`), so the coded bytes are identical regardless of thread count.
- `4` — SPRING's approach: 2-bit-pack the ACGT consensus (a hard 2 bits/base floor,
  non-ACGT bytes held in an exception list) and run the clean-room LZMA over the
  packed bytes, which reaches the long-range near-duplicate-contig repeats the
  context model cannot see. Usually the winner on a real reference.

The encoder codes `4` and `1` and keeps the smaller, so the method byte here is
also a measured outcome. Method tags `2` (LZ77) and `3` (BWT) were prototypes that
lost on the raw bases and were removed.

## Integrity

Every archive carries CRC-32C (Castagnoli) checksums — the same checksum family
BGZF/BAM and CRAM use — layered so a single flipped bit is both *detected* and
*localized* rather than silently decoded into wrong bases or quality scores.
Four checksums sit at four levels of the layout above:

- **Per-block payload CRC.** Each block frame is `[8] length · [4] CRC-32C · [ ]
  payload`; the CRC is verified *before* the payload is handed to the codecs, so
  corruption is caught and confined to one row group. This is what makes
  [`decompress --recover`](../cli/decompress.md#recovering-a-corrupted-archive)
  possible — a bad block is identified and skipped, not decoded into garbage.
- **Per-stream CRC.** The footer's per-group `(offset, len, crc32c)` triples carry
  a CRC-32C over each coded stream's bytes, so a projected single-column fetch
  ([Column projection](#footer-index-info-and-random-access)) stays verifiable —
  the per-block payload CRC covers all three coded streams at once, and the
  content digests cover the *decoded* streams, so neither can check a projected
  fetch of one coded stream in isolation.
- **Footer CRC** (`footer_crc`) guards the row-group index. It covers the footer
  body and is checked before any byte offset in the index is trusted, so a
  damaged footer can't send a seeking reader to a bogus location.
- **Whole-file CRC** (`whole_file_crc`) is a CRC-32C over the archive from byte 0
  through the `total_reads` field — a single end-to-end digest of everything but
  the two footer CRCs and the trailer.

A fifth check is not a CRC: each block payload leads with three xxh3-64 digests
over its *decoded* content — one each for names, sequence, and post-binning
quality — verified after decode. Where the CRCs catch corruption of the *stored*
bytes, these digests catch a codec that turned CRC-valid bytes into
wrong-but-in-bounds output, and localize the failure to the offending stream.

[`fqxv verify`](../cli/verify.md) runs without any codec, so it is far cheaper
than a full `decompress`: it validates the header CRC, validates `footer_crc`
before trusting an offset, and re-hashes the archive prefix against
`whole_file_crc`. That last digest covers every block payload byte, so it subsumes
the per-block CRCs rather than repeating them; a whole-file failure then triggers a
block-by-block scan to name the damaged row groups. `--quick` inverts the trade —
it walks the footer index and checks each block's stored CRC with parallel
positioned reads, skipping the whole-file digest (and so the bytes no block
covers). The per-stream CRCs are checked on demand by the projection path rather
than by either mode. The globally-clustered
[reordered layout](#reordered-archives) carries no footer, so `verify` there
drives every frame CRC by decoding into a sink instead.

## Losslessness

Read name + description, sequence, and quality bytes are preserved exactly. The
`+` separator line is normalized to a bare `+` (its optional repeated header is
not retained), matching SPRING and fqz_comp. With `--quality-bin`, quality is
mapped through the chosen binning table (Illumina, ONT, or PacBio HiFi) before
coding — an explicit, opt-in lossy transform.

## Versioning and evolution policy

The format version (`major.minor`, currently 1.0) is **independent of the crate
version** — the format is at 1.0 while the crates are still 0.x. The format
version changes only when the bytes on disk do.

The format has four escalating ways to add something, and the choice is decided by
one question: *what must an old reader do when it meets this?*

- **Minor bump alone** — the old reader may ignore it and still decode every byte
  correctly. This covers a new **non-critical** extension record and nothing else.
  Old readers skip the record; the minor is informational (nothing branches on it).
  The per-member slot labels (tag `0x01`) are the worked example: they change only
  what `decompress_split` *names* its outputs, so an archive that carries them
  decodes to identical reads on a reader that has never heard of them.
- **A `required_features` bit** — the old reader cannot decode the archive at all
  without a capability, and that is knowable before the blocks are written. The bit
  is refused at header-read (`UnsupportedFeature`, "upgrade fqxv"), before a single
  block is touched. This is the gate for whole-archive structural change; the
  whole-file `GLOBAL_REFERENCE` frame is the existing one. (An unknown *flags* bit
  is refused the same way, as `UnsupportedFlags` — every flag changes how the
  archive is read, so a bit a reader does not act on must not be ignored — but the
  flags byte has only two free bits, so it is not an evolution mechanism. A new flag
  should normally come with a feature bit; refusing unknown bits makes forgetting
  that fail closed rather than decode under the old interpretation.)
- **A critical extension tag** (tag high bit set) — the same "refuse it" outcome
  (`UnsupportedExtension`) for a *header field* whose meaning an old reader must not
  guess. Prefer a feature bit for capabilities; reserve a critical tag for metadata
  that is genuinely header-shaped and load-bearing.
- **A major bump** — the change cannot be expressed by any of the above: the fixed
  header prefix, the frame/trailer framing, or the magic itself changes, so an old
  reader cannot even parse far enough to be told no.

**The gating rule.** Any change to the footer's shape — the field set or order of
the row-group index — or to the block payload's stream layout MUST be gated by a
`required_features` bit, so an old reader refuses at the header instead of parsing
a differently-shaped structure into plausible garbage. A *per-block codec* choice
is the exception, and only because it is already self-describing: it rides the
sequence stream's leading method byte and surfaces as `UnsupportedMethod`
mid-decode. That fails loudly, but it fails later than a header refusal — so a new
method byte is acceptable for an alternative encoding of the same streams, and is
*not* a substitute for a feature bit when the block's layout itself changes.

**The rule is enforced, not merely stated.** A structure that grew without its
feature bit would otherwise still pass every checksum, because the CRCs and content
digests only attest to what a reader chose to look at. So the reader closes both
implicit extension points itself:

- a footer body whose length is not exactly `4 + 60 × n_row_groups + 8 + 4` is
  rejected, rather than walked at a stride it guessed;
- a block payload with bytes left over after its three streams is rejected, rather
  than decoded as three streams and silently shortened.

Together with the unknown-flag refusal above, that is three places where a
plausible-looking newer archive now stops at a clear error instead of being
mis-read as valid.

**On major bumps.** A major bump is a last resort. The extension mechanisms are
sized so it should not be needed: 63 unused `required_features` bits, a TLV
extension region that can hold up to 64 KiB of skippable or critical header
records, and 251 unused sequence method-byte values. Between them, foreseeable
evolution — new codecs, new metadata, new whole-file structures — fits inside
format 1. Today a differing major is refused unconditionally, with no legacy read
path; if a major bump ever does happen, the reader will keep a read path for
format 1 archives, so an archive written today does not become unreadable.
