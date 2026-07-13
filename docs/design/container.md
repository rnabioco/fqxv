# Container Format

A `.fqxv` file is a small header followed by a sequence of independently-coded
**row groups** (blocks), then a footer index and an EOF trailer. Row groups are
the unit of parallelism and of coarse random access; the footer lets a reader
seek to any of them without scanning the file.

## Layout

```text
[4]  magic "FQXV"
[2]  format version (LE)        (v1)
[1]  sequence context order (k)
[1]  quality binning tag        (0 lossless, 1 bin8, 2 bin4, 3 bin2)
[1]  flags                      (bit0 '+' normalized; bit1 reordered;
                                 bit2 keep-order; bit3 global-reorder)
[1]  group size G               (1 single-end, 2 paired, 3-4 single-cell)
repeated until the terminator:
  [8]  block payload length (LE, nonzero)
  [ ]  block payload
[8]  0                          (zero-length terminator block)
footer:
  [4]  n_row_groups (LE)
  per row group: [8] byte_offset (LE)   (points at the group's length field)
                 [4] read_count  (LE)
  [8]  total_reads (LE)
trailer (at EOF):
  [8]  footer_offset (LE)
  [4]  magic "FQXF"
```

Each block payload:

```text
[4]  n_reads (LE)
[4]  names_len (LE)  [ ] names   (fqxv-tokenizer)
[4]  seq_len   (LE)  [ ] seq     (fqxv-seq)
[4]  qual_len  (LE)  [ ] qual    (fqxv-fqzcomp)
```

## Streaming vs. seeking

The block region carries an inline `[8]` length before every payload, so a
**streaming** decoder reads row groups one after another and stops at the
zero-length terminator — it never needs to seek, so `fqxv decompress` stays
bounded-memory and pipe-friendly (`fqxv decompress x.fqxv | bwa mem -p`). A
**seeking** reader instead jumps to the EOF trailer, follows `footer_offset` back
to the footer, and reads the row-group index directly. The terminator is what
lets the same file serve both: the streaming reader stops before the footer, the
seeking reader skips straight past it.

## Footer index, `inspect`, and random access

The footer records, per row group, the byte offset of its length field and its
read count (read-start is the running sum). It covers the plain layout; reordered
archives use a distinct self-describing layout (below). This buys two things over
the old scan-everything approach:

- **`inspect` is O(row groups), not O(bytes).** It reads the footer for the read
  and block totals, then seeks to each row group and reads only its small block
  header (the `n_reads` and three stream length prefixes), skipping every coded
  payload.
- **`extract a..b` decodes only the row groups overlapping the range.** The
  footer's read counts locate the covering groups; each is decoded and the reads
  inside `[a, b)` are emitted. Granularity is the row group, not the read — every
  codec carries model state across the reads within a group, so a group is the
  smallest independently decodable unit — but only the requested reads are
  written out.

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
`--level`) or a fixed raw-sequence **byte budget**. For fixed short reads
(Illumina) the read count binds and the byte budget never triggers. For long,
ragged reads (nanopore) the byte budget binds first — otherwise a read-count
block of, say, 1M × 10 kb reads would be a ~10 GB row group that destroys
parallelism and random-access granularity and could overflow the `u32` per-stream
compressed length. Byte cuts still land on whole-spot boundaries. Boundaries
depend only on the read lengths and the two limits, never on thread scheduling,
so determinism holds.

## Grouping (paired-end and single-cell)

When `G > 1`, per-spot reads are interleaved (`m0₀, m1₀, …, m0₁, m1₁, …`). Blocks
always contain whole spots and start on member 0, so a block splits back into
its `G` files by local read index mod `G`. Interleaving lets the name tokenizer
collapse the near-identical mate names and keeps a spot's reads adjacent for the
sequence model.

`fqxv decompress` emits interleaved FASTQ (ideal for aligners); `--split`
restores the `G` separate files.

## Reordered archives

`--reorder` uses a distinct whole-file, globally-clustered layout (flag bit3)
that is self-describing and carries no footer or terminator — decode dispatches
on the flag before reading any block, so the footer index applies only to the
plain layout. Both reorder modes share this one path: with `--keep-order`
(flag bit2) names/quality are coded in original order and a permutation restores
it; without it they follow the clustered order and no permutation is written.
`extract` is not supported on reordered archives (their on-disk order is
clustered, not input order); use `fqxv decompress`. See
[Read Reordering](reordering.md).

## Losslessness

Read name + description, sequence, and quality bytes are preserved exactly. The
`+` separator line is normalized to a bare `+` (its optional repeated header is
not retained), matching SPRING and fqz_comp. With `--quality-bin`, quality is
mapped through the chosen Illumina binning table before coding — an explicit,
opt-in lossy transform.
