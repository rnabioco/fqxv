# Container Format

A `.fqxv` file is a small header followed by a sequence of independently-coded
blocks. Blocks are the unit of parallelism and of coarse random access.

## Layout

```text
[4]  magic "FQXV"
[2]  format version (LE)
[1]  sequence context order (k)
[1]  quality binning tag        (0 lossless, 1 bin8, 2 bin4, 3 bin2)
[1]  flags                      (bit0: '+' line normalized)
[1]  group size G               (1 single-end, 2 paired, 3-4 single-cell)
repeated until EOF:
  [8]  block payload length (LE)
  [ ]  block payload
```

Each block payload holds up to 256K reads:

```text
[4]  n_reads (LE)
[4]  names_len (LE)  [ ] names   (fqxv-tokenizer)
[4]  seq_len   (LE)  [ ] seq     (fqxv-seq)
[4]  qual_len  (LE)  [ ] qual    (fqxv-fqzcomp)
```

## Blocks and parallelism

The compressor reads FASTQ into blocks and hands a batch of blocks to `rayon` to
compress concurrently, writing them back in order. Decompression reads
length-prefixed blocks and decodes a batch in parallel. Because each block is
self-contained (its own name/sequence/quality models), no block depends on
another — so the archive is deterministic regardless of thread count, and a
block can be decoded without touching the rest of the file.

## Grouping (paired-end and single-cell)

When `G > 1`, per-spot reads are interleaved (`m0₀, m1₀, …, m0₁, m1₁, …`). Blocks
always contain whole spots and start on member 0, so a block splits back into
its `G` files by local read index mod `G`. Interleaving lets the name tokenizer
collapse the near-identical mate names and keeps a spot's reads adjacent for the
sequence model.

`fqxv decompress` emits interleaved FASTQ (ideal for aligners); `--split`
restores the `G` separate files.

## Losslessness

Read name + description, sequence, and quality bytes are preserved exactly. The
`+` separator line is normalized to a bare `+` (its optional repeated header is
not retained), matching SPRING and fqz_comp. With `--quality-bin`, quality is
mapped through the chosen Illumina binning table before coding — an explicit,
opt-in lossy transform.
