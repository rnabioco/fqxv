# fqxv decompress

Restore FASTQ from a `.fqxv` archive — to a file, to split mate files, or
streamed to stdout for a pipe.

## Usage

```bash
fqxv decompress <INPUT> [-o <OUTPUT> | --split <PREFIX>]
```

## Arguments

| Argument | Description |
| --- | --- |
| `<INPUT>` | Input `.fqxv` file. |

## Options

| Option | Description |
| --- | --- |
| `-o, --output <PATH>` | Interleaved FASTQ output. Omit (or use `-`) for stdout. |
| `--split <PREFIX>` | Restore separate per-spot files: `<PREFIX>_1.fastq … _G.fastq`. |
| `--recover` | Best-effort decode of a corrupted archive: skip blocks that fail their CRC and emit the rest. See [below](#recovering-a-corrupted-archive). |
| `--threads <N>` | Worker threads (0 = all cores). |

`--output` and `--split` are mutually exclusive, as are `--split` and `--recover`.

## Examples

```bash
# to a file
fqxv decompress reads.fqxv -o reads.fastq

# stream interleaved to an aligner (no temp files)
fqxv decompress sample.fqxv | bwa mem -p ref.fa -
fqxv decompress sample.fqxv | bowtie2 --interleaved - -x idx

# split a paired / single-cell archive back into its files
fqxv decompress sample.fqxv --split out
#   -> out_1.fastq, out_2.fastq   (paired)
#   -> out_1.fastq ... out_4.fastq (single-cell R1/R2/I1/I2)
```

## Notes

- For a grouped archive, the default (interleaved) output emits `m0, m1, …` per
  spot — exactly the interleaved layout aligners expect with `-p` /
  `--interleaved`.
- `--split` reads the archive's group size from its header and creates that many
  output files, in the original input order.

## Recovering a corrupted archive

If [`fqxv verify`](verify.md) reports an archive as corrupt, `--recover` salvages
everything that is still intact instead of failing outright:

```bash
fqxv decompress --recover damaged.fqxv -o recovered.fastq
```

Because blocks are independent and the footer's row-group index records each
block's absolute offset, a block that fails its CRC (or won't decode) is skipped
by seeking straight to the next one — one bad byte costs a single row group, not
the whole file. Each skipped block is logged with the reads lost, and a summary
is printed at the end:

```text
warning: recovered 812 block(s), skipped 3 corrupt block(s) — 30000 read(s) lost
```

Notes and limits:

- Output is interleaved FASTQ, exactly like a normal `decompress`; `--recover`
  cannot be combined with `--split`.
- Recovery reads the footer index, so it only applies to the **plain** layout.
  A globally-clustered [reordered](../design/container.md#reordered-archives)
  archive is all-or-nothing (its streams are mutually dependent) and returns an
  error directing you to a plain `decompress`.
- If the footer itself is unreadable (e.g. a truncated download), fall back to a
  plain streaming `decompress`, which decodes every whole block before the
  truncation point.
