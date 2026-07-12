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
| `--threads <N>` | Worker threads (0 = all cores). |

`--output` and `--split` are mutually exclusive.

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
