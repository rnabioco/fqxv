# fqxv compress

Compress one or more FASTQ files into a single `.fqxv` archive.

## Usage

```bash
fqxv compress <INPUTS>... -o <OUTPUT> [OPTIONS]
```

Give one input for single-end, two for paired-end, or three/four for single-cell
(R1/R2/I1[/I2]). Multiple inputs are interleaved per spot into one archive and
can be split back out with [`decompress --split`](decompress.md). Input order is
preserved for the split.

## Arguments

| Argument | Description |
| --- | --- |
| `<INPUTS>...` | Input FASTQ file(s), plain or gzipped (gzip auto-detected). |

## Options

| Option | Description |
| --- | --- |
| `-o, --output <PATH>` | Output `.fqxv` path (required). |
| `-l, --level <N>` | Effort 1–9; higher raises the sequence context order. Default: 5. |
| `--quality-bin <MODE>` | `lossless` (default), `bin8`, `bin4`, `bin2` (lossy). |
| `--platform <NAME>` | Sequencing platform to record: `illumina`, `nanopore`, `pacbio`, `mgi`. Auto-detected from read names by default; pass to override. |
| `--threads <N>` | Worker threads (0 = all cores). |

## Examples

```bash
# single-end, gzipped input
fqxv compress reads.fastq.gz -o reads.fqxv

# paired-end at higher effort
fqxv compress R1.fq.gz R2.fq.gz -o sample.fqxv --level 7 --threads 16

# 10x single-cell (R1 + R2 + I1 + I2)
fqxv compress R1.fq R2.fq I1.fq I2.fq -o sample.fqxv

# lossy quality (Illumina 8-level binning)
fqxv compress reads.fastq -o reads.fqxv --quality-bin bin8
```

## Notes

- **Losslessness:** read name + description, sequence, and quality are preserved
  exactly. The `+` separator line is normalized to a bare `+`.
- **Grouping:** interleaving mates shrinks the archive — near-identical mate
  names collapse to matches, and a spot's related reads sit together for the
  sequence model.
