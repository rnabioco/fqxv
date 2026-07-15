# fqxv compress

Compress one or more FASTQ files into a single `.fqxv` archive.

## Usage

```bash
fqxv compress <INPUTS>... [-o <OUTPUT>] [OPTIONS]
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
| `-o, --output <PATH>` | Output `.fqxv` path. Defaults to the first input's name with the FASTQ/gzip extension replaced by `.fqxv` (`reads.fastq.gz` → `reads.fqxv`), alongside the input. Required when reading from stdin (`-`). |
| `--max` | Maximum-compression preset: deepest sequence context plus read reordering *where it helps* (applied to short reads, auto-skipped for long reads). Overrides `--level`/`--order`. |
| `--estimate` | Predict the archive size and compression ratio from a sample of the input, then exit **without writing anything**. See [Estimating compression](#estimating-compression). |
| `--threads <N>` | Worker threads (0 = all cores). |

### Advanced options

| Option | Description |
| --- | --- |
| `-l, --level <N>` | Effort 1–9; higher raises the sequence context order and block size. Default: 5. |
| `--order <MODE>` | Read-order guarantee: `preserve` (default, restores original order), `any` (allows reordering for a better ratio; single-end order may change), or `shuffle` (like `any`, but discards order and regenerates purely positional names — reorder-lossy, single-end only). |
| `--interleaved <N>` | Interleaving of a *single* input, in members per spot (1 = single-end, 2 = paired as from `sracha get -Z`). Auto-detected from read names by default. Ignored with multiple inputs. |
| `--keep-order` | With `--order any`, force original read order to be restored (store a permutation, code names/quality in original order). Chosen automatically when it makes the archive smaller. |
| `--no-rescue` | With `--order any`, disable the adaptive assembly codecs (block-local literal-rescue and the whole-file global reference) and use the faster single-contig sequence codec only. |
| `--quality-bin <MODE>` | `lossless` (default), `bin8`, `bin4`, `bin2` (lossy). |
| `--platform <NAME>` | Sequencing platform to record: `illumina`, `nanopore`, `pacbio`, `mgi`. Auto-detected from read names by default; pass to override. |

## Examples

```bash
# single-end, gzipped input (-o defaults to reads.fqxv)
fqxv compress reads.fastq.gz

# paired-end at higher effort
fqxv compress sample_R1.fq.gz sample_R2.fq.gz -o sample.fqxv --level 7 --threads 16

# 10x single-cell (R1 + R2 + I1 + I2)
fqxv compress R1.fq R2.fq I1.fq I2.fq -o sample.fqxv

# lossy quality (Illumina 8-level binning)
fqxv compress reads.fastq -o reads.fqxv --quality-bin bin8

# maximum compression (deepest context + read reordering where it helps)
fqxv compress reads.fastq.gz --max

# estimate the ratio and archive size from a sample, writing nothing
fqxv compress reads.fastq.gz --estimate
```

## Estimating compression

`--estimate` tells you how well an input will compress **before** committing to a
full run. It codes a bounded sample of the leading reads (about one block, up to
~1M reads) with the real codecs at the chosen `--level`/`--quality-bin`, then
projects the whole-file archive and prints the result — no archive is written:

```text
reads.fastq.gz (436.93 MB)  →  estimated fqxv ~216.87 MB  (50% smaller, ~2.01x)

Estimated from a 1,048,576-read sample (1.31 GB uncompressed FASTQ):
  stream      compressed     share   rate
  names         6.42 KB     0.0%   0.05 bits/read
  sequence     17.11 MB    43.8%   1.355 bits/base
  quality      21.97 MB    56.2%   1.740 bits/base
  vs uncompressed FASTQ (~1.31 GB):  ~6.19x
```

The top line leads with the outcome: input size(s) → estimated archive size, as a
percent reduction and ratio against the input **on disk**. The per-stream table
below shows where the bytes go (names / sequence / quality), and the final line
gives the ratio against the **uncompressed** FASTQ. Multiple inputs are listed
individually under a combined headline; a streaming stdin input has no on-disk
size, so it reports the sample's own reduction and omits the whole-file
projection (pass a file to get one).

**Accuracy.** Because blocks are coded independently and the models are per-block
stationary, the sample's ratio is a faithful proxy for the whole file — in
practice within ~1–2% of a real run (every projected figure is marked `~`). The
projection scales the sample by the fraction of on-disk bytes it consumed, so it
works for plain and gzipped inputs alike.

**Reordering is a lower bound.** `--order any` / `--max` reordering is *not*
modeled — its cross-read redundancy grows with read count, so a small sample
can't capture it. With those flags the estimate is a **conservative lower bound**
(the real archive comes out that size or smaller), and the report says so.

## Notes

- **Losslessness:** read name + description, sequence, and quality are preserved
  exactly. The `+` separator line is normalized to a bare `+`.
- **Grouping:** interleaving mates shrinks the archive — near-identical mate
  names collapse to matches, and a spot's related reads sit together for the
  sequence model.
