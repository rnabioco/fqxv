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
```

## Notes

- **Losslessness:** read name + description, sequence, and quality are preserved
  exactly. The `+` separator line is normalized to a bare `+`.
- **Grouping:** interleaving mates shrinks the archive — near-identical mate
  names collapse to matches, and a spot's related reads sit together for the
  sequence model.
