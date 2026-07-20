# fqxv compress

Compress one or more FASTQ files into a single `.fqxv` archive.

![fqxv compress demo](../images/compress.gif)

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
| `-f, --force` | Overwrite the output archive if it already exists. By default `compress` refuses to clobber an existing file and errors before doing any work. |
| `--verify` | After writing the archive, re-read and fully decode it to confirm it round-trips before reporting success — recommended before deleting the source FASTQ. See [Verifying on write](#verifying-on-write). |
| `--max` | Maximum-compression preset: deepest sequence context plus read reordering *where it helps* (applied to short reads, auto-skipped for long reads). Overrides `--level`/`--order`. |
| `--estimate` | Predict the archive size and compression ratio from a sample of the input, then exit **without writing anything**. See [Estimating compression](#estimating-compression). |
| `--threads <N>` | Worker threads (0 = all cores). |

### Advanced options

| Option | Description |
| --- | --- |
| `-l, --level <N>` | Effort 1–9; higher raises the sequence context order (up to order 11, reached at level 5) and the block size, and enables a hashed high-order tier at level 8+. Default: 5. |
| `--block-reads <N>` | Reads per row group, overriding the size `--level` would pick. Decouples random-access granularity from effort: smaller groups give finer remote/parallel access and more parallelism at some ratio cost; larger groups the reverse. Sequence order still follows `--level`. Ignored by the reorder path (`--order any`/`--max`). See [Row-group sizing](#row-group-sizing). |
| `--order <MODE>` | Read-order guarantee: `preserve` (default, restores original order), `any` (allows reordering for a better ratio; single-end order may change), or `shuffle` (like `any`, but discards order and regenerates purely positional names — reorder-lossy, single-end only). |
| `--interleaved <N>` | Interleaving of a *single* input, in members per spot (1 = single-end, 2 = paired as from `sracha get -Z`). Auto-detected from read names by default. Ignored with multiple inputs. |
| `--keep-order` | With `--order any`, force original read order to be restored (store a permutation, code names/quality in original order). Chosen automatically when it makes the archive smaller. |
| `--no-rescue` | With `--order any`, disable the adaptive assembly codecs (block-local literal-rescue and the whole-file global reference) and use the faster single-contig sequence codec only. |
| `--quality-bin <MODE>` | `lossless` (default), or a lossy table: `bin8`, `bin4`, `bin2` (Illumina), `ont`, `hifi` (long read). See [Lossy quality binning](#lossy-quality-binning). |
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

# lossy quality on long reads (match the table to the platform)
fqxv compress ont_reads.fastq -o ont_reads.fqxv --quality-bin ont
fqxv compress hifi_reads.fastq -o hifi_reads.fqxv --quality-bin hifi

# finer row groups for remote/parallel access (smaller range fetches)
fqxv compress reads.fastq.gz -o reads.fqxv --block-reads 65536

# maximum compression (deepest context + read reordering where it helps)
fqxv compress reads.fastq.gz --max

# verify the archive round-trips before trusting (or deleting) the source
fqxv compress reads.fastq.gz --verify

# estimate the ratio and archive size from a sample, writing nothing
fqxv compress reads.fastq.gz --estimate
```

## Verifying on write

`--verify` closes the window between writing an archive and trusting it. After the
archive is written, it is reopened and **fully decoded** — exercising every block
CRC-32C and content digest — and the decoded read count is checked against what was
written, all before `compress` reports success:

```bash
fqxv compress reads.fastq.gz --verify
```

This catches a codec or in-flight memory error that produced a CRC-valid but wrong
archive, which the archive's own on-disk checksums cannot detect on their own. It
is the check to run before deleting the source FASTQ.

- On any failure — a decode error or a read-count mismatch — the archive is left
  in place for inspection and the command exits non-zero. A bad archive is never
  reported as done.
- It adds a full decode pass, so expect roughly double the wall time (the
  `--max` / reorder path is the heaviest).
- It validates the bytes the encoder emitted, not long-term on-disk durability;
  for the latter, re-run [`fqxv verify`](verify.md) against the stored file later.

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

## Row-group sizing

The archive is a run of independently-coded **row groups** (blocks); the row
group is the unit of parallelism, of coarse random access, and of remote column
projection. By default its size comes from `--level` — higher effort uses larger
groups, which train the sequence model on more reads for a better ratio. Pass
`--block-reads <N>` to set it directly and decouple granularity from effort.

The trade-off is real. Smaller groups mean finer random access and more
parallelism, but a worse ratio (the order-*k* sequence model has fewer reads to
train on per group) and a slightly larger footer index. Larger groups compress
better but make the smallest independently-fetchable unit coarser. When archiving
to object storage where clients issue small `Range` reads — say, fetching just
the read names, or just one row group — a smaller `--block-reads` makes those
fetches cheaper; for a write-once/read-sequentially archive the `--level` default
is the right call.

The per-group and per-stream byte offsets recorded in the footer are what make
this projection possible without reading whole blocks — see
[Container format → Column projection](../design/container.md#footer-index-info-and-random-access).
The reorder path (`--order any`/`--max`) clusters globally and does not use this
sizing, so `--block-reads` is ignored there.

## Lossy quality binning

Quality is lossless by default. `--quality-bin` maps each quality byte through a
fixed table before coding — an explicit, opt-in lossy transform that shrinks the
quality stream (usually the largest part of the archive). Sequence and read names
are never touched, and a binned read still aligns exactly where the original did.

| Mode | Levels | Calibrated for |
| --- | --- | --- |
| `lossless` | all | default; nothing is discarded |
| `bin8` | 8 | Illumina standard binning |
| `bin4` | 4 | Illumina documented 4-level (NovaSeq X / RTA4) |
| `bin2` | 2 | custom, most aggressive |
| `ont` | 4 | Oxford Nanopore (CoLoRd ONT cutpoints) |
| `hifi` | 5 | PacBio HiFi (CoLoRd HiFi cutpoints; Q93 kept exact) |

**Match the table to the platform.** The tables are not interchangeable: the
Illumina bins are absolute-Phred cutpoints that collapse HiFi's narrow high-Q
band into a single level, and `ont` applied to HiFi data folds the Q93
max-quality symbol into the 26+ bin, destroying its application meaning (measured
mean |Δ| 42.84, 99.4% of bases changed). `hifi` keeps Q93 as its own level. On
ONT data the `ont` and `hifi` tables are byte-identical, since ONT never reaches
Q93.

On the `ecoli_ont` benchmark, `--quality-bin ont` cuts the quality stream from
163.7 MB to 47.0 MB (3.5×) — the whole archive from 2.79× to 6.06× — at mean
|Δ| 3.35. Cutpoints should ultimately be
judged by downstream fidelity, not raw ratio — see
[Long-read support](../design/longread.md).

## Notes

- **Losslessness:** read name + description, sequence, and quality are preserved
  exactly. The `+` separator line is normalized to a bare `+`.
- **Grouping:** interleaving mates shrinks the archive — near-identical mate
  names collapse to matches, and a spot's related reads sit together for the
  sequence model.
