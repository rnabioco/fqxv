# Quick Start

## Compress

```bash
fqxv compress reads.fastq -o reads.fqxv
```

Gzipped input is detected automatically, and `-o` is optional — it defaults to the
input's name with the FASTQ/gzip extension replaced by `.fqxv`:

```bash
fqxv compress reads.fastq.gz            # writes reads.fqxv
```

Tune effort with `--level` (1–9; higher raises the sequence context order) and
threads with `--threads` (default 16, capped at available cores; 0 = all cores):

```bash
fqxv compress reads.fastq.gz -o reads.fqxv --level 7 --threads 16
```

Not sure how well your data will compress? Add `--estimate` to predict the
archive size and ratio from a sample of the input, writing nothing:

```bash
fqxv compress reads.fastq.gz --estimate
# reads.fastq.gz (436.93 MB)  →  estimated fqxv ~216.87 MB  (50% smaller, ~2.01x)
```

See [`compress --estimate`](../cli/compress.md#estimating-compression) for
details and accuracy.

## Inspect

```bash
fqxv info reads.fqxv
```

```text
reads.fqxv
  layout         single-end (group size 1)
  reads          1000000
  blocks         4
  sequence order 11
  quality        lossless
  plus line      normalized
  names   6189536 bytes (13.4%)
  seq    17980884 bytes (38.9%)
  qual   22083418 bytes (47.7%)
```

## Decompress

Pick a destination — a file (`-o`), split mate files (`--split`), or a stdout
stream (`-Z`). A bare `decompress` with none of these errors rather than flooding
the terminal.

```bash
fqxv decompress reads.fqxv -o reads.fastq        # plain FASTQ
fqxv decompress reads.fqxv -o reads.fastq.gz     # block-gzip (BGZF)
```

## Paired-end and single-cell

Give multiple inputs to interleave per-spot files into one archive:

```bash
fqxv compress sample_R1.fq.gz sample_R2.fq.gz -o sample.fqxv   # paired
fqxv compress R1.fq R2.fq I1.fq I2.fq -o sample.fqxv           # 10x single-cell
```

Restore the separate files, or stream interleaved straight to an aligner:

```bash
fqxv decompress sample.fqxv --split out                  # out_R1.fastq.gz, out_R2.fastq.gz, ...
fqxv decompress sample.fqxv -Z | bwa mem -p ref.fa -     # interleaved, raw, on stdout
```

`--split` writes block-gzip `.fastq.gz` with `_R1`/`_R2` labels by default; add
`--no-gzip` for plain FASTQ or `--mate-style num` for `_1`/`_2` labels.

## Lossy quality (optional)

Quality is lossless by default. Opt into binning for smaller archives when you
don't need exact quality:

```bash
fqxv compress reads.fastq -o reads.fqxv --quality-bin bin8   # or bin4 / bin2
```

The `bin8`/`bin4`/`bin2` tables are Illumina-calibrated. Long reads have their
own, and the tables are not interchangeable — pick the one matching your
platform:

```bash
fqxv compress ont_reads.fastq -o ont_reads.fqxv --quality-bin ont    # Nanopore
fqxv compress hifi_reads.fastq -o hifi_reads.fqxv --quality-bin hifi # PacBio HiFi
```

See [Lossy quality binning](../cli/compress.md#lossy-quality-binning) for what
each table costs.
