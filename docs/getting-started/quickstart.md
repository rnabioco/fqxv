# Quick Start

## Compress

```bash
fqxv compress reads.fastq -o reads.fqxv
```

Gzipped input is detected automatically:

```bash
fqxv compress reads.fastq.gz -o reads.fqxv
```

Tune effort with `--level` (1–9; higher raises the sequence context order) and
threads with `--threads` (default 16, capped at available cores; 0 = all cores):

```bash
fqxv compress reads.fastq.gz -o reads.fqxv --level 7 --threads 16
```

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

```bash
fqxv decompress reads.fqxv -o reads.fastq
```

## Paired-end and single-cell

Give multiple inputs to interleave per-spot files into one archive:

```bash
fqxv compress R1.fq.gz R2.fq.gz -o sample.fqxv           # paired
fqxv compress R1.fq R2.fq I1.fq I2.fq -o sample.fqxv     # 10x single-cell
```

Restore the separate files, or stream interleaved straight to an aligner:

```bash
fqxv decompress sample.fqxv --split out                  # out_1.fastq, out_2.fastq, ...
fqxv decompress sample.fqxv | bwa mem -p ref.fa -         # interleaved on stdout
```

## Lossy quality (optional)

Quality is lossless by default. Opt into Illumina-style binning for smaller
archives when you don't need exact quality:

```bash
fqxv compress reads.fastq -o reads.fqxv --quality-bin bin8   # or bin4 / bin2
```
