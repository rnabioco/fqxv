# 🗜️ fqxv

A fast, reference-free **FASTQ archiver** for short-read sequencing data,
written in Rust. `fqxv` compresses each part of a FASTQ record with a codec
tuned to it — a context model for quality, an order-k model for sequence, a
positional tokenizer for names — and composes them into one parallel,
block-based container.

## Why fqxv

- **Reference-free & lossless** — no genome required; read name + description,
  sequence, and quality are preserved exactly (the redundant `+` line is
  normalized to a bare `+`, as SPRING and fqz_comp do).
- **Strong ratios** — clean-room implementations that match or beat the C
  reference tools stream-for-stream (see [Benchmarks](benchmarks.md)).
- **Parallel** — blocks compress and decompress across cores with `rayon`;
  output is deterministic regardless of thread count.
- **Paired & single-cell aware** — interleave R1/R2 (and 10x I1/I2) into one
  archive; split them back out, or stream interleaved to an aligner.
- **One crate per algorithm** — every codec is an independently usable,
  independently published Rust crate.

## Quick look

```bash
# single-end (gzip input auto-detected; -o defaults to reads.fqxv)
fqxv compress reads.fastq.gz
fqxv decompress reads.fqxv -o reads.fastq

# paired-end / single-cell: one archive, split back or stream to an aligner
fqxv compress sample_R1.fq.gz sample_R2.fq.gz -o sample.fqxv
fqxv decompress sample.fqxv --split out            # out_R1.fastq.gz, out_R2.fastq.gz
fqxv decompress sample.fqxv -Z | bwa mem -p ref.fa -  # interleaved, raw, to stdout

fqxv info sample.fqxv                              # layout, reads, per-stream sizes (--tsv/--json)
```

## How it works

A FASTQ record splits into three streams that compress very differently, and
`fqxv` gives each its own codec:

| Stream | Codec crate | Approach |
| --- | --- | --- |
| Quality scores | `fqxv-fqzcomp` | context model (prev-quals + position), range-coded |
| Sequence | `fqxv-seq` | order-k adaptive base model over a 2-bit alphabet |
| Read names | `fqxv-tokenizer` | positional tokens (match / delta / literal), rANS-coded |

The entropy backends — `fqxv-rans` (rANS Nx16, with an AVX2 path) and
`fqxv-range` (a Subbotin range coder) — are themselves standalone crates.

## Where to go next

- [Installation](getting-started/installation.md) — build the CLI and crates
- [Quick Start](getting-started/quickstart.md) — compress, inspect, decompress
- [CLI Reference](cli/index.md) — every command and flag
- [Design](design/index.md) — codecs, container format, and reordering
- [Benchmarks](benchmarks.md) — how fqxv stacks up against gzip, fqz_comp, SPRING
- [Acknowledgments](acknowledgments.md) — the prior work fqxv builds on
