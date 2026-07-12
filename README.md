# 🗜️ fqxv

A Rust toolkit for archiving short-read FASTQ, built as a workspace of
one-crate-per-algorithm codecs plus a container format and CLI.

> Status: **early development.** Nothing is stable yet. See
> [the plan](#milestones) below.

## Why

Illumina FASTQ splits into three streams that compress very differently, and the
wins are additive:

| Stream | Share of a lossless archive | What moves it |
| --- | --- | --- |
| Quality scores | ~50–70% | context-model entropy coder (fqzcomp-class) |
| Sequence (bases) | most of the rest | read reordering / assembly (PgRC2/SPRING-class) |
| Read names | small | positional tokenizer |

No mature FASTQ-domain compressor exists as a Rust crate today. `fqxv` fills that
niche with clean-room implementations from the [CRAM 3.1 codecs
spec](https://samtools.github.io/hts-specs/CRAMcodecs.pdf) and the source papers.

## Workspace

| Crate | Role |
| --- | --- |
| `fqxv-rans` | rANS Nx16 entropy coder — scalar + SSE4.2 + AVX2, runtime dispatch |
| `fqxv-range` | binary range coder + adaptive models (serial) |
| `fqxv-fqzcomp` | quality context model over `fqxv-range`; opt-in lossy binning |
| `fqxv-tokenizer` | positional read-name tokenizer; entropy backend = `fqxv-rans` |
| `fqxv-seq` | 2-bit / variable-length base packing + order-k model |
| `fqxv-reorder` | PgRC2/SPRING-class read reordering engine |
| `fqxv` | container format library; composes the above |
| `fqxv-cli` | the `fqxv` command-line binary |

Every crate is independently publishable to crates.io and dual-licensed
**MIT OR Apache-2.0**.

## Design principles

- **Clean-room.** Implemented from specs/papers, not translated from C. See
  [`THIRD-PARTY-NOTICES.md`](THIRD-PARTY-NOTICES.md).
- **Rayon everywhere.** Block/record parallelism throughout; output is
  byte-identical regardless of thread count.
- **Benchmark-driven.** The CLI ships only codecs that beat the field in the
  `bench/` harness (fqzcomp5, SPRING, PgRC2, zstd -19, gzip baselines).

## Usage

```bash
# single-end
fqxv compress reads.fastq.gz -o reads.fqxv          # gzip input auto-detected
fqxv decompress reads.fqxv -o reads.fastq

# paired-end / single-cell: interleave per-spot files into one archive
fqxv compress R1.fq.gz R2.fq.gz -o sample.fqxv                 # paired
fqxv compress R1.fq R2.fq I1.fq I2.fq -o sample.fqxv          # 10x single-cell

# restore the separate files, or stream interleaved to an aligner
fqxv decompress sample.fqxv --split out                       # out_1.fastq, out_2.fastq, ...
fqxv decompress sample.fqxv | bwa mem -p ref.fa -             # interleaved to stdout

fqxv info sample.fqxv                                          # layout, reads, per-stream sizes
```

Combining mates/index reads shrinks the archive (near-identical mate names
collapse; the sequence model sees a spot's related reads together) and keeps one
file per sample. `compress`/`decompress` are `rayon`-parallel (`--threads`).
Lossless by default; `--quality-bin {bin8,bin4,bin2}` opts into lossy quality.

## Milestones

- **M0** — benchmark harness + baselines (`bench/`)
- **M1** — `fqxv-rans` (the SIMD showcase)
- **M2** — `fqxv-range` + `fqxv-fqzcomp` (quality)
- **M3** — `fqxv-tokenizer` + `fqxv-seq`
- **M4** — `fqxv-reorder`
- **M5** — `fqxv` container + CLI

## License

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or
[MIT license](LICENSE-MIT) at your option.
