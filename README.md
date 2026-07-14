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
# single-end (gzip input auto-detected; -o defaults to reads.fqxv)
fqxv compress reads.fastq.gz
fqxv decompress reads.fqxv -o reads.fastq

# paired-end / single-cell: interleave per-spot files into one archive
fqxv compress sample_R1.fq.gz sample_R2.fq.gz -o sample.fqxv   # paired
fqxv compress R1.fq R2.fq I1.fq I2.fq -o sample.fqxv           # 10x single-cell

# restore the separate mate files, or stream interleaved to an aligner
fqxv decompress sample.fqxv --split sample                     # sample_R1.fastq.gz, sample_R2.fastq.gz
fqxv decompress sample.fqxv --split s --no-gzip --mate-style num  # s_1.fastq, s_2.fastq (plain)
fqxv decompress sample.fqxv -Z | bwa mem -p ref.fa -          # interleaved, raw, to stdout

fqxv info sample.fqxv                                          # layout, reads, per-stream sizes (--tsv/--json)
fqxv info sample.fqxv --stats                                  # + read-length, GC%, quality dist (decodes)
fqxv verify sample.fqxv                                        # CRC integrity check (exit non-zero if corrupt)
```

Combining mates/index reads shrinks the archive (near-identical mate names
collapse; the sequence model sees a spot's related reads together) and keeps one
file per sample. `compress`/`decompress` are `rayon`-parallel (`--threads`).
Lossless by default; `--quality-bin {bin8,bin4,bin2}` opts into lossy quality.

`decompress` needs an explicit destination — `-o FILE`, `--split PREFIX`, or
`-Z/--stdout` — so a bare invocation never floods the terminal. `--split` writes
block-gzip (BGZF) `*_R1.fastq.gz`/`*_R2.fastq.gz` by default; `--no-gzip` emits
plain FASTQ and `--mate-style num` uses `_1`/`_2` labels. For `-o FILE` the
compression follows the extension (`.gz` → BGZF); `-Z` streams raw FASTQ to stdout.
Every decode verifies each block's CRC and content digest, and a file `decompress`
also confirms the decoded read count against the archive footer, so a truncated
archive fails loudly instead of yielding a short, silent result.

## Milestones

- **M0** — benchmark harness + baselines (`bench/`)
- **M1** — `fqxv-rans` (the SIMD showcase)
- **M2** — `fqxv-range` + `fqxv-fqzcomp` (quality)
- **M3** — `fqxv-tokenizer` + `fqxv-seq`
- **M4** — `fqxv-reorder`
- **M5** — `fqxv` container + CLI

## Acknowledgments

`fqxv` stands on a large body of prior work. Everything here is a clean-room
implementation from public specifications and papers — no third-party source is
vendored — but these projects and their authors made it possible, and we
cross-checked against several of them for correctness:

- **htscodecs** ([samtools/htscodecs](https://github.com/samtools/htscodecs),
  James Bonfield / Genome Research Ltd) and the [CRAM 3.1 codecs
  spec](https://samtools.github.io/hts-specs/CRAMcodecs.pdf) — the reference for
  our rANS Nx16 coder, fqzcomp quality model, and name tokenizer.
- **fqzcomp** (James Bonfield) — the quality-score context model our
  `fqxv-fqzcomp` codec is modeled on.
- **noodles** ([zaeleus/noodles](https://github.com/zaeleus/noodles),
  Michael Macias) — Rust CRAM codec implementation we cross-checked test vectors
  against.
- **rANS / ryg_rans** — Jarek Duda's asymmetric numeral systems and Fabien
  Giesen's `ryg_rans` (public domain / CC0), plus Eugene Shelwien's range-coder
  design, underpin our entropy coders.
- **SPRING** (Chandak et al., *Bioinformatics* 2019) and **PgRC2** (Kowalski &
  Grabowski, *Bioinformatics* 2025) — the algorithmic references for the
  read-reordering engine.

See [`THIRD-PARTY-NOTICES.md`](THIRD-PARTY-NOTICES.md) for licenses and full
attribution.

## License

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or
[MIT license](LICENSE-MIT) at your option.
