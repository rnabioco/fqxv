# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The on-disk `.fqxv` format is **not yet stable**: `FORMAT_VERSION` is `0` and
archives are not guaranteed to be readable across releases until a `1.0.0`.

## [Unreleased]

### Added

- **Grouped (paired / single-cell) read reordering.** The global-cluster reorder
  path now accepts interleaved input: reads are clustered ignoring mate
  structure, the group size is recorded in the header, and a stored permutation
  reconstructs the original spot interleaving on `decompress` and `--split`.
  Grouped reorder is therefore always order-preserving.

### Changed

- **CLI: `--reorder`/`--keep-order` replaced by `--order preserve|any`** (an
  Advanced option; default `preserve`). The flag now names the guarantee the user
  cares about — whether original read order survives — rather than the reordering
  mechanism. `any` allows reordering for a better ratio (single-end order may
  change); grouped input still round-trips in order. The library `Params.reorder`
  / `Params.keep_order` fields are unchanged.

## [0.1.0] - 2026-07-12

Initial release of `fqxv`, a Rust toolkit for lossless (opt-in lossy) archiving
of short-read FASTQ. Codecs are clean-room implementations from specs and papers
(CRAM 3.1 codecs, fqzcomp/SPRING/PgRC2); see `THIRD-PARTY-NOTICES.md`.

### Added

- **`fqxv` container format and library** — the `.fqxv` on-disk layout: a header
  followed by independent, parallel-codable **row groups** (blocks), a footer
  index, and an EOF trailer. Each row group splits FASTQ into three streams
  (names, sequence, quality) handled by three codecs and is byte-budgeted. The
  footer index makes `inspect` O(row groups) rather than O(bytes) and enables
  coarse random access (seek to and decode only the row groups overlapping a read
  range). The terminator lets the same file serve both a streaming reader and a
  seeking reader. Public API: `compress`, `compress_auto`, `compress_multi`,
  `compress_interleaved`, `decompress`, `decompress_split`, `inspect`, `peek`,
  and the `Info`/`Params`/`Stats` types.
- **`fqxv` CLI** — clap front-end (`fqxv` binary) with `compress`,
  `decompress`, and `info` subcommands. Reads gzip-compressed or plain FASTQ,
  supports stdin/stdout streaming, and auto-detects interleaved paired streams.
  `--level 1-9` maps to sequence context order and block size; `--threads`
  defaults to 16 and is clamped to available cores.
- **`fqxv-rans`** — rANS Nx16 entropy coder (CRAM 3.1) with 32 interleaved
  states and order-0/order-1 models. Backends are selected at runtime via CPU
  feature detection: scalar (correctness reference, all orders) plus SIMD
  order-0 encode and decode on **AVX2** and **AVX-512**, dispatched to the
  widest available path. Every backend produces byte-identical output.
- **`fqxv-range`** — serial binary range coder with adaptive bit models, the
  arithmetic-coding primitive underlying the quality and sequence codecs.
- **`fqxv-fqzcomp`** — fqz_comp-class quality-score context model. Owns
  `QualityBinning`: lossless by default, with opt-in lossy `Bin8`/`Bin4`/`Bin2`
  modes (re-exported as `fqxv::QualityBinning`).
- **`fqxv-tokenizer`** — positional read-name tokenizer with per-column delta
  bucketing and per-role payload streams; rANS entropy backend.
- **`fqxv-seq`** — order-k adaptive context model over 2-bit ACGT symbols
  (range-coded, variable read lengths); non-ACGT bytes go to an exception list.
- **`fqxv-reorder`** — PgRC2/SPRING-class read reordering via canonical-minimizer
  clustering (reverse-complement aware) to exploit cross-read redundancy;
  clustered differential codec for single-end and paired container modes.
- **Spot grouping / interleaving** — N-way grouping for paired mates and 10x
  R1/R2/I1/I2; blocks always hold whole spots and start on member 0 so they
  split cleanly.
- **Parallelism** — block-level `rayon` parallelism on both compress and
  decompress paths, plus parallel FASTQ parsing and a pipelined reader. Output
  is byte-identical regardless of thread count.
- **Logging** — `tracing`-based logging with `-v`/`-vv`/`-vvv` verbosity.
- **Documentation site** — zensical-based GitHub Pages site.
- **Benchmark harness** — separate `bench/` harness (pixi env, Slurm) for
  comparing against fqz_comp, SPRING, PgRC2, zstd, and gzip, plus criterion
  microbenchmarks for the rANS hot paths.

### Security

- Guarded the fqzcomp and seq decoders against length-header allocation bombs
  (untrusted length prefixes can no longer trigger unbounded allocation).

### Notes

- **`+` normalization** — the optional repeated `+` header line is dropped (as
  SPRING and fqz_comp do). Name, sequence, and quality are otherwise preserved
  exactly; this is the one documented deviation from byte-losslessness.

[Unreleased]: https://github.com/rnabioco/fqxv/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/rnabioco/fqxv/releases/tag/v0.1.0
