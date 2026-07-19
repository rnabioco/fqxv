# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The on-disk `.fqxv` format is **not yet stable**: `FORMAT_VERSION` is bumped
freely and archives are not guaranteed to be readable across releases until a
`1.0.0` (each build reads only its own format version).

## [Unreleased]

### Added

- **Sequence-conditioned quality for long reads** — the quality coder can now
  condition each score on the read's bases (current base, next base, and the
  homopolymer run-length) instead of the read position, which carries no signal
  on long reads. The mode is chosen automatically by mean read length and recorded
  in a self-describing header byte, so short-read archives are byte-identical to
  before and keep decoding sequence and quality in parallel; long-read archives
  decode the sequence first and feed it to the quality decoder. On PacBio HiFi
  this cuts the quality stream ~10% (lossless), the dominant share of a HiFi
  archive. New `fqzcomp` API: `encode_seq`/`decode_seq`/`needs_sequence`; new
  random-access projection helpers `decode_quality_with_seq`/
  `quality_needs_sequence`. Bumps the `fqzcomp` stream version to 2.

- **Remote / parallel column projection** — the footer row-group index now
  records, per group, a `(offset, len, crc32c)` triple for each of the three
  coded streams (names, sequence, quality). A client can fetch the archive tail,
  parse the index, and issue a single range request for just one stream — read
  names are <1% of the archive, so an ID-only client fetches ~100× less than
  before, and a sequence-only client (k-mer screening, classification) skips the
  quality stream. The joint block content digest can't verify a single fetched
  stream, so each stream carries its own CRC-32C.
- **Random-access API** — a public, IO-free `Index` (parse from a seekable
  reader or a fetched suffix buffer via `Index::from_suffix`), `Index::byte_ranges`
  to turn `(groups, stream)` into byte ranges to `GET`, `Index::verify_stream`,
  and per-stream (`decode_names`/`decode_sequence`/`decode_quality`) and
  whole-block (`decode_block_contents`) decoders. The caller drives fetching
  (local `File`, `object_store`, async HTTP, …).
- **`compress --block-reads N`** — set the reads-per-row-group directly,
  decoupling random-access granularity from the `--level` effort knob. Smaller
  groups give finer remote access and more parallelism at some ratio cost.
- **Per-stream content digests** — the plain block payload now carries three
  xxh3-64 digests (names, sequence, quality) in place of the single joint digest.
  A post-decode mismatch (a codec round-tripping CRC-valid bytes into
  wrong-but-in-bounds output) now names the offending stream instead of only the
  block. Each digest still folds in `n_reads` and its stream's per-read lengths,
  so boundary pinning is unchanged; cost is 16 extra bytes per block.

### Changed

- **`FORMAT_VERSION` → 5** across the Unreleased changes (extended footer index,
  long-read overlap sequence codec, and the per-stream block digests above). As
  always in the pre-1.0 format, a build reads only its own version.
- **`inspect`** now sums per-stream sizes straight from the footer index instead
  of seeking to each block header — one footer read is the whole metadata cost.

## [0.2.0] - 2026-07-15

### Added

- **`compress --estimate`** — predict the compression ratio and archive size
  for a FASTQ without writing an archive. `--estimate tsv` emits machine-readable
  output, and `info`/`verify` now accept multiple files or directories for batch
  reporting.
- **Reference sequence coder for `--max`/reorder** — a SPRING-style codec that
  2-bit-packs the assembled global reference and entropy-codes it with a
  clean-room LZMA, adopted only when it beats the raw representation.
- **`--order shuffle` renumber mode** — a true SPRING-style read renumbering that
  discards the input permutation, reaching SPRING-competitive ratios on datasets
  where read order carries no information.
- **CLI run feedback** — a TTY-aware progress spinner and human-readable
  compress/decompress run summaries (suppressed under `--quiet` and when stderr
  is not a terminal).

### Changed

- **Quality coding** — fqzcomp now models quality over the symbol alphabet that
  actually occurs in the file rather than a fixed `0..QMAX` range, shrinking the
  quality stream on data with sparse quality alphabets. Byte-identical
  round-trips are preserved.

### Performance

- **Reorder merge** — the overlap-merge k-mer index now uses a rolling hash with
  sharding (~10% faster `--max`, byte-identical output), and the
  `merge_reference` vote-scatter is parallelized.

### Documentation

- VHS-rendered terminal demos on the docs site, refreshed benchmarks (including
  an `fqxv-shuffle` row), and clearer positioning.

## [0.1.0] - 2026-07-14

Initial release of `fqxv`, a Rust toolkit for lossless (opt-in lossy) archiving
of FASTQ. Codecs are clean-room implementations from specs and papers
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
- **`fqxv` CLI** — clap front-end (`fqxv` binary) with `compress`, `decompress`,
  `info`, and `verify` subcommands. Reads gzip-compressed or plain FASTQ, supports
  stdin/stdout streaming, and auto-detects interleaved paired streams. `--level
  1-9` maps to sequence context order and block size; `--max` is the best-ratio
  preset (deepest context plus read reordering, applied where it helps);
  `--order preserve|any|shuffle` names the order guarantee (`--keep-order`,
  `--no-rescue`, `--interleaved` refine it); `info --stats` reports content
  statistics and `verify --quick` does a fast footer-index CRC check. `--threads`
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
  clustering (reverse-complement aware) to exploit cross-read redundancy.
  Includes a whole-file **global-reference** sequence codec: an assembler builds
  a shared reference over the clustered reads and codes each read as a position
  on it, an **overlap-merge** pass consolidates the reference, and it is adopted
  only when it beats the block-local codecs (never worse). Assembly (windowed),
  the overlap-merge, and the reference coder are all parallel.
- **`fqxv-bytes`** — shared byte-serialization primitives (LEB128 varints,
  zig-zag) used across the codec crates.
- **Spot grouping / interleaving** — N-way grouping for paired mates and 10x
  R1/R2/I1/I2; blocks always hold whole spots and start on member 0 so they
  split cleanly. Grouped reorder records the group size and stores a permutation,
  so it is always order-preserving.
- **Parallelism** — block-level `rayon` parallelism throughout compress and
  decompress (assembly, overlap-merge, reference coding, per-block codecs, FASTQ
  parsing), all with **fixed, thread-independent boundaries**, so output is
  byte-identical regardless of thread count. No external/C compressor dependency.
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

[0.2.0]: https://github.com/rnabioco/fqxv/releases/tag/v0.2.0
[0.1.0]: https://github.com/rnabioco/fqxv/releases/tag/v0.1.0
