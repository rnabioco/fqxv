# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The on-disk `.fqxv` format is **not yet stable**: `FORMAT_VERSION` is bumped
freely and archives are not guaranteed to be readable across releases until a
`1.0.0` (each build reads only its own format version).

## [Unreleased]

### Added

- **Sequence-conditioned, context-mixed quality for long reads** — the quality
  coder now conditions each score on the read's bases and recent qualities instead
  of read position (which carries no signal on long reads). It mixes several
  context models of increasing richness (coarse/mid/rich) with adaptive,
  confidence-gated weights — a logistic mixer — because a per-block adaptive model
  can't exploit a richer *single* context but *can* blend a well-trained coarse one
  with a sparse rich one. On PacBio HiFi this takes the quality stream (the dominant
  share of a HiFi archive) below CoLoRd, lossless. The mode is chosen automatically
  by mean read length and recorded in a self-describing header byte, so short-read
  archives are byte-identical to before; long-read archives decode the sequence
  first and feed it to the quality decoder. The mixer is fixed-point/integer
  throughout, so archives are bit-identical across platforms. New `fqzcomp` API:
  `encode_seq`/`decode_seq`/`needs_sequence`; new random-access projection helpers
  `decode_quality_with_seq`/`quality_needs_sequence`.

- **Long-read (ONT / PacBio) support** — platform-aware compression for
  Nanopore and PacBio reads: `--platform illumina|nanopore|pacbio` (with header
  auto-detection), long-read quality binning (`--quality-bin ont|hifi`, matching
  CoLoRd's cutpoints), and a dedicated long-read overlap sequence codec that
  drives its minimizer sketch from the detected platform.
- **Long-read overlap sequence codec (`fqxv-lroverlap`)** — a new cross-read
  overlap codec (minimizers → overlaps → layout → consensus → per-read banded
  edit script → rANS) wired into the container as the sequence path for long-read
  blocks, auto-selected and kept only when it beats order-k. It reaches CoLoRd
  parity on ONT/HiFi (e.g. ~0.653 → ~0.067 bits/base at depth) and codes its edit
  streams (substitutions, ops, insertions) with per-stream context models rather
  than a flat order-0 code — worth a further −4.3% on the ONT sequence stream,
  losslessly. A WFA aligner (HiFi) and an AVX2 anti-diagonal aligner accelerate
  encoding, both byte-identical to the scalar reference.
- **Closed-syncmer seeding for ONT** — the overlap codec seeded from window
  minimizers, which a base error in *any* k-mer of the window can deselect; at
  Nanopore's ~10% error that is the dominant loss of shared anchors. ONT now seeds
  with closed syncmers (a k-mer is an anchor when its minimal `s`-mer sits at its
  first or last position), so selection depends only on the k-mer's own bases and
  an intact shared k-mer is co-selected regardless of neighbouring errors. Density
  is `2/(w+1)` as before, so `k`, anchor count and specificity are unchanged — only
  conservation improves: ONT overlap-coded sequence goes 1.631 → 1.559 bits/base
  against an order-k baseline of 1.808, widening with per-block coverage. Seeding
  is encode-only, so archives stay decodable by any reader. PacBio keeps window
  minimizers, which are already near-optimal below ~1% error.
- **Content-based platform detection** — platform was detected from read-name
  grammar alone, so SRA-reformatted runs (bare `SRR…` headers) recorded `unknown`
  and were handed the Nanopore sketch. When names carry no platform signal, fqxv
  now classifies long-read runs by mean per-base quality, which separates the
  platforms cleanly (measured across 24 corpus accessions with ENA ground truth:
  Nanopore 6.9–23.5, PacBio 36.9–84.5). Ambiguous data stays `unknown`, so a wrong
  platform is still never recorded. On a real HiFi run this restores the
  low-divergence WFA path and cuts encode time ~31%.
- **Shared whole-file reference for long reads** — the overlap codec reaches
  CoLoRd parity *within* a block, but the container re-assembled and re-stored the
  same consensus reference in every 256 MiB block. It now assembles one consensus
  over the whole file and stores it **once** in a framed region between the header
  and the first block, coding every block's reads against that frozen frame — so
  the genome is stored once, not once per block. Auto-selected for long-read input,
  behind a whole-file never-worse gate against order-k, and gated on the
  `GLOBAL_REFERENCE` feature bit (which is set in the plain layout now as well as
  the reorder layout, so a reader without it refuses rather than misreads). On
  PacBio HiFi this drops the sequence stream ~0.102 → ~0.084 bits/base at two
  blocks, widening with block count. New `fqxv-lroverlap` API: `Reference`,
  `build_reference`, `encode_against`/`decode_against`.
- **Python bindings (`fqxv` via PyO3 / maturin)** — read-only access from
  Python: a streaming record iterator and column projection (fetch just names,
  sequence, or quality) over an existing archive.
- **`compress --verify`** — an opt-in read-after-write check that round-trips the
  freshly written archive back to the original records before exiting, backed by a
  public `verify_roundtrip` in the library.
- **Block sync markers + footer-independent recovery** — each block carries a
  sync marker so a reader can resynchronize and recover blocks even when the
  footer index is missing or truncated.
- **`compress -f`/`--force`** — compression now refuses to overwrite an existing
  output unless `--force` is given.
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
- **Crash-safe compress output** — compression writes to a sibling temp file and
  atomically renames it into place only once the whole archive (header, blocks,
  *and* footer trailer) is on disk. Interrupting a stream mid-run (Ctrl-C on
  `sracha get -Z | fqxv compress -`, say) or hitting any error no longer leaves a
  corrupt, footer-less `.fqxv` at the destination: the partial temp is removed on
  a `?` bail and by a SIGINT/SIGTERM/SIGHUP handler, and the destination path only
  ever holds a complete archive.
- **Live compress progress** — the compress indicator now reports how much data
  has been processed and the rate. With a known input size (a file) it renders a
  percentage bar; for a stdin stream of unknown length it shows a bytes + rate
  readout, so a pause waiting on an upstream producer reads as `0 B` rather than a
  hang.

### Changed

- **On-disk format versioning is now `FORMAT_MAJOR.FORMAT_MINOR` (currently
  `1.0`)**, replacing the single monotonic `FORMAT_VERSION` integer. A reader
  refuses a differing major and tolerates a newer minor (backward-compatible
  additions), a documented forward-compatibility contract for the pre-stable
  format. This release's format carries the long-read overlap codec, the extended
  per-stream footer index, and the per-stream block digests below.
- **`platform` has its own header byte** — the platform tag no longer shares bits
  with the flags byte (a prior collision made Illumina `--order any` archives
  undecodable).
- **`inspect`** now sums per-stream sizes straight from the footer index instead
  of seeking to each block header — one footer read is the whole metadata cost.
- **`fqxv-dna` primitives crate** — the 2-bit ACGT lookup and reverse-complement
  helpers are extracted into a shared leaf crate, and the `fqxv-reorder` monolith
  is split into focused modules.
- **Default log verbosity is `warn`** — routine per-run `info` diagnostics now
  require `-v`, so they no longer interleave with and smear the live compress
  indicator on a shared stderr (`-vv` = debug, `-vvv` = trace with targets).
- **`--verify` verifies before publishing** — the read-after-write check now runs
  against the temp file, so a failed verification leaves the unverified output off
  the destination path (kept aside for inspection) and exits non-zero, rather than
  leaving a suspect archive in place.

### Performance

- **Single-end compress** now streams through a block pipeline instead of
  buffering the whole input.
- **Reorder read storage** flattens `cl_reads` into an arena, cutting ~168 MB of
  peak memory with byte-identical output.
- **Quality models** are sized to the alphabet actually present in the file
  (~17% faster quality coding, ~45 MB less memory).

### Fixed

- **Silent quality corruption** on certain inputs and a `--block-reads`
  compression-bomb path, both surfaced by CLI stress testing.
- **fqzcomp** no longer rejects legitimately compressible quality (which could
  drop data).
- **Read-name headers** are preserved byte-exactly.
- **Long-read memory** — the overlap codec's overlap/layout/placement and the
  banded aligner's DP matrix are now bounded, fixing OOMs on high-error ONT and
  large amplicon inputs (with O(n) placement).
- Assorted CLI fixes: `--order any` no longer aborts on a truncated FASTQ that
  `preserve` rejects, `--estimate` handles empty input, `--interleaved 0` is
  rejected, and `--verify` honors `--threads`.
- **Python binding tests** now build their own tiny `.fqxv` fixture with the CLI
  instead of latching onto whatever archive sits at the repo root (the suite ran
  for minutes against a large sample and could fail on a half-written one); they
  now finish in well under a second. The `Info` repr also shows the container
  version as `format=major.minor` instead of the raw packed integer.

### Security

- **Decode-path hardening** — a fuzzing-driven pass (cargo-fuzz targets over the
  public decode entries, plus a weekly fuzz schedule) closed multiple
  decompression bombs across the rANS, quality, sequence, name, reorder, and
  container decoders: untrusted length/size headers can no longer trigger
  unbounded allocation.

### Documentation

- Long-read benchmark tables, a format-comparison page, a Python API reference,
  and `--quality-bin ont|hifi` documentation.

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
