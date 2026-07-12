# Testing & Robustness

> Developer note. A running map of test coverage, decoder-robustness
> guarantees, and roadmap items — derived from the issue trackers and fuzzing
> harnesses of the tools `fqxv` reimplements (**htscodecs** for rANS Nx16, the
> fqzcomp quality model, and the name tokenizer; **fqzcomp**, **SPRING**, and
> **PgRC/PgRC2** for reordering) plus a local robustness probe of our own
> decoders. Each item cites the upstream bug or behaviour it guards against.

Status legend: ✅ done · ▶ recommended next · ○ backlog.

## What already exists

Every codec crate has inline `#[test]` round-trips plus one `proptest!` block;
`fqxv/src/container.rs` has integration coverage (pairing, interleaving,
single-cell, reorder modes, `+`-normalization, truncation). All decoders return
`Result` and guard header/length reads with `Malformed`. The gaps below are
where round-trip-only testing structurally can't reach: malformed input,
degenerate distributions, and cross-tool edge cases.

## 1. Decoder robustness on untrusted input (highest priority)

htscodecs fuzzes every one of these codecs via OSS-Fuzz, feeding arbitrary bytes
straight to the *decompressor*. We have no equivalent yet. A local probe (30k
random inputs per codec) found no panics, but **targeted adversarial headers did
crash the process**:

- ✅ **fqzcomp/seq unbounded allocation from length headers.** A ~13-byte
  `fqzcomp` stream declaring a huge length *count* aborted the process
  (`memory allocation of 288230376151711740 bytes failed`, SIGABRT); a 12-byte
  stream declaring `1000 reads × u32::MAX` length spun a ~4.3 TB decode
  (>120 s timeout). Fixed by using `try_reserve` instead of
  `with_capacity`/`resize`, a checked sum on `total`, and a decompression-bomb
  guard that bounds declared output by `payload_len × max-symbols-per-byte` — a
  principled ceiling from the range model's `1<<13` frequency cap, with wide
  margin so legitimate streams never trip it. Regression tests:
  `rejects_huge_length_count`, `rejects_huge_total_length`,
  `decode_never_aborts_on_garbage` in both `fqxv-fqzcomp` and `fqxv-seq`.
  Mirrors the htscodecs fqzcomp 1-bp buffer-overrun fix (NEWS 1.5.2).
- ▶ **Extend `decode_never_aborts_on_garbage` to every codec** (rans,
  tokenizer, reorder, and the `fqxv` container). Cheap, high-value.
- ▶ **Adopt `cargo-fuzz`** with one target per decoder (`*_fuzz`) plus
  round-trip targets (`*_fuzzrt`), matching htscodecs' split. No fuzz harness
  exists today.
- ○ **tokenizer/reorder speculative allocations** cap at `1<<20`/`1<<22`, so the
  worst case is ~96 MB rather than an abort — lower risk, but converting them to
  `try_reserve` would remove the last alloc-abort paths for consistency.

## 2. Degenerate-distribution & raw-path cases (entropy coders)

Round-trip tests over random data almost never generate these:

- ▶ **Empty stream through compress *and* decompress** — htscodecs #99 was a
  divide-by-zero (SIGFPE) on empty input. We have decode-empty tests; add
  compress-empty for rans/range.
- ▶ **Single distinct symbol / all-same-byte with 12-bit frequency** — htscodecs
  1.5.2 fixed a SIMD mis-decode here. Directly exercises our "SIMD ≡ scalar"
  invariant; add an explicit AVX2-vs-scalar proptest seeded with degenerate
  distributions (all-same, 2-symbol at 4095/4096), not just uniform random.
- ▶ **Raw/CAT passthrough path**, if/when added — htscodecs #144 fell through
  into the order-0 encoder *because the test suite omitted the raw path*.
- ○ **Malformed frequency table** (doesn't sum to total) → decode must reject,
  not read uninitialised entries (htscodecs NEWS 1.2.2).
- ○ **1–2 byte inputs** and **block-boundary-straddling large inputs**
  (htscodecs #128: >1.25 MiB output-sizing bug).
- ○ **`output.len() <= bound(input)` property** for every encoder — guards the
  class of max-growth underestimates behind several htscodecs security fixes.

## 3. Name tokenizer

- ▶ **Names with bytes ≥ 0x80** — htscodecs #105 was a signed-char crash on
  `[0x80, 0x0a]`. Rust avoids the UB, but index/sign logic can still misbehave.
- ▶ **Exactly 128 and 129 tokens** — htscodecs 1.5.2 overflowed writing a 129th.
  We have `MAX_COL`; test a name that tokenizes to exactly the boundary.
- ○ Empty name, single-char, all-numeric, all-hex; varying column counts
  row-to-row (cold-start delta state).
- ○ **UUID4 read-name corpus** (Nanopore) as a benchmark + roadmap: htscodecs
  #130/#131 show a transpose strategy beats the tokenizer for these; worth a
  hex-aware / leading-zero-preserving token class (see §5).

## 4. Sequence content & read reordering (SPRING / PgRC)

- ▶ **Exact A/C/G/T/N count *and* N-position preservation across reorder** —
  PgRC #6 drifted per-base counts because N-reads go to a separate subset. Our
  strongest reorder guard: assert exact base composition and N positions survive.
- ▶ **Compress-twice, byte-identical, independent of `--threads`** — SPRING #35
  produced non-deterministic archives; this directly asserts our determinism
  invariant.
- ▶ **`--threads 0`** — SPRING #44 infinite-looped. We clamp (commit `93eb5ff`);
  add a regression test at the CLI/container boundary.
- ▶ **Mixed/variable read lengths across a block boundary** — PgRC #2 is
  fixed-length-only (≤255 bp); variable length is the recurring top user ask and
  something we support, so we should own the test.
- ○ **RC palindromes and a read + its own reverse complement** — core of the
  reorder engine; assert clustering doesn't corrupt orientation.
- ○ **N↔quality coupling and IUPAC codes** — fqzcomp silently rewrites quality
  for N bases and collapses non-ACGTN to N (both lossy). Test whether we
  preserve base/quality independently and IUPAC losslessly, and document the
  stance either way (relates to the `+`-normalization deviation).
- ○ Structural fixtures: empty read; single read/single spot; all-N read;
  511/512-bp boundary (SPRING switches algorithms >511 bp); line count not a
  multiple of 4; seq/qual length mismatch → clean error.

## 5. Feature-request candidates (from upstream trackers)

- **Selectable lossy quality tiers.** SPRING offers lossless / Illumina-8-bin /
  QVZ / binary-threshold; we have `QualityBinning` (Bin8/4/2). A **binary
  2-level threshold** mode is a common explicit ask.
- **Hex-aware / UUID name tokenization.** htscodecs #131 (columnar by token
  type, leading-zero preservation) and #130 (transpose for UUID4). Adaptive:
  detect name style and switch strategy; regresses on Illumina if non-adaptive.
- **fqzcomp READ1/READ2 + strand selectors.** CRAM 3.1 uses multiple parameter
  sets at higher levels and per-record strand/mate hints; ties into our
  interleaved-spots invariant and `--level` mapping.
- **Archive introspection / partial decode.** SPRING #38/#42 show demand for
  documented, inspectable internals and subset/random-access decode — extends
  our existing `inspect`.
- **Streaming (`-`) stdin/stdout.** Recurring cross-tool need; confirm/support.

## 6. Interop conformance

- ○ Decode a handful of **htscodecs `tests/dat` / `tests/names` fixtures** and
  assert byte-identical plaintext — cross-implementation decode conformance
  against the reference, separate from internal round-trips. `CRAMcodecs.pdf` is
  the normative reference. The fqzcomp `.test` fixture format (one line of
  qualities per record, optional trailing READ1/READ2 flag) is a convenient
  interop format to adopt.

## 7. CLI

- ○ `fqxv-cli` has **zero tests**. Add end-to-end round-trips over small FASTQ
  fixtures (compress→decompress, byte-exact at FASTQ-content level), the
  `--level` mapping, and error paths (bad input, `--threads 0`).
