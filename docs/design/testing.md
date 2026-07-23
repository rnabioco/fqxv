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
straight to the *decompressor*. We now have an equivalent — coverage-guided
`cargo-fuzz` targets over every decode entry point (see below). A local probe (30k
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
- ✅ **`cargo-fuzz` harness over every decoder.** `fuzz/fuzz_targets/` has seven
  coverage-guided targets — `container`, `rans`, `tokenizer`, `seq`, `fqzcomp`,
  `reorder`, `lroverlap` — each feeding arbitrary bytes straight to a decode entry
  point, the invariant being that a decoder returns `Ok`/`Err` on *any* input and
  never panics/aborts (release is `panic = "abort"`). This complements the in-tree
  proptest corruption harness (`crates/fqxv/tests/corruption.rs`), which mutates
  valid encodings in normal CI. See `fuzz/README.md`.
- ▶ **Extend the inline `decode_never_aborts_on_garbage` test to every codec.**
  Still open. The cheap inline random-garbage regression test lives in
  `fqxv-fqzcomp` and `fqxv-seq` (and, for the packed sub-path only, `fqxv-reorder`'s
  `refpack`); **`fqxv-rans`, `fqxv-tokenizer`, `fqxv-range`, and the top-level
  `fqxv-reorder` decode path still lack it.** Mirror it in those so a random-garbage
  smoke check runs without the nightly fuzz toolchain. (`fqxv-rans` has
  `decode_survives_mutation`, which mutates valid encodings rather than feeding
  arbitrary bytes — not equivalent.)
- ○ **tokenizer/reorder speculative allocations** cap at `1<<20`/`1<<22`, so the
  worst case is ~96 MB rather than an abort — lower risk, but converting them to
  `try_reserve` would remove the last alloc-abort paths for consistency.

## 2. Degenerate-distribution & raw-path cases (entropy coders)

Round-trip tests over random data almost never generate these:

- ✅ **Empty stream through compress *and* decompress** — htscodecs #99 was a
  divide-by-zero (SIGFPE) on empty input. `roundtrip_empty` in both `fqxv-rans`
  and `fqxv-range` now encodes *and* decodes empty input (not just decode-empty).
- ✅ **Single distinct symbol / all-same-byte with 12-bit frequency** — htscodecs
  1.5.2 fixed a SIMD mis-decode here. `avx2_encode_matches_scalar_skewed`
  (`fqxv-rans`) exercises the "SIMD ≡ scalar" invariant on degenerate cases:
  all-same (`vec![b'Q'; 40_000]`, freq == TOTFREQ), a dominant symbol (freq >
  2048), and a 2-symbol skew. Residual (▶): the AVX2-vs-scalar *proptests* are
  still seeded with a small random alphabet, not degenerate distributions — a
  proptest seeded with all-same / 2-symbol-at-4095/4096 would close the gap.
- ▶ **Raw/CAT passthrough path**, if/when added — htscodecs #144 fell through
  into the order-0 encoder *because the test suite omitted the raw path*. (No raw
  path exists today; this fires only if one is added.)
- ✅ **Malformed frequency table** (doesn't sum to total) → decode must reject.
  `decode_rejects_malformed_freq_table_without_aborting` and
  `from_freqs_rejects_malformed_tables` (`fqxv-rans`) assert `Err(Malformed)` on
  sum ≠ TOTFREQ and over-cap tables.
- ▶ **1–2 byte inputs** ✅ (`roundtrip_single_byte`, `roundtrip_short_and_odd_lengths`,
  range `roundtrip_single`) but **block-boundary-straddling large inputs** remain
  ○: tests cover the 32-state interleave boundary, not the htscodecs #128
  >1.25 MiB output-sizing boundary specifically.
- ○ **`output.len() <= bound(input)` property** for every encoder — guards the
  class of max-growth underestimates behind several htscodecs security fixes.
  Not yet written (the nearest test asserts a compression *ratio* on favorable
  data, not a worst-case growth bound).

## 3. Name tokenizer

- ▶ **Names with bytes ≥ 0x80** — htscodecs #105 was a signed-char crash on
  `[0x80, 0x0a]`. Rust avoids the UB, but index/sign logic can still misbehave.
  Still open: the arbitrary-name generators sample only from `b"AB:._-0129 "`
  (all < 0x80); the `0xFF` bytes in the corruption tests mutate *encoded* output,
  not name input.
- ▶ **Exactly `MAX_COL` and `MAX_COL`+1 tokens** — htscodecs 1.5.2 overflowed
  writing a 129th token. `MAX_COL` here is **63** (not 128); no test constructs a
  name that tokenizes to exactly the column-cap boundary. Still open.
- ○ Empty name ✅, single-char ✅, all-numeric ✅, varying column counts row-to-row
  ✅ (all in `roundtrip_varying_structure`); **all-hex** still missing (ties into
  the hex-aware token class, §5).
- ○ **UUID4 read-name corpus** (Nanopore) as a benchmark + roadmap: htscodecs
  #130/#131 show a transpose strategy beats the tokenizer for these; worth a
  hex-aware / leading-zero-preserving token class (see §5).

## 4. Sequence content & read reordering (SPRING / PgRC)

- ▶ **Exact A/C/G/T/N count *and* N-position preservation across reorder** —
  PgRC #6 drifted per-base counts because N-reads go to a separate subset. Partly
  covered: `reorder_free_preserves_records_as_a_set` /
  `reorder_rescue_preserves_records_as_a_set` assert the full record *multiset*
  survives (implying N positions), and `content_stats_and_metadata` asserts exact
  base counts — but only on the *plain* path, and the reorder fixture has no N
  bases. Still open: an explicit base-count + N-position assertion on an
  N-containing input *through reorder*.
- ✅ **Compress-twice, byte-identical, independent of `--threads`** — SPRING #35
  produced non-deterministic archives. `archive_is_deterministic_across_threads`
  (1 vs 4 threads, byte-identical), plus `reorder_paired_is_thread_count_deterministic`
  and `merge_roundtrips_and_is_deterministic`.
- ✅ **`--threads 0`** — SPRING #44 infinite-looped. We clamp (commit `93eb5ff`);
  `resolve_threads_zero_is_all_cores_and_explicit_is_clamped` (`fqxv-cli`) is the
  regression test.
- ✅ **Mixed/variable read lengths across a block boundary** — PgRC #2 is
  fixed-length-only (≤255 bp). `ragged_lengths_roundtrip_multiblock` (30 reads,
  10–310 bp, `block_reads: 5`) and `paired_split_spans_multiple_blocks` own it.
- ○ **RC palindromes and a read + its own reverse complement** — read + own
  revcomp clustering is ✅ (`revcomp_duplicate_flips_to_match`); a true palindrome
  (sequence == its own revcomp) is not exercised through clustering.
- ○ **N↔quality coupling and IUPAC codes** — IUPAC losslessness is covered at the
  codec level (round-trip proptests seeded with `RYKM` etc. in `fqxv-seq`,
  `refpack`, `fqxv-lroverlap`), but there is no container-level IUPAC test, no
  N↔quality independent-preservation test, and no test documenting the
  collapse-to-N stance.
- ○ Structural fixtures: empty read ✅ / all-N read ✅ (`handles_short_and_n_reads`);
  single read/single spot ✅; line count not a multiple of 4 ✅
  (`a_truncated_fastq_errors_in_both_order_modes`); seq/qual length mismatch →
  clean error ✅ (`compensating_seq_qual_mismatch_is_rejected`). Still missing:
  **511/512-bp boundary** (SPRING switches algorithms >511 bp).

## 5. Feature-request candidates (from upstream trackers)

- ✅ **Selectable lossy quality tiers.** SPRING offers lossless / Illumina-8-bin /
  QVZ / binary-threshold; we have `QualityBinning` (`Lossless`/`Bin8`/`Bin4`/`Bin2`
  plus `BinOnt`/`BinHifi`). The **binary 2-level threshold** ask is shipped as
  `Bin2` (`--quality-bin bin2`; fixed Q25 cutpoint — a user-selectable threshold
  value would extend it).
- **Hex-aware / UUID name tokenization.** htscodecs #131 (columnar by token
  type, leading-zero preservation) and #130 (transpose for UUID4). Adaptive:
  detect name style and switch strategy; regresses on Illumina if non-adaptive.
- **fqzcomp READ1/READ2 + strand selectors.** CRAM 3.1 uses multiple parameter
  sets at higher levels and per-record strand/mate hints; ties into our
  interleaved-spots invariant and `--level` mapping.
- **Archive introspection / partial decode.** *Partly shipped.* Introspection is
  done (`fqxv info`, and the `FQXF` footer index gives per-stream
  `(offset, length, crc32c)`); random-access column projection over HTTP `Range`
  exists in the Python `fqxv.remote` API. Remaining: no CLI subcommand for
  subset/row-range decode — partial decode is library/Python-only.
- ✅ **Streaming (`-`) stdin/stdout.** Shipped in #242: decompress reads an archive
  from stdin (`-`) and writes FASTQ to stdout (`-Z`/`-o -`); compress reads FASTQ
  from stdin. Carve-outs enforced: `--split`/`--recover` from stdin are rejected
  (no rewind), compress-from-stdin requires `-o`, and the archive itself is not
  written to stdout (needs a seekable footer). Test: `stdin_stream.rs`.

## 6. Interop conformance

- ○ Decode a handful of **htscodecs `tests/dat` / `tests/names` fixtures** and
  assert byte-identical plaintext — cross-implementation decode conformance
  against the reference, separate from internal round-trips. `CRAMcodecs.pdf` is
  the normative reference. The fqzcomp `.test` fixture format (one line of
  qualities per record, optional trailing READ1/READ2 flag) is a convenient
  interop format to adopt.

## 7. CLI

- ✅ `fqxv-cli` now has integration tests (`crates/fqxv-cli/tests/`):
  `batch.rs` (`info`/`verify`/`estimate` batch shapes, `--verify`, bad-input and
  `--interleaved 0` error paths), `quality_bin.rs` (compress→decompress with
  `--quality-bin`), and `stdin_stream.rs` (stdin round-trip + `--split`-from-stdin
  rejection). Plus the `resolve_threads` unit test for `--threads 0`.
- ▶ Remaining CLI gaps: no test exercises the **`--level` → order/block mapping**,
  and no end-to-end CLI test drives `--threads 0` (only the `resolve_threads`
  unit test does).
