# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

The on-disk `.fqxv` format is **stable at 1.0**. That version is independent of the
release versions below — the format is at 1.0 while the crates are still 0.x, and
the two numbers move separately. Archives written by any format-1.x writer remain
readable by later releases: a reader accepts its own format major version and
tolerates newer minors, skips additions it can safely ignore (a non-critical header
extension record), and refuses the ones it cannot — a required-feature bit or a
per-block codec method it does not know — with an "upgrade fqxv" error rather than
misreading. Compatibility fails loudly, never silently. A format major bump would
be announced as a breaking change; see the
[evolution policy](docs/design/container.md#versioning-and-evolution-policy).

## [Unreleased]

### Fixed

- **Single-cell interleaving is no longer archived as paired.** Layout detection
  only ever compared records pairwise, so a 4-member spot whose members share one
  read name — 10x R1/R2/I1/I2 as `fasterq-dump` and `sracha` emit it, where the
  members differ only in `length=` — looked like a mate pair at every offset and
  was recorded as *paired*. `fqxv info` reported twice the real spot count, and
  `--split` wrote two mate files that each interleaved two different slots at two
  different read lengths, with no warning. Detection now cuts the peeked records
  into runs of consecutive same-spot names and infers the common run length (up to
  4), so 3- and 4-member layouts are detected and restored to one file per slot.
  Archives written before this are unaffected: the stream itself was always
  byte-lossless, only the recorded grouping was wrong.
- **Interleaved SRA deflines are detected as paired.** `@RUN.5 5 length=150` repeats
  the *spot name* on both mates, and its leading digit was read as a Casava mate
  number — so both mates looked like member 5, and since members must be distinct
  that vetoed the pairing and the stream was archived single-end, leaving `--split`
  unable to restore the mate files. A leading description digit now counts as a mate
  field only when followed by `:` (`1:N:0:…`), which is what distinguishes it from a
  numeric spot name. This is the defline the documented
  `sracha get -Z … | fqxv compress -` pipeline produces.
- **Out-of-step inputs warn instead of silently mis-pairing.** With `G > 1` member
  assignment is positional, so a stream missing one member mid-file shifts every
  later read into the wrong mate file while the read count stays a clean multiple of
  `G` — invisible to the spot-multiple check and to every CRC. Compression now
  compares the members' names within each spot and warns with the count and the
  first offending spot. This is the shape `fasterq-dump --split-files` produces for a
  spot whose `READ_LEN` slot is 0: it writes no record at all for the empty slot.

### Added

- Container-level tests for the FASTQ shapes SRA conversion produces: a zero-length
  read inside a `G > 1` spot surviving `--split` in its own slot, and `.` plus the
  full IUPAC ambiguity set and lowercase bases round-tripping (SRA's 4na→text map is
  `.ACMGRSVTWYHKDBN`, so converted FASTQ carries more than `N`).

## [0.6.0] - 2026-07-31

### Added

- **Golden archive fixtures pin the on-disk format.** `crates/fqxv/tests/fixtures/`
  now holds one small archive per on-disk layout — plain/order-k, grouped (with the
  header extension record), the whole-file reorder layout, and the long-read overlap
  codec — plus a manifest of what each must decode to. Until now *every* test in the
  workspace compared a build against itself: the round-trips, the proptests, the
  fuzz targets, the thread-determinism byte-comparisons, and the accession corpus
  (which recompresses from FASTQ each run) are all invariant to a change that moves
  the encoder and decoder together, which is exactly the change that breaks archives
  already written. Regenerate with
  `cargo run --release -p fqxv --example make_fixtures -- crates/fqxv/tests/fixtures`,
  but note that regenerating is almost always the wrong response to a failure.
- **CI checks cross-release compatibility.** A new `compat` job downloads the
  previous release's CLI, compresses with it and decompresses with the PR build
  (a hard failure if that breaks), then does the reverse and requires the old binary
  to either reproduce the bytes exactly or refuse loudly — a silent success with
  different output fails the job.
- **The changelog is on the documentation site.** Included verbatim from this file,
  so the two cannot drift.
- **A documented format evolution policy** (`docs/design/container.md`): which
  mechanism a given change belongs in, the rule that any change to the footer's
  shape or the block payload's stream layout must be gated by a `required_features`
  bit, and the project's intent on major bumps (a last resort; a format-1 read path
  is retained if one ever happens).

### Changed

- **Three implicit parts of the layout now fail closed.** Each was a structural
  convention living in a constant rather than on disk, and each would accept an
  archive from a future writer and get it wrong:
  - A block payload carrying a **fourth stream** decoded its first three and dropped
    the rest with no error at all — and the per-block content digests, which cover
    only the streams that were read, still matched, so the archive looked healthy.
    Trailing bytes after the three streams are now rejected.
  - A **footer written at a different per-group stride** (a future fourth stream
    location) passed its CRC and was then walked at the wrong stride, producing an
    index that could satisfy every range check: `fqxv info` reported a 1,000-read
    archive as 46,048 reads. The footer body length must now match its row-group
    count exactly.
  - **Unknown header flag bits** (6 and 7) were ignored, so a future flag whose
    meaning is purely semantic would decode under the old interpretation with every
    CRC and digest still matching. Unknown bits are now refused with a new
    `Error::UnsupportedFlags`.

  All three are backward compatible: no archive any release has written trips them.
- **Original per-slot member labels are recorded and restored.** Member identity in
  the container is positional, so `decompress_split` could only number its outputs
  `_R1.._RG` — compressing a run's `_2.fastq` + `_4.fastq` and restoring them
  renamed them `_R1`/`_R2`, and a 10x `I1` came back as `_R3`. The CLI now derives
  each member's slot token from the input file names (all-or-nothing: every input
  must yield a distinct token, or nothing is recorded) and stores them in the
  header. `--mate-style auto` (the new default) restores the stored labels, falling
  back to `_R1`,`_R2`,…; explicit `r`/`num` stay positional. This is the **first
  user of the header extension region**: a non-critical TLV record (tag `0x01`), so
  a reader that predates the tag skips it and decodes byte-identical reads under
  positional names. The format version deliberately stays at **1.0** — the record is
  skippable, and nothing about decoding depends on it. Archives compressed without
  labels still write an empty extension region.

### Fixed

- **`--estimate`'s help text described the wrong mechanism.** It claimed to code the
  leading reads with the real codecs; it measures the sample's empirical entropy and
  projects from that, which is why it is so much faster than a real run. The
  documentation was right and the shipped `--help` string was not. A documentation
  audit against the actual CLI, bindings, and codecs corrected this and a good deal
  more — the `fqxv info` sample output was internally impossible, `-v` was
  documented as debug when it is info, batch mode over several archives or a
  directory was undocumented, and the container layout diagram omitted the
  whole-file reference frame that sits between the header and the first block.
- **`fqxv.Index.groups()` no longer panics** on an allocation failure; it raises.
- **Unequal input read counts are rejected in both directions.** `compress_multi`
  treated member 0's EOF as a clean end of input, so a *short* member 0 ended the
  archive early and silently dropped the surplus reads from the longer members —
  yielding a well-formed, CRC-valid archive that was quietly missing data and still
  reported success (with an inflated ratio, since it is computed against the full
  input size). All three interleaving sites now confirm members 1..G are also spent.
  The reverse direction already errored.
- **A trailing partial spot is rejected in the reorder layout too.** The plain
  layout refuses an interleaved stream whose record count is not a multiple of the
  group size, but the reorder layout only checked it at one call site, which
  `compress_auto` bypasses — so `--order any` / `--max` on a mate-named stream with
  an odd record count recorded `group_size = 2` over a stream that is not
  spot-aligned and split it into mismatched files. The check moved into
  `encode_reordered`, the choke point every reorder entry point passes through.
- **The streaming parser validates the `+` separator line.** `read_raw_record`
  compared only the sequence and quality lengths, so a header line followed by EOF
  (`"@name\n"`) parsed as a zero-length record — silently repairing a truncated file
  into a valid archive — and a garbage third line was consumed and discarded. Both
  were reachable only through the multi-input path (a single-file compress rejects
  them at the noodles-based peek). The `+` line's *content* is still dropped; `+`
  normalization is unchanged.
- **A non-empty header extension region no longer corrupts the footer.**
  `FooterIndex::new` seeded block offsets from the extension-empty header length and
  the recovery scan started there, so any extension record shifted every recorded
  offset and failed the footer CRC. Both now use the header's actual length. The
  region had never been exercised before the member-label record.

## [0.5.2] - 2026-07-29

### Added

- **`fqxv.estimate()` accepts paired-end (and other grouped) inputs.** Pass a list
  or tuple of sources — paired mates, 10x R1/R2/I1/I2, sharded lanes — and each is
  sampled with a split of the `sample_reads` cap, with the per-stream sizes summed
  into one aggregate `Estimate`, matching `fqxv compress … --estimate` byte for
  byte. Single-source behavior is unchanged. `Estimate` also gains a `platform`
  field (`"illumina"`, `"nanopore"`, `"pacbio"`, `"mgi"`, `"unknown"`).
- **Cross-platform guard on grouped estimates.** A group whose inputs resolve to
  two different *known* platforms is now a hard error in both the Python API and
  the CLI's `--estimate`, instead of silently summing two differently-calibrated
  numbers. An `Unknown` input (SRA-renamed, no content signal) stays compatible
  with anything, so real SRA data never trips a false positive.

### Changed

- **New `fqxv-align` crate.** The two edit-distance alignment implementations —
  banded Needleman-Wunsch (AVX2 + scalar) and the wavefront aligner — moved out of
  `fqxv-lroverlap` into a zero-dependency leaf crate, so an external consumer can
  reach them without pulling in the whole long-read codec cone. The move is a pure
  rename: `fqxv-lroverlap` re-exports all seven public names at their existing
  paths, and its public surface is unchanged. `fqxv-lroverlap` is now
  `#![forbid(unsafe_code)]`, since the AVX2 backend that motivated the weaker
  `deny` has left the crate.
- Dependency bumps: `noodles-bgzf` 0.48 → 0.49, `clap` 4.6.2 → 4.6.3,
  `xxhash-rust` 0.8.17 → 0.8.18, `serde_json` 1.0.150 → 1.0.151.

## [0.5.1] - 2026-07-23

### Performance

- **`fqxv compress --estimate` is now 15–60× faster** — always subsecond on a
  block-sized sample. It predicts the archive from the sample's measured *entropy*
  (what the real entropy coders converge to) instead of coding it: static order-k
  sequence entropy, a k-mer duplication sketch that captures long-read cross-read
  redundancy, order-1 quality entropy, and the real tokenizer for names — each
  scaled by a small per-platform calibration factor. The measurement passes run in
  parallel across read chunks, and the sample is capped by read count for short
  reads and by one block of bases for long reads (so the sketch sees a real block's
  coverage). Predicted archive size stays within ~1% of the previous coding-based
  estimate across Illumina, PacBio HiFi, and Nanopore.

## [0.5.0] - 2026-07-23

### Added

- **`fqxv decompress -` reads from stdin,** streaming the archive straight into the
  decoder (it only reads forward and stops at the terminator frame). This is how you
  read a remote archive — pipe a transfer tool in, so all the auth, retries, and
  resume stay in the tool built for it and fqxv gains no HTTP dependency:
  `aws s3 cp s3://bucket/reads.fqxv - | fqxv decompress - -Z | bwa mem ref.fa -`
  (or `curl` on a presigned URL, `gsutil cat`, …). A truncated stream still fails
  (premature EOF) rather than yielding a short file. `--recover` and `--split` need
  a seekable file (a stream can't be rewound) and refuse stdin with a clear message.
- **Python: read archives over the network.** The `fqxv` wheel gains a
  dependency-free `fqxv.remote` module (standard-library `urllib`). The streaming
  entry points — `fqxv.open`, `decompress_to_*` — now accept **any file-like
  object**, so a `boto3`/`urllib`/`httpx` response streams straight in
  (`fqxv.open(s3.get_object(...)["Body"])`); `fqxv.remote.stream(url)` /
  `download(url, dest)` wrap that for a URL. `fqxv.remote.RemoteArchive` adds
  **column projection** over HTTP byte-range requests: fetch just the footer index
  from the archive tail, then only the names (~1% of the file) or the sequence,
  CRC-verified per stream. It rests on IO-free primitives —
  `parse_index_suffix`, per-stream ranges/CRCs on `Index`, and
  `decode_{names,sequences,qualities}_bytes` — that a custom async client can drive
  directly for concurrent range fetches.

### Performance

- **Reorder compress: SIMD merge scan, malloc-free placement, incremental
  `cast_vote`.** Three byte-identical, zero-ratio-cost speedups on the hot
  short-read reorder assembly/merge phases (`--order any`/`shuffle`/`--max`).
  (1) The overlap-merge successor scan now compares 16 bytes at a time
  (`_mm_cmpeq_epi8` + movemask popcount, budget checked per block) behind a
  runtime `is_x86_feature_detected!("sse2")` dispatch with a byte-identical
  scalar fallback (no raised global baseline) — the count matches the scalar
  loop exactly when within budget and exceeds it exactly when the scalar loop
  would break, so `best_key` and the archive are unchanged. (2) `place_on_contig`
  and the rescue assembler's `try_place`/`place` no longer allocate a
  `Vec<usize>` of mismatch positions per candidate: candidates are scored with a
  count-only pass and only the winner's positions are materialized (into a
  caller-reused scratch buffer on the clustered path). (3) `cast_vote` updates
  the per-column plurality in O(1) by comparing the just-incremented base against
  the current winner (lowest-index tie rule preserved) instead of `max_by_key`
  over all four counts. About **2% faster single-thread** (54.87 s → 53.83 s,
  1.02× ± 0.00, on 1M-read NovaSeq `--order shuffle` compress); within noise at
  `--threads 16` (1.00× ± 0.04), where the wall clock is dominated by the serial
  assembly prelude. **Byte-identical** archives on the validation datasets
  (NovaSeq/GAIIx/MiSeq/MGI, both `--order shuffle` and `--max`), and
  thread-deterministic (`--threads 1` == `--threads 16`). (#219)
- **Skip the quality-quantizer probe on Nanopore.** The per-block quantizer trial
  (`MODE_SEQ_BINMIX_Q`) only wins on skewed PacBio HiFi/Revio quality; on the flatter
  Nanopore distribution it never does, yet the bounded prefix probe still ran and cost
  ~5% of ONT compress for 0 bytes saved. `encode_seq` now takes a `try_quantizer` flag
  and the container passes `false` on Nanopore, skipping the histogram build and the
  probe entirely. **Output is byte-identical** (the kept stream was the baseline either
  way); only ONT compress speed changes. HiFi/Revio still trial and keep the quantizer.

### Internal

- **Optimize the workspace codec crates in dev/test builds.** The `[profile.dev]`
  `package."*"` override optimized external dependencies but not workspace members, so
  the codec crates compiled at opt-level 0 for tests — and the test suite runs real
  compression, so the long-read round-trip tests took 30–45 s each. Per-package
  `opt-level = 2` overrides for the compute crates cut the CI test suite roughly 4–5×
  with no behavior change (SIMD is runtime-detected; determinism holds across opt-levels).
  The CLI stays at opt-level 0 to keep the edit-compile loop fast.

- **Byte-identical speed in the long-read overlap search (chaining).** The minimizer
  overlap search (`fqxv-lroverlap`) is the top self-cost of Nanopore compress, and a
  overlap search (`fqxv-lroverlap`) is the top self-cost of Nanopore compress, and a
  query's per-anchor chaining dominates it (a fresh profile put `Chainer::chain` at
  ~23% and `find_overlaps` at ~14% of ONT compress). Three output-preserving changes
  strip wasted work without moving a single archive byte:
  - **The chainer's redundant re-sort is dropped.** `find_overlaps` already sorts a
    query's whole anchor set by `(target, strand, tpos, qpos)`, so each
    `(target, strand)` group it hands the chainer arrives already in `(tpos, qpos)`
    order — exactly what the chainer's own `sort_unstable` would produce. A new
    `Chainer::chain_presorted` entry point skips that re-sort (it still dedups),
    saving hundreds of sorts per read at high coverage.
  - **The chaining DP's scratch buffers are reused.** The per-anchor score,
    predecessor, order, and used arrays — plus the per-chain path — were
    heap-allocated on every group (hundreds of small allocations per read). They are
    now thread-local buffers, cleared and resized per group, never read stale — so
    the chain set stays a pure function of the input.
  - **The hot anchor-bucket sort is a radix sort.** Each anchor's group-and-position
    key `(target, strand, tpos, qpos)` packs losslessly into a `u128` whose
    ascending order is exactly the tuple order, so the per-query sort (~7% of ONT
    compress) is now a deterministic LSD radix sort producing the identical total
    order rather than a comparison sort.
  About **15% faster single-thread** ONT compress (`ecoli_ont`/DRR205413,
  559 s → 475 s), where chaining is the largest self-cost. At 16 threads long-read
  compress is block-bound — a handful of large blocks gate the wall clock, not the
  per-anchor work — so the gain is dataset-shaped: ~noise on ONT (few, large blocks)
  but about **10% faster** on high-coverage HiFi (`ecoli_hifi`, `--platform pacbio`,
  254 s → 228 s at 16 threads), where the many smaller blocks keep the cores fed.
  Archives are **byte-identical** on ONT (default and `--max`) and HiFi, and
  thread-count invariant (`--threads 1` == `--threads 16`). Proptests pin the radix
  order to the comparison sort's and the presorted chain path to the sorting path.
  The levers mirror minimap2's own presorted chaining and `radix_sort_128x`. (#151)
- **Per-block quality context quantizer, trialled and kept only when it wins
  (HiFi/Revio).** The long-read quality coder used to build its recent-quality
  context with one fixed quantization (`q1>>1`/`q2>>3`/`q3>>4`), which merges
  adjacent Phred values — costly on HiFi/Revio, where quality is packed at the top
  of the scale and neighbouring high-Q values need to be told apart. The encoder now
  also builds a per-block quantizer from the block's quality histogram (the fqzcomp
  `qtab` / CoLoRd platform-quantizer idea): full context resolution where quality
  actually varies, equal-population folding elsewhere. Both quantizers code the
  block and the **smaller** is kept, with the choice recorded in a self-describing
  header mode byte (`MODE_SEQ_BINMIX_Q`) and the small table transmitted so decode is
  unambiguous — so a block can only match or shrink (never-worse by construction).
  Quality-stream savings, measured against a clean baseline: **PacBio HiFi (ecoli,
  Sequel II) −1.17%**, **Revio amplicon −3.25%**, **Revio WGS −0.25%**; **Nanopore
  0%** and **short reads 0% (byte-identical)** — neither is touched. A bounded
  prefix probe decides whether the second full encode is worth coding, so the trial
  costs ≈ +5% compress on Nanopore (where it never wins) and is paid in full only on
  the skewed long-read data where it does. Lossless and thread-deterministic;
  archives round-trip byte-for-byte and are identical regardless of thread count.

- **Anchor-restricted long-read tile coding (CoLoRd-style).** The multi-reference
  ONT tiler used to re-align each tile with one banded DP over the whole
  `read × reference` window, re-deriving the exact-match stretches its own minimizer
  chain had already proven identical. It now walks that chain, emits each anchor as a
  free `Match` copy-run, and runs the aligner only on the short inter-anchor gaps and
  the two flanks — so alignment work scales with divergence, not read length. About
  **2.5× faster** single-thread ONT tiler compress, with an aggregate ratio
  improvement across the ONT corpus (most accessions smaller — e.g. DRR205413 −2.9%,
  DRR351396 −2.8%, DRR424350 −2.6%; worst real-data case ≈ +0.3%). Lossless and
  thread-deterministic; the edit-op stream keeps the same semantics, so the decoder
  is unchanged and archives still round-trip byte-for-byte. A tile whose chain is not
  recovered falls back to the whole-window DP; `FQXV_TILE_NO_ANCHORGAP=1` forces the
  DP everywhere for A/B. (#226)

- **Anchor-restricted coding for the shared-reference consensus codec (HiFi).** The
  same lever as #226, now on the non-tiler path the high-coverage HiFi/PacBio blocks
  take: each read used to be re-aligned against its consensus with one banded DP over
  the whole `read × consensus` window (`align_banded`, the dominant HiFi compress
  self-cost). It now recovers the read↔consensus exact-k-mer chain, emits each anchor
  as a free `Match` copy-run, and aligns only the short inter-anchor gaps and two
  flanks — so alignment work scales with the read's divergence from its consensus,
  not its length. On low-error HiFi the read shares nearly all of its k-mers with the
  consensus, collapsing the DP area. About **1.23× faster** at-scale (16-thread)
  full-file `ecoli_hifi` compress (313s → 254s, two trials), at an equal-or-slightly-better ratio
  (−0.05% archive on the full file; byte-identical on a single-block subset, where the
  anchor path emits the same equal-cost edit stream). Lossless and thread-deterministic
  — the edit-op stream keeps the same semantics, so the decoder is unchanged and
  archives still round-trip. A read whose chain is not recovered falls back to the
  whole-window DP;
  `FQXV_LRO_NO_ANCHORGAP=1` forces the DP everywhere for A/B. The `anchor_chain` /
  `anchorgap_build` machinery is now shared by both coders (new `anchorgap` module).
  (#220)

- **Reorder drops the redundant v2 single-contig sequence candidate per block.** On
  the adaptive `rescue` path (`--order any`/`--max` default) each block coded both
  the v2 single-contig codec and the v3 literal-rescue codec and kept the smaller.
  v3 generalizes v2 — it attaches the reads v2 strands as literals and otherwise
  degenerates to the same coding — so it is the block-local floor on its own;
  coding v2 too was near-redundant. v2 is now coded only under `--no-rescue` (its
  intended fast single-contig path). Output is **byte-identical** on the validation
  datasets (NovaSeq/GAIIx/MiSeq/MGI, `--order shuffle` and `--max`) — v3 was never
  larger than v2 on any block — for one fewer per-block sequence encode
  (~1.08× faster single-thread `--max` on a 2M-read NovaSeq subset).

- **Single-end reorder codes quality and names once, not twice.** The never-worse
  gate on `--order any`/`--order shuffle`/`--max` used to code the whole file both
  clustered and plain and keep the smaller *full* archive — coding the quality and
  name streams a second time. The keep/skip decision is really "did clustering +
  the permutation beat plain order-k *on the sequence*", so the gate now codes each
  candidate's sequence + names + permutation, decides on that non-quality total,
  and only the winning layout codes quality — once, in the order it emits (the
  clustered encode is split into a prepare/finish pair; the plain candidate reuses
  its coded names+sequence via a new `write_plain_layout` path). Output is
  **byte-identical** to the previous build on the validation datasets
  (NovaSeq/GAIIx/MiSeq/MGI, both `--order shuffle` and `--max`) and
  thread-deterministic; it removes a redundant whole-file quality/name encode
  (a small single-thread win, since quality coding is cheap next to the reorder
  assembly and the still-required plain sequence coding).

- **Lower peak memory on Nanopore compression — the whole-input buffer is gone.**
  `compress_auto` routed every long-read input to a buffered path that held the
  entire file in memory to build the whole-file shared reference (#168). #211
  disabled that reference on high-error Nanopore, leaving the buffer as dead
  weight there. Explicit-Nanopore single-end input now takes the streaming path
  instead (one block in flight, not the whole file): heap-profiling put the buffer
  at ~25% of peak, and on a 600 MB ONT file peak RSS drops **3.72 GB → 2.91 GB**
  (−22%). The archive is **byte-identical** — streaming cuts blocks at the same
  `MAX_BLOCK_SEQ_BYTES` budget and, with no shared reference, each Nanopore block
  codes with the same plain per-block codec either way (verified by `cmp` and the
  determinism round-trips). HiFi still buffers (its shared reference pays off); an
  auto-detected ONT set without `--platform` also still buffers, since the
  streaming path detects platform from read names (issue #225).
  
### Changed

- **Long-read compress is faster at identical output.** Two changes cut compress
  time without moving any ratio or byte of the archive: the overlap-consensus
  candidate that Nanopore always discarded is no longer built (#223, ~42% faster
  ONT default), and the banded-DP traceback is de-packed to one byte per cell
  (#222, ~25% faster ONT / ~13% faster HiFi). Combined, default-mode ONT compress
  is ~2.3× faster than v0.4.0. A full benchmark rerun on this build reproduced
  every ratio and per-stream size byte-for-byte.

### Added

- **Python bindings expose `estimate()` and `verify()`** (#218).

## [0.4.0] - 2026-07-21

### Added

- **Best-of-N reference selection takes the ONT tiler to CoLoRd parity** — the
  multi-reference tiling codec (`SEQ_METHOD_TILE`) now weighs several earlier-read
  references per tile and keeps the cheapest edit script, instead of blindly taking
  the furthest-reaching neighbour. At ONT coverage many earlier reads span the same
  region with independent error patterns, so the lowest-cost reference agrees with
  the query at more positions — the way CoLoRd picks its anchor, applied per tile.
  With a wider alignment band on top this is the dominant ONT sequence-ratio lever:
  on the `ecoli_ont` archive the sequence stream drops from **43.3 MB to 34.5 MB
  (1.150 → 0.915 bits/base)**, matching CoLoRd's ~34 MB, and the whole archive goes
  **2.92× → 3.05×** — losslessly (`--verify`). It is encoder-only (the tile block
  self-describes, so the decoder is unchanged and single-reference coding is
  byte-for-byte identical to before) and compute-heavy (~linear in the fan-out on
  the tiler's already-vectorised alignment), so it is gated to the top effort
  levels: the default keeps the single-reference cover and today's ONT speed,
  `-l7`/`-l8` enable best-of-2/-4, and `--max` reaches the parity operating point.
  Nanopore long reads only. `FQXV_TILE_BAND` / `FQXV_TILE_REFS` override the band
  and fan-out for rebuild-free A/B measurement.

- **Raw-LZMA sequence codec for ordinary-coverage long reads** — a new per-block
  sequence method (`SEQ_METHOD_LZMA`) codes the block's ASCII bases with a
  large-window LZ, kept only when it beats the existing candidates. At ordinary
  long-read coverage — a real genome, not 300× of one organism — overlapping HiFi
  reads share long *exact* substrings that neither the within-read order-k model
  nor the consensus-edit overlap codec captures, but a large-window LZ finds
  directly. On the full PacBio Revio WGS run the archive goes from **9.72× to
  13.35×** — from second-worst in the lossless field, below plain `zstd19`/`xz9`,
  to above both — losslessly. The overlap codec still wins on
  high-coverage/amplicon data, and the per-block `keep_smaller` gate picks the
  right one automatically. `FQXV_SEQ_NO_LZMA` disables the candidate for A/B
  measurement.

### Changed

- **`--version` reports git provenance on development builds.** A binary built
  from a clean checkout of the release tag still prints just `fqxv 0.4.0`;
  anything else — commits past the tag, a dirty tree, an untagged branch —
  appends the git description (`fqxv 0.4.0 (v0.4.0-7-gab12cd34-dirty)`), so a
  bug report identifies the exact build. Source trees with no git (a crates.io
  tarball) print the plain version.

- **Real progress bars for paired compress and decompress.** Both paths showed a
  bare indeterminate spinner for the whole run — minutes of no information on a
  paired sample at a few cores. Compress picked its indicator on input *count*, so
  paired mates always fell through to the spinner even though the summed input
  length the bar needs was already computed; interleaving reads every input in
  lockstep, so one shared counter against that sum drives the same bar the
  single-input path already used. Decompress had no counter at all and now reports
  reads finished against the footer's read count (already read up front for the
  truncation check) — counting output records, not archive bytes consumed, which in
  `--threads`-sized batches runs ahead of the work and would pin at 100% while the
  last batch is still decoding. Both bars now carry an ETA. The reorder and stdin
  carve-outs still degrade to an indeterminate readout rather than report a
  percentage that would stall or lie, and `--quiet` and non-TTY output are
  unchanged (#189).

### Performance

- **Lower peak memory on Nanopore compression.** The minimizer-occurrence record
  (`Occ`) that dominates the long-read overlap index is packed from 24 bytes to 16
  — the read index and the strand flag now share one `u32` (read in the low 31
  bits, strand in the high bit) instead of a `u32` beside a `bool` that padded the
  struct out to a third 8-byte word. That is a one-third cut of the index's largest
  allocation, on the order of a gigabyte on a large ONT block, so it directly eases
  the Nanopore memory pressure. Output is byte-identical — a hand-written `Ord`
  reproduces the previous `(hash, read, pos, strand)` sort order exactly — and that
  was verified both by the determinism round-trips and by an archive-level diff
  against the previous build.

- **Nanopore compression is ~2.3x faster** — the shared whole-file reference layout
  (#168) is now skipped on high-error Nanopore, where it never pays off. It helps
  low-error PacBio HiFi (a clean consensus stored once), but on noisy ONT the
  reference frame always costs more than it saves, so the whole-file gate (#184)
  rejects it every time — *after* building the whole-file reference and coding every
  block against it, a second full long-read assembly per block that profiling put at
  ~45% of ONT compress CPU, all discarded. Skipping it on Nanopore goes straight to
  the plain layout the gate would fall back to regardless, so the output is
  **byte-identical** — measured 7:20 → 3:11 on a 600 MB E. coli ONT file. HiFi is
  unaffected. Mirrors the Nanopore LZMA gate above.

- **PacBio HiFi is smaller *and* ~2x faster to compress** — the raw-sequence LZMA
  codec (`SEQ_METHOD_LZMA`, #197) now seeds its match finder on **12 bytes**
  instead of 4. On raw ASCII a 4-byte seed keys only ~256 distinct DNA 4-grams, so
  the hash chains collided catastrophically and the depth cap truncated them —
  which both cost time (millions of pointer-chases through a capped chain) and
  *lost matches* beyond the cap. A 12-base seed keys ~16.7M grams, so a 12-mer
  recurs only ~a dozen times: chains are short, the cap never binds, and the
  finder sees every candidate. On a 40k-read Revio WGS subset the sequence stream
  drops **0.791 → 0.638 b/base** and the whole archive **29.4 → 25.2 MB
  (−14%)**, while compress time **falls ~2x** (329s → 166s, below even the
  original 4-byte path). Validated lossless and on a second HiFi dataset. The
  match finder is tuned per call site, so the 2-bit packed-reference path keeps
  its 4-byte seed and byte-identical output; decode is seed-agnostic, so no format
  change. The match-extension inner loop is also AVX2-vectorized (CPU-feature
  gated, byte-identical, ~8% on top; `FQXV_LZMA_NO_SIMD` forces scalar).

- **Nanopore compression is ~5.5x faster** — the raw-LZMA sequence candidate
  (`SEQ_METHOD_LZMA`) is now skipped on high-error Nanopore, where it cannot win.
  LZMA pays off only where reads share long *exact* substrings — ordinary-coverage
  low-error PacBio HiFi (0.60 vs the overlap codec's 1.39 b/base); on Nanopore a
  base error every ~10 bases chops those matches short, so it reliably *loses* to
  the overlap/order-k codecs (2.2 vs 1.79 b/base) while costing the most encode
  time of any candidate. Because it runs on the serial block-0 probe as well as
  the fan-out blocks, coding-and-discarding it dominated the ONT wall-clock: on a
  600 MB E. coli ONT run compress drops **1362s → 246s with byte-identical
  output** (LZMA was ~82% of the time and changed nothing). The gate keys off the
  detected platform, so PacBio still codes LZMA and keeps its win, and any
  platform where the loss can't be ruled out up front (Unknown) still codes it —
  the never-worse ratio guarantee is unchanged. `keep_smaller` continues to floor
  every block regardless.

- **Long-read compression is ~40% faster where the shared reference wins
  outright.** Making the whole-file gate exact meant coding both layouts for
  every block, which costs a second long-read assembly per block — on PacBio
  HiFi that roughly doubled compress time (225s → 461s on `ecoli_hifi`) to
  evaluate a candidate that never had a chance. The encoder now probes only the
  **first** block both ways and, when the reference wins by a wide enough margin
  to cover the reference frame three times over, skips the plain candidate for
  the remaining blocks. `ecoli_hifi` drops 461s → 275s with byte-identical
  output; `ecoli_ont`, where the plain layout genuinely wins, sees a negative
  margin and still runs the exact gate, so the ONT ratio is unchanged. The probe
  is integer arithmetic over block 0, so the decision stays thread-count
  invariant, and a shortcut block is still floored at order-k.

### Fixed

- **Higher effort no longer produces a larger archive.** On a low-redundancy
  library (BGISEQ-500) `--max` and `-l9` compressed *worse* than the default,
  because two effort knobs were applied unconditionally even where they cost more
  than they saved (#196). Both are now never-worse: (1) the level-8+ hashed
  high-order sequence tier is coded alongside plain order-k and the smaller kept,
  so enabling it can only help; (2) single-end `--order any`/`--max` codes the
  clustered and plain layouts and keeps the smaller, so reordering is used only
  when it pays — the same "code both, keep the smaller" rule the long-read shared
  reference (#192) and the overlap-vs-order-k choice already use. On the BGISEQ
  case `--max` drops from a 340 KB regression to matching the default; on a
  reorder-favourable NovaSeq set `--max` keeps its 2.4x win. Grouped (paired /
  10x) reorder is unchanged: it pays the permutation for spot reconstruction
  regardless, so the tradeoff differs.

- **The long-read shared-reference gate now measures against the layout it falls
  back to.** Adoption of the whole-file reference was gated on beating the plain
  *order-k* total, but the fallback path keeps the smaller of the per-block
  overlap codec and order-k — a stronger result than the bar being tested. A
  reference that lost to the per-block overlap codec could therefore still be
  adopted. Both layouts are now coded in the first pass and the smaller wins, so
  the never-worse property holds against the real alternative; the fallback
  reuses those streams instead of re-coding them. No archive changes on the
  benchmark corpus — ONT and HiFi output are byte-identical — so this closes a
  latent hole rather than changing current results (#184). **Correction:** this
  entry originally claimed coding both ways "costs no measurable wall time". That
  was wrong; it came from an A/B in which the build had not picked up the change,
  so two identical binaries were compared. Measured properly, coding the second
  candidate roughly doubles long-read compress time (`ecoli_hifi` 225s → 461s,
  six blocks). The correctness argument for the gate is unchanged, but the cost
  is real and is the price of an exact never-worse comparison.

- **ONT seeding is chosen by index coverage, recovering 2.79 MB.** Closed
  syncmers conserve anchors better than window minimizers at ONT error rates, but
  only once coverage is deep enough for the extra anchors to find partners. The
  long-read encoder builds two indexes on opposite sides of that crossover, and
  using syncmers for both put the ONT archive behind. Measured on `ecoli_ont`
  (block 1, 268 Mbase):

  | index | syncmer | minimizer |
  |---|---:|---:|
  | whole-file reference | **1.280 b/base** | 1.416 |
  | per-block overlap | 1.517 | **1.243 b/base** |

  Each index now takes the scheme that suits its coverage, so both land on their
  better number instead of one giving up ~18%. `ecoli_ont` goes 216,381,249 →
  213,587,759 bytes (ratio 2.791 → 2.827), verified content-lossless. Both
  schemes share `w` and `k`, so anchor density and index cost are unchanged.
  PacBio is deliberately left unsplit — at <1% error minimizers are already
  near-optimal at either coverage — and HiFi output is byte-identical (#184).

- **Interrupting a decompress no longer leaves a partial FASTQ at the
  destination.** Ctrl-C during a decode published truncated output at its real
  name, and a short FASTQ is undetectable — BGZF blocks are self-contained so
  `zcat` exits 0 with no warning, and FASTQ carries no terminator or record
  count, so most of a run could sit there looking like a clean decode. Compress
  was already protected (sibling temp, rename on completion, registered with the
  signal handler); decompress created split mates and `-o FILE` directly at their
  final paths, so the handler had nothing to clean up and `process::exit(130)`
  never unwound to run `Drop`. Decompress outputs now route through the same
  `AtomicOutput` as compress, and the commit is deferred until *after* the footer
  read-count check, so a truncated archive also leaves the destination untouched
  instead of publishing a short file and then erroring. Recover commits too — its
  output is deliberately partial but still a finished artifact. Exit code is still
  130, and a full decode round-trips byte-identically (#188).

### Security

- **Decode-path allocation guards for the long-read overlap and reorder codecs.**
  Several decoders sized allocations straight from untrusted on-disk length fields,
  before the backing bytes were known to exist — the #142 allocation-bomb class the
  tiler and the reorder reference decoder already guard against, ported to the
  siblings that had missed it. In `fqxv-lroverlap`'s consensus overlap decoder
  (sequence methods `SEQ_METHOD_OVERLAP` / `SEQ_METHOD_OVERLAP_REF`) a ~10-byte
  crafted block could pre-reserve a read-count or contig-count vector from a raw
  header varint (a capacity-overflow *abort* under `panic = "abort"`) or drive a
  multi-terabyte `total_bases` zero-fill; it now rejects an implausible
  `total_bases` up front against a shared `MAX_BASES_PER_BYTE` ceiling, leaves the
  count-driven vectors unreserved (each is bounded by the stream it reads from),
  and caps every insertion run at the read length. The reorder clustered / rescue /
  global block decoders bounded each input stream but not the *aggregate*
  reconstructed output, so a few kilobytes of `MATCH`-clone and `CONTIG`-copy ops
  could expand into unbounded memory; they now accumulate output bases against the
  same per-block ceiling the streams use. All of these are reachable from
  `fqxv decompress` on a crafted or corrupt archive (frame CRCs do not help — an
  attacker recomputes them). Decode-only and byte-identical on every valid archive;
  covered by a new regression test and the existing decode round-trips.

## [0.3.0] - 2026-07-19

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

- **The on-disk format is stable at 1.0.** Versioning moved from a single
  monotonic `FORMAT_VERSION` integer to `FORMAT_MAJOR.FORMAT_MINOR`: a reader
  refuses a differing major and tolerates a newer minor, and additive features are
  gated behind required-feature bits so a reader that predates a feature refuses
  the archive outright rather than misreading it. That contract is now a stability
  guarantee — archives written by any format-1.x writer stay readable by later
  releases, and a major bump would be announced as a breaking change. The 1.0 format carries the
  long-read overlap codec, the extended per-stream footer index, and the
  per-stream block digests below.
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
- **`--help` layout** — the flags most runs need are listed first and the
  rarely-touched knobs are grouped under an `Advanced` heading, so the default
  help is short enough to read. `-h` prints a one-line summary per flag and
  `--help` expands the detail (when a knob changes the data or the guarantees,
  the long form says so).
- **Rust edition 2024** — the workspace moved from edition 2021. The MSRV is
  unchanged at 1.95, so this is invisible to downstream builds on a supported
  toolchain.

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

[0.6.0]: https://github.com/rnabioco/fqxv/releases/tag/v0.6.0
[0.5.2]: https://github.com/rnabioco/fqxv/releases/tag/v0.5.2
[0.5.1]: https://github.com/rnabioco/fqxv/releases/tag/v0.5.1
[0.5.0]: https://github.com/rnabioco/fqxv/releases/tag/v0.5.0
[0.4.0]: https://github.com/rnabioco/fqxv/releases/tag/v0.4.0
[0.3.0]: https://github.com/rnabioco/fqxv/releases/tag/v0.3.0
[0.2.0]: https://github.com/rnabioco/fqxv/releases/tag/v0.2.0
[0.1.0]: https://github.com/rnabioco/fqxv/releases/tag/v0.1.0
