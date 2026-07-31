# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`fqxv` is a Rust toolkit for lossless (opt-in lossy) archiving of FASTQ. It's a Cargo workspace of one-crate-per-algorithm codecs plus a container
format library and CLI. All codecs are **clean-room** implementations from specs
and papers (CRAM 3.1 codecs spec, fqzcomp/SPRING/PgRC2), never translated from C
— see `THIRD-PARTY-NOTICES.md`. The on-disk container format is `1.0`
(`FORMAT_MAJOR`.`FORMAT_MINOR` in `crates/fqxv/src/lib.rs`), a number independent
of the crate version (0.5.x): a reader refuses a differing major and tolerates a
newer minor (backward-compatible additions).

## Commands

```bash
cargo nextest run --workspace          # run tests (CI uses --profile ci)
cargo nextest run -p fqxv-rans         # tests for one crate
cargo nextest run -p fqxv-rans decode  # tests matching a substring
cargo test --doc --workspace           # doctests (nextest does NOT run these)
cargo clippy --workspace --all-targets --features fqxv-rans/bench
cargo fmt --all --check
cargo check --workspace                # MSRV is 1.95; keep it building on 1.95
cargo run -p fqxv-cli -- compress reads.fastq.gz -o reads.fqxv
cargo bench -p fqxv-rans               # criterion microbenchmarks
```

**Run every compute command through Slurm — never bare, never in the small
interactive session.** On the Bodhi cluster the login node and the default
`sinteractive` session are tiny (~2 cores); `cargo build`/`check`/`nextest`/
`clippy`/`bench` and any analysis must go to the `rna` partition (88+c/754G) via
`srun`/`salloc`/`sbatch`. Prefer a longer-lived allocation you `srun --overlap`
into over relaunching. Example:

```bash
salloc --no-shell -p rna -c 32 --mem 64G -J fqxv-build   # once
srun --overlap --jobid=<JOBID> -- cargo nextest run -p fqxv   # per command
```

CI (`.github/workflows/ci.yml`) runs check, nextest, doctests, fmt, clippy, and
an MSRV-1.95 check — all with `RUSTFLAGS=-Dwarnings`, so warnings fail the build.
The `check`/`clippy` jobs pass `--features fqxv-rans/bench`; match that when
reproducing CI locally. Each crate keeps its unit tests and `proptest!` blocks
inline in `src/` under `#[cfg(test)]`, and ships a runnable `examples/*.rs`.

The `bench/` directory is a **separate** benchmarking harness (pixi env, Slurm,
`$SCRATCH` data) for comparing releases against the field (fqz_comp, SPRING,
PgRC2, zstd, gzip) — it is not part of the Cargo build. See `bench/README.md`.

## Architecture

The crates form a strict dependency DAG; lower layers never depend on higher
ones. Build understanding bottom-up:

- **`fqxv-bytes`** — leaf crate of on-disk byte primitives shared by every codec:
  LEB128 varints, zig-zag, the `Reader<'a, E>` bounds-checked cursor (generic
  over each crate's `Error` via the `ReaderError` trait), and the read-length
  array codec (`write_lens`/`read_lens`). No dependencies; the single source of
  truth for these encodings.
- **`fqxv-dna`** — leaf crate of nucleotide primitives shared by the sequence
  codecs: the 2-bit ACGT lookup (`BASE_LUT`, `code_strict` case-sensitive /
  `code_fold` case-insensitive, `SYM2BASE`, `base_of_sym`, `is_acgt`) and reverse
  complement (`revcomp`/`revcomp_into` complement both cases; `revcomp_acgt`
  complements uppercase only, passing lowercase through). No dependencies; the
  single source of truth so the previously copy-pasted variants can't drift.
- **`fqxv-rans`** — rANS Nx16 entropy coder (CRAM 3.1). 32 interleaved states;
  order-0/order-1 models. Backends live behind one API and are chosen at runtime
  via `is_x86_feature_detected!`: **scalar** (all orders, the correctness
  reference), **AVX2**, and **AVX-512** — the two vector backends cover order-0
  *both* encode and decode (`avx2.rs`, `avx512.rs`), widest path wins. Order-1 is
  scalar only, and `Backend::Sse42` is a reported CPU tier with no SSE module
  behind it, so it runs scalar too. **Every backend must produce byte-identical
  output** — the cross-backend tests also decode each backend's stream with the
  others, and they self-skip when the CPU lacks the feature, so a green run on a
  non-AVX-512 machine has not actually exercised that path. The `bench` feature
  exposes internal entry points for microbenchmarks only.
- **`fqxv-range`** — serial binary range coder + adaptive bit models. The
  arithmetic-coding primitive that `fqxv-fqzcomp` and `fqxv-seq` build on.
- **`fqxv-fqzcomp`** (→ range) — quality-score context model; owns
  `QualityBinning` (lossless default; `Bin8/Bin4/Bin2` lossy). Re-exported as
  `fqxv::QualityBinning`.
- **`fqxv-tokenizer`** (→ rans) — positional read-name tokenizer with per-column
  delta bucketing; entropy backend is rANS.
- **`fqxv-seq`** (→ dna, range) — order-k adaptive context model over 2-bit ACGT
  symbols (range-coded, variable read lengths); non-ACGT bytes go to an
  exception list. Not a raw 2-bit *packing* path — every base is context-coded.
- **`fqxv-reorder`** (→ dna, rans, seq) — PgRC2/SPRING-class read reordering
  (minimizer clustering, reverse-complement aware) for cross-read redundancy.
  `lib.rs` holds only the crate-common core (`Error`, `IntMap`, decode limits);
  the codec lives in sibling modules (`plan`, `column`, `clustered`, `rescue`,
  `global`, `merge`, plus `refpack`) re-exported flat from the root.
- **`fqxv-align`** — leaf crate of edit-distance alignment primitives, extracted
  from `fqxv-lroverlap` (#252): `align_banded` (banded Needleman-Wunsch, AVX2
  anti-diagonal backend + scalar fallback) and `wfa_align`/`wfa_align_opt`
  (clean-room wavefront, work scales with score not length; `_opt` abandons a pair
  once the score passes a cap). Same unit-edit cost model, so distances agree —
  but the two pick different equal-cost paths, so their edit scripts are *not*
  byte-identical. No dependencies.
- **`fqxv-lroverlap`** (→ align, dna, rans, range, seq) — long-read cross-read overlap codec
  (minimizers → overlaps → layout → consensus → per-read banded edit script →
  rANS). `encode`/`decode` are the container's sequence path for long-read blocks
  (auto-selected, kept only when it beats order-k). Sibling of `fqxv-reorder`;
  never depends on it.
- **`fqxv`** — the `.fqxv` container format; composes all codec crates into
  `compress`/`compress_multi`/`decompress`/`decompress_split`/`inspect`. This is
  where the on-disk layout lives (`src/container/`).
- **`fqxv-cli`** — thin clap front-end over the `fqxv` library (`fqxv` binary).
- **`fqxv-python`** — read-only PyO3 bindings published to PyPI as `fqxv`
  (maturin mixed layout, abi3 wheels). Wraps decode, `inspect`/`estimate`/`verify`,
  and the footer-index projection primitives that `fqxv.remote` drives over HTTP
  range requests. Compression stays in the CLI. Not part of the codec DAG.

The container (`crates/fqxv/src/container/`) is a variable-length header followed
by independent, parallel-codable blocks. The header is a 21-byte fixed prefix, a
TLV extension region (`ext_len` bytes — no longer always empty: tag `0x01` carries
per-member slot labels), and a 4-byte CRC, so anything positioned after it must use
the header's actual length, never a hardcoded 25. Each block splits FASTQ into three
streams handled by three codecs: **names** (tokenizer), **sequence** (order-k seq,
reorder, or the long-read overlap codec — a leading method byte per block picks
one), **quality** (fqzcomp). The exact byte layout is documented in the module doc
comment at the top of `container/mod.rs` — read it before touching the format; the
evolution policy (what warrants a minor bump vs a feature bit vs a major bump) is in
`docs/design/container.md`.

## Invariants to preserve

- **Determinism.** Output must be byte-identical regardless of thread count.
  Blocks are the unit of `rayon` parallelism; keep per-block work order-free.
- **SIMD ≡ scalar.** Any vector backend must match the scalar reference exactly.
  Add proptest round-trips when adding a backend or codec path.
- **Interleaving/spots.** With group size `G > 1`, reads from one spot are
  interleaved (paired mates, 10x R1/R2/I1/I2). Blocks always hold whole spots
  and start on member 0 so a block splits cleanly by `local_index % G`.
- **`+` normalization is intentional** — the optional repeated `+` header line is
  dropped (as SPRING/fqz_comp do); name, sequence, and quality are otherwise
  preserved exactly. This is the one documented deviation from byte-losslessness.
- **CLI effort mapping.** `--level 1-9` maps to sequence context order and block
  size via `level_to_order`/`level_to_block` in `fqxv-cli/src/main.rs`.

## Conventions

- Every crate is independently publishable to crates.io (path + version deps in
  the workspace `Cargo.toml`); dual-licensed MIT OR Apache-2.0. Releasing is
  bump-merge-tag: land a `chore(release): vX.Y.Z` PR (workspace version, the
  internal path-dep pins, `Cargo.lock`, and a dated `CHANGELOG.md` section), then
  tag the squash commit and push. `.github/workflows/release.yml` does the rest —
  it creates the GitHub Release with that changelog section as the notes and
  attaches the CLI binaries; `wheels.yml` publishes the Python wheels. Distribution
  is via bioconda and PyPI, not crates.io.
- Workspace lints forbid `unsafe_op_in_unsafe_fn` and warn on missing docs /
  debug impls / unreachable pub. `cast_possible_truncation` and
  `needless_range_loop` are allowed workspace-wide because the coders rely on
  intentional wrapping arithmetic and parallel-array indexing.
- `.cargo/config.toml` pins `target-cpu=x86-64-v3` for local/bench builds only
  (not for published crates). Do not raise the global baseline for SIMD; use
  `#[target_feature]` + runtime detection like `fqxv-rans`.
- `.cargo/config.toml` also links with `lld` (`-C link-arg=-fuse-ld=lld`). This
  is a no-op on fat-LTO release/bench builds (all the cost is in the compiler)
  but speeds the `cargo nextest run` edit-test loop, which links ~10 per-crate
  test binaries. The toolchain ships lld as `rust-lld`; make `cc` find it once
  with `ln -sf "$(rustc --print sysroot)"/lib/rustlib/*/bin/gcc-ld/ld.lld
  ~/.local/bin/`. CI `apt-get install`s the system `lld` in every compiling job
  — the flag also applies to proc-macro/build-script links, so all of them (not
  just the test jobs) need it.
- Release profile is fat-LTO, single codegen unit, `panic = "abort"`. Use the
  `profiling` profile for samply/perf (keeps symbols, no LTO).

## Profiling this codebase

Hard-won on the reorder audit (#113); re-deriving these costs a run each.

- **Profile the `profiling` profile, never release.** `perf --call-graph=dwarf`
  aborts on the fat-LTO release binary (thin/stripped debug) but works fine on
  `cargo build --profile profiling` (`debug=2`, `strip=false`, full
  `.debug_frame`) — exit 0, clean stacks.
- **`--call-graph=fp` is worse than useless here.** The hot loops clobber `rbp`,
  so the frame-pointer walk reads *sequence data* as return addresses. If a call
  graph shows addresses like `0x4341544141434141` (that is the ASCII `ACAATACA`),
  discard it. This is what manufactured the phantom "22% unattributed".
- **`call_mut` is not a mystery, it is the rayon closure shim.**
  `core::ops::function::impls::<impl FnMut for &F>::call_mut` is where any fully
  inlined closure-body code gets billed. Decompose it by srcline rather than
  treating it as one symbol — half of it was hashbrown's get-side probe.
- **Trust perf's own `--full-source-path` srcline resolution** over hand-rolled
  `addr2line`: against a running PIE that needs the exact mmap base from
  `perf script --show-mmap-events`, and guessing it silently resolves hot samples
  into cold symbols (it once "found" 9% in `tabled::Table::fmt`).
- LBR is unavailable on this cluster (`--call-graph=lbr` exits 255) — the PMU is
  virtualised. Use `-e cpu-clock -F 999`.
- Use hyperfine for wall-clock deltas, and **verify the binary you are timing
  actually contains your change** (`strings`, or a probe log line). Cargo can
  skip a rebuild after a `git stash pop`, and an A/B of two identical binaries
  reports a confident, meaningless zero — that mistake shipped a 2x long-read
  encode regression as "no measurable cost".
