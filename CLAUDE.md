# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

`fqxv` is a Rust toolkit for lossless (opt-in lossy) archiving of FASTQ. It's a Cargo workspace of one-crate-per-algorithm codecs plus a container
format library and CLI. All codecs are **clean-room** implementations from specs
and papers (CRAM 3.1 codecs spec, fqzcomp/SPRING/PgRC2), never translated from C
— see `THIRD-PARTY-NOTICES.md`. Status is early development; `FORMAT_VERSION` is
`2` and nothing on disk is stable yet (each build reads only its own version).

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

- **`fqxv-rans`** — rANS Nx16 entropy coder (CRAM 3.1). 32 interleaved states;
  order-0/order-1 models. Backends live behind one API and are chosen at runtime
  via `is_x86_feature_detected!`: **scalar** (all orders, the correctness
  reference) and **AVX2** (order-0 decode only). SSE4.2 and vectorized
  order-1/encode are unimplemented — `Backend` reports the detected CPU tier but
  anything past AVX2 order-0 decode runs the scalar path. **Every backend must
  produce byte-identical output.** The `bench` feature exposes internal entry
  points for microbenchmarks only.
- **`fqxv-range`** — serial binary range coder + adaptive bit models. The
  arithmetic-coding primitive that `fqxv-fqzcomp` and `fqxv-seq` build on.
- **`fqxv-fqzcomp`** (→ range) — quality-score context model; owns
  `QualityBinning` (lossless default; `Bin8/Bin4/Bin2` lossy). Re-exported as
  `fqxv::QualityBinning`.
- **`fqxv-tokenizer`** (→ rans) — positional read-name tokenizer with per-column
  delta bucketing; entropy backend is rANS.
- **`fqxv-seq`** (→ range) — order-k adaptive context model over 2-bit ACGT
  symbols (range-coded, variable read lengths); non-ACGT bytes go to an
  exception list. Not a raw 2-bit *packing* path — every base is context-coded.
- **`fqxv-reorder`** (→ rans, seq) — PgRC2/SPRING-class read reordering
  (minimizer clustering, reverse-complement aware) for cross-read redundancy.
- **`fqxv`** — the `.fqxv` container format; composes all codec crates into
  `compress`/`compress_multi`/`decompress`/`decompress_split`/`inspect`. This is
  where the on-disk layout lives (`src/container.rs`).
- **`fqxv-cli`** — thin clap front-end over the `fqxv` library (`fqxv` binary).

The container (`crates/fqxv/src/container.rs`) is a 10-byte header followed by
independent, parallel-codable blocks. Each block splits FASTQ into three streams
handled by three codecs: **names** (tokenizer), **sequence** (seq or reorder),
**quality** (fqzcomp). The exact byte layout is documented in the module doc
comment at the top of `container.rs` — read it before touching the format.

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
  the workspace `Cargo.toml`); dual-licensed MIT OR Apache-2.0. Releases go
  through release-plz (`.github/workflows/release-plz.yml`, gated on
  `CARGO_REGISTRY_TOKEN`).
- Workspace lints forbid `unsafe_op_in_unsafe_fn` and warn on missing docs /
  debug impls / unreachable pub. `cast_possible_truncation` and
  `needless_range_loop` are allowed workspace-wide because the coders rely on
  intentional wrapping arithmetic and parallel-array indexing.
- `.cargo/config.toml` pins `target-cpu=x86-64-v3` for local/bench builds only
  (not for published crates). Do not raise the global baseline for SIMD; use
  `#[target_feature]` + runtime detection like `fqxv-rans`.
- Release profile is fat-LTO, single codegen unit, `panic = "abort"`. Use the
  `profiling` profile for samply/perf (keeps symbols, no LTO).
