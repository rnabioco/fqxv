# fqxv fuzzing

Coverage-guided [`cargo-fuzz`](https://rust-fuzz.github.io/book/cargo-fuzz.html)
targets over every public **decode** entry point. The invariant under test: a
decoder must return `Ok`/`Err` on *any* bytes — never panic or abort (the release
profile is `panic = "abort"`, so a panic is a crash).

This complements the in-tree proptest corruption harness
(`crates/fqxv/tests/corruption.rs`), which runs in normal CI: that harness
encodes valid data and mutates it (reaching deep quickly), while these targets
explore the raw byte space with coverage feedback (finding cases mutation won't).

## Targets

| Target | Entry point(s) |
| --- | --- |
| `container` | `fqxv::{decompress, decompress_recover, inspect, expected_reads}` |
| `rans` | `fqxv_rans::decode` |
| `tokenizer` | `fqxv_tokenizer::decode` |
| `seq` | `fqxv_seq::decode` |
| `fqzcomp` | `fqxv_fqzcomp::decode` |
| `reorder` | `fqxv_reorder::decode_clustered_auto` + `GlobalReference::{decode,decode_blocked,decode_lzma,decode_packed}` |

## Running

Needs a nightly toolchain and `cargo-fuzz`:

```bash
rustup toolchain install nightly
cargo install cargo-fuzz
```

Then, from the repo root:

```bash
export ASAN_OPTIONS=allocator_may_return_null=1
cargo +nightly fuzz run rans -- -malloc_limit_mb=0                     # until a crash
cargo +nightly fuzz run container -- -max_total_time=120 -malloc_limit_mb=0
cargo +nightly fuzz build                                              # just build all targets
```

**The allocation flags are load-bearing.** The decoders use `try_reserve` to
reject a corrupt length gracefully, but that still *asks* the allocator for the
huge size — `-malloc_limit_mb=0` and `ASAN_OPTIONS=allocator_may_return_null=1`
let the allocator return null so the guard can do its job instead of the sanitizer
aborting. A genuine *infallible* over-allocation still aborts via
`handle_alloc_error`, and a hang is caught by libFuzzer's per-input timeout — both
are real crashes.

**Known: decompression-bomb OOMs.** The `seq`, `fqzcomp`, and `reorder` targets
can OOM: a header may declare a large output the *standalone* codec API allocates
for (up to the RSS limit) before erroring. That is resource exhaustion, not a
memory-safety bug — and the **container** (the real trust boundary) bounds it
structurally, so its target is bomb-resistant. `rans` and `tokenizer` have their
count/length-driven allocations bounded by input size, so they are clean too.
Bounding the remaining three the same way (issue #90) is what lets them join the
weekly schedule; until then they run on manual dispatch only.

A crash writes a reproducer under `fuzz/artifacts/<target>/`; replay it with:

```bash
cargo +nightly fuzz run <target> fuzz/artifacts/<target>/crash-<hash>
```

## Seeding the corpus

The `container` target sits behind a per-frame CRC, so raw mutation rarely gets
past the magic + checksum. Seed it with real archives to make it effective:

```bash
mkdir -p fuzz/corpus/container
cp *.fqxv fuzz/corpus/container/          # the sample archives in the repo root
```

The codec targets carry no CRC, so coverage feedback reaches their internals from
an empty corpus; seeding still speeds them up (drop any encoded stream in
`fuzz/corpus/<target>/`).

## CI

`.github/workflows/fuzz.yml` runs the fuzz targets for a short, bounded time. The
**weekly schedule** (Mondays) runs only the verified OOM-clean targets
(`container`, `rans`, `tokenizer`); **manual dispatch** runs all six (for hunting,
where the decompression-bomb OOMs above are expected findings to triage). It is
**non-blocking** — a separate workflow from `ci.yml`, so it never gates a PR. A
crash uploads its reproducer as an artifact and fails that run for visibility.
