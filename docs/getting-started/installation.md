# Installation

## Prerequisites

- Rust 1.95 or later (the workspace MSRV)
- Cargo (comes with Rust)

## Building the CLI

```bash
git clone https://github.com/rnabioco/fqxv.git
cd fqxv
cargo build --release
```

The binary is at `target/release/fqxv`. Copy it onto your `PATH`:

```bash
cp target/release/fqxv ~/.local/bin/
```

Verify:

```bash
fqxv --version
fqxv --help
```

## Using the crates

`fqxv` is a Cargo workspace of one-crate-per-algorithm codecs plus the `fqxv`
container library. Depend on whichever layer you need:

```toml
[dependencies]
# the whole archiver (container + all codecs)
fqxv = { git = "https://github.com/rnabioco/fqxv.git" }

# or an individual codec
fqxv-rans     = { git = "https://github.com/rnabioco/fqxv.git" }  # rANS Nx16
fqxv-range    = { git = "https://github.com/rnabioco/fqxv.git" }  # range coder
fqxv-fqzcomp  = { git = "https://github.com/rnabioco/fqxv.git" }  # quality model
fqxv-seq      = { git = "https://github.com/rnabioco/fqxv.git" }  # sequence model
fqxv-tokenizer= { git = "https://github.com/rnabioco/fqxv.git" }  # read-name tokenizer
fqxv-reorder  = { git = "https://github.com/rnabioco/fqxv.git" }  # read clustering
fqxv-bytes    = { git = "https://github.com/rnabioco/fqxv.git" }  # shared byte primitives
```

(`fqxv-bytes` is a leaf crate of the LEB128/zig-zag primitives the codec crates
share; the codecs pull it in transitively, so you rarely depend on it directly.)

Every crate is dual-licensed **MIT OR Apache-2.0**.

## Development

```bash
cargo test --workspace          # unit + property tests
cargo test --doc --workspace    # doctests
cargo clippy --workspace --all-targets
cargo fmt --all
```

Benchmarks (against gzip / zstd / xz / fqz_comp / fqzcomp5 / SPRING) live under
`bench/` and use [pixi](https://pixi.sh); see the repository `bench/README.md`.
