# Installation

## Prebuilt binaries

Every release attaches a static `fqxv` binary per platform to its [GitHub
Release](https://github.com/rnabioco/fqxv/releases), plus a `SHA256SUMS.txt`:

| Asset | Platform |
| --- | --- |
| `fqxv-vX.Y.Z-x86_64-unknown-linux-musl.tar.gz` | Linux x86-64 (static, any distro) |
| `fqxv-vX.Y.Z-aarch64-unknown-linux-musl.tar.gz` | Linux arm64 (static) |
| `fqxv-vX.Y.Z-x86_64-apple-darwin.tar.gz` | macOS Intel |
| `fqxv-vX.Y.Z-aarch64-apple-darwin.tar.gz` | macOS Apple silicon |
| `fqxv-vX.Y.Z-x86_64-pc-windows-msvc.zip` | Windows x86-64 |

```bash
VER=v0.6.0   # the latest release tag
curl -LO https://github.com/rnabioco/fqxv/releases/download/$VER/fqxv-$VER-x86_64-unknown-linux-musl.tar.gz
tar xzf fqxv-$VER-x86_64-unknown-linux-musl.tar.gz
mv fqxv ~/.local/bin/
```

The binaries are built for each target's generic baseline; `fqxv-rans` picks its
AVX2/AVX-512 paths at runtime, so one binary runs on old and new CPUs alike.
`fqxv` is not yet on [bioconda](https://bioconda.github.io/) — until it is, use
these binaries or build from source below.

## Prerequisites (building from source)

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

Or install it into `~/.cargo/bin` without keeping a checkout:

```bash
cargo install --git https://github.com/rnabioco/fqxv fqxv-cli
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
fqxv-lroverlap= { git = "https://github.com/rnabioco/fqxv.git" }  # long-read overlap codec
fqxv-align    = { git = "https://github.com/rnabioco/fqxv.git" }  # banded alignment / WFA
fqxv-bytes    = { git = "https://github.com/rnabioco/fqxv.git" }  # shared byte primitives
fqxv-dna      = { git = "https://github.com/rnabioco/fqxv.git" }  # shared nucleotide primitives
```

(`fqxv-bytes` and `fqxv-dna` are leaf crates of the LEB128/zig-zag and 2-bit
ACGT/revcomp primitives the codec crates share; the codecs pull them in
transitively, so you rarely depend on them directly.)

The crates are not published to crates.io — distribution is the CLI binaries
above and the Python package below — so depend on them by git.

Every crate is dual-licensed **MIT OR Apache-2.0**.

## Python

A read-only Python package reads `.fqxv` archives directly — see the
[Python API](../python/index.md):

```bash
uv pip install fqxv
```

## Development

```bash
cargo nextest run --workspace   # unit + property tests (CI uses --profile ci)
cargo test --doc --workspace    # doctests (nextest does not run these)
cargo clippy --workspace --all-targets --features fqxv-rans/bench
cargo fmt --all
```

CI runs the same set with `RUSTFLAGS=-Dwarnings`, plus a build against the 1.95
MSRV.

Benchmarks (against gzip / zstd / xz / fqz_comp / fqzcomp5 / SPRING / CoLoRd)
live under `bench/` and run in the `bench` [pixi](https://pixi.sh) environment
declared in the root `pixi.toml` (`pixi install -e bench`); see the repository
`bench/README.md`.
