# Installation

[![Bioconda](https://img.shields.io/conda/vn/bioconda/fqxv?label=bioconda)](https://anaconda.org/bioconda/fqxv)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/fqxv?label=downloads)](https://anaconda.org/bioconda/fqxv)
[![PyPI](https://img.shields.io/pypi/v/fqxv)](https://pypi.org/project/fqxv/)

## Bioconda

`fqxv` is packaged on [Bioconda](https://bioconda.github.io/), which is the
easiest way to get it into a project environment.

With [pixi](https://pixi.sh):

```bash
pixi add bioconda::fqxv
```

With `conda` (or `mamba`/`micromamba`) — Bioconda needs the `conda-forge`
channel alongside it:

```bash
conda install -c conda-forge -c bioconda fqxv
```

Pin the version for reproducibility:

```bash
pixi add "bioconda::fqxv==0.7.0"
```

Or declare it in a `pixi.toml` / `environment.yml`:

```toml
# pixi.toml
[dependencies]
fqxv = { version = "==0.7.0", channel = "bioconda" }
```

```yaml
# environment.yml
channels: [conda-forge, bioconda]
dependencies:
  - fqxv=0.7.0
```

The recipe carries a `run_exports` pin on the minor version, so while `fqxv` is
pre-1.0 an environment solved against it stays on the minor series it was built
with. That pin is about the *CLI surface*, which is still 0.x — the on-disk
format is stable at 1.0 and archives stay readable across releases regardless.

## Containers

Because fqxv is on Bioconda, [BioContainers](https://biocontainers.pro/)
automatically publishes a Docker/Singularity image for every release — no local
build required.

```bash
# Docker / Podman
docker run --rm -v "$PWD:/data" -w /data \
  quay.io/biocontainers/fqxv:0.7.0--hfa8f182_0 fqxv compress reads.fastq.gz

# Singularity / Apptainer
singularity run \
  https://depot.galaxyproject.org/singularity/fqxv:0.7.0--hfa8f182_0 fqxv --help
```

quay.io publishes no `latest` tag for Bioconda-derived images, so a tag has to
name a concrete `<version>--<build>` that really exists. The pins here track the
current Bioconda release and are refreshed weekly by a CI job
(`.github/workflows/bioconda-sync.yml`). BioContainers builds an image a day or
two behind a new release, so just after a release the tag here may still name
the previous version — see
[quay.io](https://quay.io/repository/biocontainers/fqxv?tab=tags) for every
published `<version>--<build>` tag.

The image carries only the `fqxv` CLI. Mount your working directory (the `-v`
above) so `fqxv` can read the FASTQ and write the archive back out; Singularity
bind-mounts `$PWD` by default.

## Nextflow

In [Nextflow](https://www.nextflow.io/), point a process at the image directly
or let the `conda` directive resolve it:

```groovy
process FQXV_COMPRESS {
    container 'quay.io/biocontainers/fqxv:0.7.0--hfa8f182_0'
    // or: conda 'bioconda::fqxv=0.7.0'

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.fqxv")

    script:
    """
    fqxv compress ${reads} -o ${meta.id}.fqxv --verify --threads ${task.cpus}
    """
}
```

Pin the version for reproducibility, as above. Dropping the version from the
`conda` directive (`conda 'bioconda::fqxv'`) resolves to whatever is current in
Bioconda instead, which is convenient for ad-hoc runs but makes the pipeline
non-reproducible.

Pass `--threads ${task.cpus}` so `fqxv` respects the executor's allocation
rather than its default of 16 workers. Output is deterministic regardless of
thread count, so the same input still produces a byte-identical archive when the
allocation changes. `--verify` re-decodes the fresh archive and only commits it
on a clean round-trip, which is worth the time in a pipeline that deletes its
FASTQ afterwards.

The same shape works for Snakemake (`conda:` / `container:` directives) and for
WDL/CWL (`docker:` runtime).

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
VER=v0.7.0   # the latest release tag
curl -LO https://github.com/rnabioco/fqxv/releases/download/$VER/fqxv-$VER-x86_64-unknown-linux-musl.tar.gz
tar xzf fqxv-$VER-x86_64-unknown-linux-musl.tar.gz
mv fqxv ~/.local/bin/
```

The binaries are built for each target's generic baseline; `fqxv-rans` picks its
AVX2/AVX-512 paths at runtime, so one binary runs on old and new CPUs alike.
Reach for these when you want a single static file with no environment manager
around it — Windows is binary-only, since Bioconda does not target it.

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

It ships abi3 wheels on [PyPI](https://pypi.org/project/fqxv/) and is separate
from the Bioconda CLI package: install both if you want to compress from the
shell and read archives from Python.

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
