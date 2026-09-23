# 🗜️ fqxv

[![Bioconda](https://img.shields.io/conda/vn/bioconda/fqxv?label=bioconda)](https://anaconda.org/bioconda/fqxv)
[![Bioconda downloads](https://img.shields.io/conda/dn/bioconda/fqxv?label=downloads)](https://anaconda.org/bioconda/fqxv)
[![PyPI](https://img.shields.io/pypi/v/fqxv)](https://pypi.org/project/fqxv/)
[![License](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](#license)

A Rust toolkit for archiving FASTQ, built as a workspace of
one-crate-per-algorithm codecs plus a container format and CLI.
[Documentation](https://rnabioco.github.io/fqxv/).

![fqxv compress and decompress demo](docs/images/readme.gif)

**The on-disk format is stable at 1.0** — a number independent of the crate
version (currently 0.7.x); the format and the software are versioned separately.
Archives written today stay readable by later releases: a reader accepts its own
format major version and tolerates newer minors, additions it can safely ignore
(a skippable header record) it skips, and anything it cannot — a required feature
bit or a codec it does not know — it refuses with an "upgrade fqxv" error rather
than misreading. Compatibility fails loudly, never silently. Every archive is
**deterministic** (byte-identical regardless of thread count), checksummed, and
**verified lossless** on decode.

## Results

Best lossless ratio per platform, against the strongest alternative for that data
— whole-file runs, from [`bench/RESULTS.md`](bench/RESULTS.md)
([full benchmarks](docs/benchmarks.md)):

| Platform | fqxv | best alternative | gzip |
| --- | ---: | ---: | ---: |
| Illumina NovaSeq (binned quality) | **29.2×** | SPRING 25.3× | 5.0× |
| Illumina GAIIx (full-range quality) | **11.6×** | SPRING 10.0× | 3.7× |
| Illumina MiSeq (*E. coli* WGS) | **7.4×** | SPRING 7.3× | 2.9× |
| PacBio HiFi (*E. coli*, ~300×) | **4.8×** | CoLoRd 4.4× | 2.3× |
| Oxford Nanopore (MinION) | **3.06×** | CoLoRd 3.05× | 1.9× |

The Illumina rows use `--order shuffle`, which reorders reads and renumbers names —
the same trade SPRING's own mode makes, so the comparison is like-for-like; the
Nanopore row uses `--max`. Fully order-preserving and default-effort numbers are in
the [benchmarks](docs/benchmarks.md); those tables measure 4M-read subsets, so they
report lower ratios than the whole-file runs above. The full per-dataset matrix is
in [`bench/RESULTS.md`](bench/RESULTS.md), including the hardest MGI/BGISEQ case,
where SPRING is still smaller.

fqxv now **edges ahead of CoLoRd on Nanopore** (3.06× vs 3.05×): its ONT *quality*
stream is already the smaller of the two, and a best-of-N tiling sequence codec with
anchor-restricted coding (engaged at `-l9`/`--max`) closes the cross-read sequence
gap that used to trail. On HiFi at ~300× fqxv leads on every stream; on a modern
Revio WGS run at ordinary coverage, a new raw-LZMA sequence path nearly doubled
fqxv's ratio (9.7× → 17×), landing just behind CoLoRd — see the
[benchmarks](docs/benchmarks.md).

Pure Rust, no external or C compressor.

Beyond ratio, fqxv is the only reference-free FASTQ archiver that both streams
forward off a pipe (`… | fqxv decompress -`) and supports seekable per-stream
random access — see the [format comparison](docs/format-comparison.md) for how it
stacks up against gzip/zstd, SPRING, CoLoRd, CRAM, and SRA on determinism,
losslessness, integrity, and remote access.

## Install

Install via [Bioconda](https://bioconda.github.io/):

```bash
pixi add bioconda::fqxv
# or: conda install -c conda-forge -c bioconda fqxv
```

Or download a prebuilt binary — every release attaches a static `fqxv`
(Linux x86-64/arm64, macOS Intel/Apple silicon, Windows x86-64) plus a
`SHA256SUMS.txt` to its
[GitHub Release](https://github.com/rnabioco/fqxv/releases):

```bash
VER=v0.7.0   # the latest release tag
curl -LO https://github.com/rnabioco/fqxv/releases/download/$VER/fqxv-$VER-x86_64-unknown-linux-musl.tar.gz
tar xzf fqxv-$VER-x86_64-unknown-linux-musl.tar.gz && mv fqxv ~/.local/bin/
```

Or build from source with Cargo (Rust 1.95+):

```bash
cargo install --git https://github.com/rnabioco/fqxv fqxv-cli
```

### Containers

Because fqxv is on Bioconda, [BioContainers](https://biocontainers.pro/)
automatically publishes a Docker/Singularity image for every release — no
local build required.

```bash
# Docker / Podman
docker run --rm quay.io/biocontainers/fqxv:0.7.0--hfa8f182_0 fqxv --help

# Singularity / Apptainer
singularity run \
  https://depot.galaxyproject.org/singularity/fqxv:0.7.0--hfa8f182_0 fqxv --help
```

quay.io publishes no `latest` tag for Bioconda-derived images, so a tag has to
name a concrete `<version>--<build>`; the pins here track the current Bioconda
release and are refreshed weekly by a CI job. BioContainers builds an image a
day or two behind a new release, so just after a release the tag here may still
name the previous version — see
[quay.io](https://quay.io/repository/biocontainers/fqxv?tab=tags) for every
published tag.

### Python

Any of the above gives you the `fqxv` binary. A read-only Python package reads
`.fqxv` archives directly (compression stays in the CLI):

```bash
uv pip install fqxv
```

Using fqxv in a workflow manager? See
[Nextflow and workflow integration](https://rnabioco.github.io/fqxv/getting-started/installation/#nextflow)
in the docs. See
[Installation](https://rnabioco.github.io/fqxv/getting-started/installation/)
for the full asset list and the crate-level dependencies.

## Usage

The [demo above](docs/images/readme.gif) shows the basic round trip:

```bash
fqxv compress reads.fastq.gz    # gzip input auto-detected; writes reads.fqxv
fqxv decompress reads.fqxv -o reads.fastq
```

Paired-end and single-cell inputs interleave into one archive
(`fqxv compress R1.fq.gz R2.fq.gz -o sample.fqxv`). Lossless by default;
`--quality-bin {bin8,bin4,bin2,ont,hifi}` opts into lossy quality (`ont`/`hifi`
are the long-read tables) and `--max` chases the best ratio. Add `--verify` to
re-decode the new archive and confirm it round-trips before you trust (or delete)
the source. Run `fqxv --help` for the full option set.

Long reads (ONT/PacBio) get a cross-read overlap sequence codec and long-read
quality models automatically, based on the detected platform — see [long-read
support](docs/design/longread.md).

## Acknowledgments

`fqxv` stands on a large body of prior work. Everything here is a clean-room
implementation from public specifications and papers — no third-party source is
vendored, with one deliberate exception for the cryptographic primitives behind
optional archive encryption (see [`THIRD-PARTY-NOTICES.md`](THIRD-PARTY-NOTICES.md))
— but these projects and their authors made it possible, and we cross-checked
against several of them for correctness:

- **htscodecs** ([samtools/htscodecs](https://github.com/samtools/htscodecs),
  James Bonfield / Genome Research Ltd) and the [CRAM 3.1 codecs
  spec](https://samtools.github.io/hts-specs/CRAMcodecs.pdf) — the reference for
  our rANS Nx16 coder, fqzcomp quality model, and name tokenizer.
- **fqzcomp** (James Bonfield) — the quality-score context model our
  `fqxv-fqzcomp` codec is modeled on.
- **noodles** ([zaeleus/noodles](https://github.com/zaeleus/noodles),
  Michael Macias) — Rust CRAM codec implementation we cross-checked test vectors
  against.
- **rANS / ryg_rans** — Jarek Duda's asymmetric numeral systems and Fabien
  Giesen's `ryg_rans` (public domain / CC0), plus Eugene Shelwien's range-coder
  design, underpin our entropy coders.
- **SPRING** (Chandak et al., *Bioinformatics* 2019) and **PgRC2** (Kowalski &
  Grabowski, *Bioinformatics* 2025) — the algorithmic references for the
  read-reordering engine.
- **CoLoRd** (Kokot et al., *Nature Methods* 2022) and **minimap2** / **miniasm**
  (Heng Li, *Bioinformatics* 2018 / 2016) — the references for the long-read
  quality-binning tables and the long-read overlap work.

See [`THIRD-PARTY-NOTICES.md`](THIRD-PARTY-NOTICES.md) for licenses and full
attribution.

## License

Licensed under either of [Apache License, Version 2.0](LICENSE-APACHE) or
[MIT license](LICENSE-MIT) at your option.
