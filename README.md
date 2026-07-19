# 🗜️ fqxv

A Rust toolkit for archiving FASTQ, built as a workspace of
one-crate-per-algorithm codecs plus a container format and CLI.

![fqxv compress and decompress demo](docs/images/readme.gif)

> [!WARNING]
> **Not production-ready — the format is still stabilizing.** `fqxv` is in
> **early development** (v0.2.0). The library and CLI work end-to-end and are
> [benchmarked against the field](docs/benchmarks.md), but **do not use `fqxv` as
> the only copy of data you care about.**
>
> - **The on-disk `.fqxv` format (version 1.0) is not frozen and may change
>   incompatibly before a stability guarantee.** A reader accepts any archive with
>   its own major version (newer minors are tolerated) but refuses a different
>   major, so a future major bump can make today's archives unreadable. Pin an
>   exact version if you need archives to survive an upgrade — and prefer
>   re-compressing from the original FASTQ.
> - **Bugs are still being found.** Correctness work is ongoing, including
>   defects that produced archives which could not be decompressed. Archives are
>   checksummed and `compress --verify` reads one back to confirm, but **keep
>   your original FASTQ** until the format is frozen.
>
> A frozen format and a stability guarantee will be announced before `fqxv` is
> recommended for archival use.

## Why

FASTQ splits into three streams that compress very differently, and the
wins are additive:

| Stream | Share of a lossless archive | What moves it |
| --- | --- | --- |
| Quality scores | ~50–70% | context-model entropy coder (fqzcomp-class) |
| Sequence (bases) | most of the rest | read reordering / assembly (PgRC2/SPRING-class) |
| Read names | small | positional tokenizer |

`fqxv` handles each with a clean-room codec implemented from the [CRAM 3.1 codecs
spec](https://samtools.github.io/hts-specs/CRAMcodecs.pdf) and the source papers —
no mature FASTQ-domain compressor exists as a Rust crate today. On 4M-read RNA-seq
subsets it is the **smallest lossless compressor of the field**, beating SPRING,
`fqz_comp`, `zstd -19`, `xz -9`, and `gzip` (see [benchmarks](docs/benchmarks.md)):

| NovaSeq (binned), 4M reads | ratio | lossless |
| --- | ---: | :---: |
| **`fqxv --order shuffle`** | **23.9×** | seq+qual (renumbered) |
| SPRING | 21.9× | reordered |
| **`fqxv --max`** | 20.2× | **yes (order-preserving)** |
| fqz_comp | 9.6× | no (fails round-trip) |
| zstd -19 / xz -9 | 9.4× / 8.9× | yes |

Every archive is **deterministic** (byte-identical regardless of thread count) and
**verified lossless** on decode. Pure Rust, no external/C compressor.

## Install

Until `fqxv` lands on [bioconda](https://bioconda.github.io/), install from source
with Cargo (Rust 1.95+):

```bash
cargo install --git https://github.com/rnabioco/fqxv fqxv-cli
```

This builds the `fqxv` binary.

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

Long reads (ONT/PacBio) compress correctly today, and the quality stream is
already at parity with CoLoRd; the sequence stream is not yet competitive at high
coverage. See [long-read support](docs/design/longread.md) for the measurements.

## Acknowledgments

`fqxv` stands on a large body of prior work. Everything here is a clean-room
implementation from public specifications and papers — no third-party source is
vendored — but these projects and their authors made it possible, and we
cross-checked against several of them for correctness:

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
