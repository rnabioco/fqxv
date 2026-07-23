# Format comparison

Where `fqxv` sits among the formats people use to store FASTQ. This page compares
**capabilities and guarantees**; for compression ratios and throughput on real
data, see [Benchmarks](benchmarks.md).

FASTQ archiving splits into a few families:

- **General-purpose byte compressors** — `gzip` (`.fastq.gz`), `zstd`, `xz`. They
  compress the FASTQ *bytes* without parsing them. Everything round-trips because
  nothing is interpreted, but they gain nothing from FASTQ structure (the three
  streams, the small quality alphabet, cross-read redundancy).
- **FASTQ-specialized compressors** — `fqz_comp`, `SPRING`. They split names /
  sequence / quality and model each; SPRING additionally reorders reads to expose
  cross-read redundancy.
- **Long-read specialized** — `CoLoRd`. Codes each long read against a similar
  earlier read (overlap + edit script), the dominant lever at high coverage.
- **Reference / repository formats** — `CRAM` (needs the alignment and a reference
  genome) and NCBI's native `.sra` (the columnar archive SRA distributes). Different
  category: not something you run over a raw FASTQ file.

`fqxv` is a FASTQ-specialized archiver that also carries the long-read overlap
lever, is reference-free, and adds container-level guarantees (determinism,
verified round-trips, random access, integrity) that the others do not make.

## Capability matrix

| | Reference-free | Lossless seq+qual | Preserves order + names | Deterministic[^det] | Forward streaming[^stream] | Seekable random access[^ra] | Integrity + recovery[^crc] | Long-read overlap |
| --- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`fqxv`** | ✓ | ✓ | ✓ (default / `--max`) | **✓** | ✓ | ✓ | ✓ | ✓ |
| `gzip` (`.fastq.gz`) | ✓ | ✓[^blind] | ✓[^blind] | ✓ | ✓ | — | CRC only | — |
| `zstd` / `xz` | ✓ | ✓[^blind] | ✓[^blind] | ✓ | ✓ | —[^xz] | checksum only | — |
| `fqz_comp` | ✓ | configurable[^fqz] | ✓ | — | ✓ | — | — | — |
| `SPRING` | ✓ | ✓ (reordered) | **✗ (reorders)** | — | ✗[^spring] | — | — | — |
| `CoLoRd` | ✓ | ✓ (`-q org`) | ✓ | — | ✗ | — | — | ✓ |
| `CRAM` | **✗ (needs reference)** | configurable | ✓ | — | ✗[^seekonly] | ✓ | ✓ | n/a (aligned) |
| `SRA` / VDB | ✓ | **✗ (`.sralite`; not byte-exact)**[^sra] | ✓[^sra] | — | ✗[^seekonly] | ✓ | ✓ | n/a |

[^det]: **Deterministic** = the archive is byte-identical regardless of thread
    count. `fqxv` guarantees this end to end (fixed, thread-independent block
    boundaries); the others either run single-threaded, or do not promise
    bit-reproducible output across thread counts.
[^stream]: **Forward streaming** = decode as a single forward pass over a
    *non-seekable* pipe (e.g. `aws s3 cp s3://bkt/reads.fqxv - | fqxv decompress -`),
    with no local copy and constant memory — the first records come out before the
    last byte arrives. General byte compressors stream this way too, but cannot
    project a FASTQ stream. `SPRING`/`CoLoRd` need scratch disk and multiple passes;
    `CRAM`/`SRA` are random-access database formats that require a *seekable* source
    (a full download or many Range seeks), so they cannot decode off a forward pipe.
[^ra]: **Seekable random access** = seek to a read group or project a single stream
    (e.g. just the read names) with a range request, given a *seekable* source, without
    decoding the whole archive. `fqxv` carries a footer row-group index and per-stream
    CRCs; `CRAM` (`.crai`) and `SRA`/VDB (columnar) also index for random access.
[^seekonly]: Random-access only: needs a seekable backing store (a local file or HTTP
    Range) plus an index, and cannot decode incrementally off a forward pipe. SRA's
    recommended workflow (`prefetch`) downloads the whole archive before decoding.
[^sra]: NCBI's repository format (VDB), consumed through the SRA Toolkit. Distributed
    in a full-quality **Normalized** form (`.sra`) and a **Lite** form (`.sralite`) that
    collapses quality to a single value per read; neither guarantees byte-exact FASTQ
    round-trip. A different category — the archive NCBI distributes, not a tool you run
    on your own FASTQ.
[^spring]: `SPRING` compresses in a multi-pass pipeline that needs a scratch directory
    (~10–30% of input, more when lossy), so it cannot stream from a pipe.
[^crc]: **Integrity + recovery** = a checksum on every coded payload *and* the
    ability to resynchronize to intact blocks when the index or a length prefix is
    corrupt. `fqxv` carries a CRC-32C per payload plus `BLOCK_MAGIC` sync markers
    for `decompress_recover`; general compressors carry a whole-stream checksum
    that detects corruption but cannot localize or recover past it.
[^blind]: General byte compressors preserve order and names because they never
    parse the FASTQ — they compress the raw bytes. That is also why they cannot
    exploit the quality alphabet or cross-read redundancy.
[^xz]: `xz` can carry a block index and `zstd` supports seekable framing, but
    neither gives stream-level (names / sequence / quality) projection of a FASTQ.
[^fqz]: `fqz_comp` can be run losslessly, but its quality model is lossy by
    default and it fails the order-independent content round-trip in the `bench/`
    harness on both Illumina sets — see [Benchmarks](benchmarks.md).

## What is unique to `fqxv`

Three guarantees in the matrix are `fqxv`'s alone among the FASTQ tools here:

- **Deterministic output.** Blocks are the unit of parallelism and every stage
  uses fixed, thread-independent boundaries, so the same input always produces the
  same bytes — whether you run on 1 core or 88. Reproducible archives are a
  first-class property, not an accident of a single-threaded run.
- **Verified losslessness.** Every block prepends xxh3-64 digests of its *decoded*
  content (names, sequence, post-binning quality), checked after decode, and
  `compress --verify` round-trips the whole archive before it is written. A codec
  bug that produced CRC-valid but wrong bytes is caught at runtime.
- **Remote access — both flavors — and integrity in the container.** `fqxv` is the
  only tool here that both *streams forward* off a non-seekable pipe (decode as bytes
  arrive, `… | fqxv decompress -`, constant memory) and offers *seekable* per-stream
  random access (the footer row-group index + per-stream CRCs project a single stream
  with one range request). `CRAM` and `SRA` provide the seekable half only, and only
  from a seekable source; the rest offer neither. The per-payload CRC-32C plus
  `BLOCK_MAGIC` sync markers additionally detect, localize, and recover from corruption
  rather than failing the whole file.

It also spans **both** short- and long-read levers in one format: read reordering
for short-read cross-read redundancy (SPRING-class) and the `fqxv-lroverlap`
overlap codec for long reads (CoLoRd-class), auto-selected per block, with the
smaller of overlap-vs-order-*k* kept so a block never regresses.

## Operating points

`fqxv` exposes the ratio/fidelity trade as three presets rather than a single
mode (see [Benchmarks](benchmarks.md) for numbers):

| Preset | Ratio | Read order | Names | Guarantee |
| --- | --- | --- | --- | --- |
| `fqxv` (default) | fast, modest | preserved | preserved | fully lossless |
| `--max` | best fully-lossless | preserved | preserved | fully lossless (stores a permutation) |
| `--order shuffle` | smallest | renumbered | discarded | seq+qual multiset (the trade SPRING makes) |

Lossy quality (`--quality-bin bin8|bin4|bin2` for short reads, `ont` / `hifi` for
long reads) is a separate tier layered on any preset — compare those against each
other, not against the lossless rows.

## Choosing a format

- **Reproducible, fully lossless archive you can seek into** → `fqxv` (default for
  speed, `--max` for the best fully-lossless ratio). The only option here that is
  deterministic *and* offers per-stream random access.
- **Stream or project straight from object storage (S3/HTTP)** → `fqxv`. It decodes
  off a pipe (`aws s3 cp … - | fqxv decompress -`) and can fetch a single stream with a
  Range request. `CRAM` and `SRA` support remote random access too, but only from a
  seekable source (or a full local copy) and in a different category — reference-based
  alignments or NCBI's columnar repository.
- **Smallest lossless short-read archive** → `fqxv --order shuffle` or `SPRING`;
  `fqxv --order shuffle` is smaller on both benchmark sets under the same
  reorder+renumber rules.
- **Long reads** → `fqxv` (cross-read sequence codecs) or `CoLoRd`. `fqxv`'s
  quality *leads* CoLoRd on both platforms (a sequence-conditioned, context-mixed
  coder), and the sequence stream now closes the gap too: a whole-file overlap
  reference on high-coverage HiFi, raw-LZMA on ordinary-coverage HiFi, and
  best-of-N tiling with anchor-restricted coding on Nanopore bring fqxv **past
  CoLoRd on both ONT and HiFi Sequel II**. See [Long-read support](design/longread.md).
- **Aligned reads with a reference on hand** → `CRAM` is purpose-built for that and
  reference-based; `fqxv` is for the raw, reference-free FASTQ.
- **Maximum portability, no tooling** → `.fastq.gz`. Universally readable, but the
  largest and structure-blind.

!!! note
    Capabilities are stable properties; ratios are not. Which format is *smallest*
    depends on platform, coverage, and quality regime — run the `bench/` harness on
    your own libraries before deciding.
