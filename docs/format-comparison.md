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

| | Reference-free | Lossless seq+qual | Preserves order + names | Deterministic[^det] | Random access[^ra] | Integrity + recovery[^crc] | Long-read overlap |
| --- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`fqxv`** | ✓ | ✓ | ✓ (default / `--max`) | **✓** | ✓ | ✓ | ✓ |
| `gzip` (`.fastq.gz`) | ✓ | ✓[^blind] | ✓[^blind] | ✓ | — | CRC only | — |
| `zstd` / `xz` | ✓ | ✓[^blind] | ✓[^blind] | ✓ | —[^xz] | checksum only | — |
| `fqz_comp` | ✓ | configurable[^fqz] | ✓ | — | — | — | — |
| `SPRING` | ✓ | ✓ (reordered) | **✗ (reorders)** | — | — | — | — |
| `CoLoRd` | ✓ | ✓ (`-q org`) | ✓ | — | — | — | ✓ |
| `CRAM` | **✗ (needs reference)** | configurable | ✓ | — | ✓ | ✓ | n/a (aligned) |

[^det]: **Deterministic** = the archive is byte-identical regardless of thread
    count. `fqxv` guarantees this end to end (fixed, thread-independent block
    boundaries); the others either run single-threaded, or do not promise
    bit-reproducible output across thread counts.
[^ra]: **Random access** = seek to a read group or project a single stream (e.g.
    just the read names) without decoding the whole archive. `fqxv` carries a
    footer row-group index and per-stream CRCs, so a remote client can fetch one
    stream with a range request.
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
- **Random access and integrity in the container.** The footer row-group index +
  per-stream CRCs allow projecting a single stream with one range request; the
  per-payload CRC-32C plus `BLOCK_MAGIC` sync markers detect, localize, and recover
  from corruption rather than failing the whole file.

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
