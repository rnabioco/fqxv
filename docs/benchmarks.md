# Benchmarks

Field comparison against general-purpose and FASTQ-specific compressors on two
real Illumina RNA-seq runs spanning both quality regimes. Produced by the
harness in `bench/` (a Slurm job-array pipeline: one dataset × tool cell per
node), 4M-read subsets, 32 threads per tool.

!!! note
    Ratios and throughput vary with data type (WGS vs RNA-seq), coverage, and
    quality regime (full-range vs binned). Run the `bench/` harness on your own
    libraries before drawing conclusions.

`fqxv` has two operating points: the default (fast, block-local) and **`--max`**
(deepest context plus read reordering — the best-ratio preset). Both are
lossless and, uniquely among the tools here, **deterministic**: byte-identical
output regardless of thread count.

## NovaSeq 6000, binned quality (DRR174812) — 4M reads, 572 Mbase

| Tool | Ratio | Compress | Decompress | Peak RSS | Lossless | Deterministic |
| --- | ---: | ---: | ---: | ---: | :---: | :---: |
| SPRING | **21.9×** | 66 MB/s | 151 MB/s | 5.1 GB | yes | — |
| **`fqxv --max`** | **19.4×** | 40 MB/s | 133 MB/s | 9.8 GB | **yes** | **yes** |
| `fqxv` (default) | 9.7× | **115 MB/s** | 102 MB/s | 3.4 GB | **yes** | **yes** |
| fqz_comp | 9.6× | 101 MB/s | 69 MB/s | 75 MB | **no** | — |
| zstd -19 | 9.4× | 10 MB/s | 396 MB/s | 3.6 GB | yes | — |
| xz -9 | 8.9× | 9 MB/s | 326 MB/s | 6.3 GB | yes | — |
| gzip | 5.0× | 207 MB/s | 180 MB/s | 15 MB | yes | — |

## GAIIx, full-range quality (SRR453566) — 4M reads, 385 Mbase

| Tool | Ratio | Compress | Decompress | Lossless | Deterministic |
| --- | ---: | ---: | ---: | :---: | :---: |
| SPRING | 9.82× | 45 MB/s | 119 MB/s | yes | — |
| **`fqxv --max`** | 9.66× | **51 MB/s** | 131 MB/s | **yes** | **yes** |
| `fqxv` (default) | 7.09× | **116 MB/s** | 114 MB/s | **yes** | **yes** |
| fqz_comp | 7.05× | 91 MB/s | 70 MB/s | **no** | — |
| xz -9 | 6.07× | 6 MB/s | 265 MB/s | yes | — |
| zstd -19 | 5.87× | 11 MB/s | 502 MB/s | yes | — |
| gzip | 3.52× | 191 MB/s | 160 MB/s | yes | — |

## Reading the numbers

- **`fqxv --max` is the clear #2 on ratio**, within ~6–12% of SPRING on both
  datasets, and it beats fqz_comp, zstd -19, xz -9, and gzip decisively. On the
  full-range set it also **compresses faster than SPRING** (51 vs 45 MB/s) at
  parity ratio.
- **The default mode trades ratio for speed cleanly**: still ahead of fqz_comp /
  zstd / xz on ratio, at ~115 MB/s — an order of magnitude faster to compress
  than zstd -19 or xz -9.
- **Determinism and verified losslessness are `fqxv`'s alone here.** Every `fqxv`
  archive round-trips its exact content and is byte-identical across thread
  counts; `fqz_comp` fails the content round-trip in this harness on both sets.
- **`--max` stream split** (NovaSeq): names 0.1 MB, sequence 42.2 MB (58%),
  quality 30.9 MB (42%). The remaining ratio gap to SPRING is in the sequence
  stream's read-order handling, not the entropy coders.

## Round-trip fidelity (alignment level)

A compression ratio only counts if the reads survive intact. Beyond `fqxv`'s
internal round-trip digest, the `bench/bam_identity.sh` harness proves fidelity
*through a real analysis*: it aligns the original reads and the `fqxv`
round-tripped reads with `bwa mem` and compares the BAMs with an order-independent
multiset digest (`bench/bamcmp.rs`).

On E. coli MiSeq (SRR2627175, 2.19 M reads, GCF_000005845.2):

| Mode | BAM vs original | Detail |
| --- | --- | --- |
| lossless (default) | **byte-identical** | every field, including the coordinate-sorted file |
| `--order any` / `--order shuffle` | **byte-identical** | output order preserved on real SRA data |
| `--quality-bin bin8` | reads unmoved; QUAL only | mean \|Δ\| 1.10, max 4 Phred (72.7% of bases) |
| `--quality-bin bin4` | reads unmoved; QUAL only | mean \|Δ\| 2.79, max 10 Phred (96.7%) |
| `--quality-bin bin2` | reads unmoved; QUAL only | mean \|Δ\| 1.90, max 13 Phred (78.0%) |

- **Lossless is byte-identical at the BAM level**, not just the FASTQ: placement
  (FLAG/POS/MAPQ/CIGAR/SEQ), quality, and read names all match.
- **Lossy quality binning never moves a read.** Placement is identical for every
  bin; only the quality string changes, by a bounded amount.

!!! note "Read order is load-bearing"
    A control that shuffles the *identical* read set (same multiset, verified by
    `fqdigest`) and realigns shows `bwa mem` itself placing **~1.2%** of reads
    differently — deterministically, independent of thread count and unaffected
    by `-K`. So a reproducible BAM depends on preserving read order, which `fqxv`
    does by default and in its reorder/`--max` modes on real data. This is an
    aligner property, not a `fqxv` effect: the reads going in are byte-identical.

## Codecs

Every codec is clean-room, pure Rust — there is **no external/C compressor
dependency**. The reordered sequence path assembles a shared reference over the
clustered reads and codes each read as a position on it; the reference itself is
coded by the order-*k* `fqxv-seq` model, split into fixed blocks compressed in
parallel. Assembly, the overlap-merge, the reference coder, and the per-block
codecs all fan out with `rayon`, and every stage uses fixed, thread-independent
boundaries so output never depends on the thread count.
