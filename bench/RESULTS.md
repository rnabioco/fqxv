# Benchmark results — 2026-08-12

Snapshot from the unified parallel harness (`submit_parallel.sh`) on the Bodhi
`rna` partition: 9 datasets × platform-appropriate tools = 150 cells, each fanned
out one-tool-per-node, all COMPLETED. Covers short-read Illumina (MiSeq,
NovaSeq6000, GAIIx) and MGI/BGISEQ, plus long-read ONT (MinION) and PacBio HiFi
on **both** Sequel II and Revio. Ratios are deterministic (thread-count
independent); `rt=yes` means the round-trip content digest matched (lossless, or
lossy-expected for the binned points). Reproduce with `bash submit_parallel.sh`.

Two builds, stated per row set. The matrix ran at `4d6b434` (pre-mode-5); every
**fqxv long-read cell was then re-measured at `102fa86`** — the v0.7.0 codebase
with chunked quality (#279) — and those are the fqxv numbers below. Short-read
fqxv cells and all field-tool cells stand from the `4d6b434` run: #279 leaves
short-read output **byte-identical** (verified), and tool cells do not depend on
the fqxv build. Every short-read ratio reproduced **byte-for-byte** from the
2026-07-22 tables; of the long-read datasets only `hifi_revio_wgs` drifted at
the third digit (uniformly, every tool including gzip — input restaging, not
codec).

What moved this cycle is a deliberate trade: **default-level long-read archives
now give up a sliver of ratio to become parallel-decodable.** The Nanopore block
budget dropped 256 → 64 MiB (#275; blocks are the unit of decode parallelism),
and long-read quality — 92–98% of decode — is now coded in within-block chunks
(#279, context mode 5/6): together **ONT default 3.01× → 2.96×** and **HiFi
4.77× → 4.76×**, for full decode **~1.9× faster at 16+ threads** (the matrix
cells measured ONT decode 26.3 s → 13.5 s at 16 threads) and HiFi compress
**17% faster**. `-l9`/`--max` pins the serial maximal-block layout — those
archives are **byte-identical** to the previous release, so every headline
best-ratio number (ONT 3.06×, HiFi 4.77×) is unchanged. See
`docs/decode-scaling.md` for the decode story.

## Short-read, lossless

`fqxv --order shuffle` (rename-and-renumber, lossless) beats SPRING on the three
Illumina datasets and on ultra-short miRNA reads, but **loses to SPRING on the
hard MGI case**. `fqxv-max` is the order-preserving lossless point.

| dataset          | fqxv | fqxv-max | fqxv-shuffle | spring | zstd19 | xz9 | gzip |
|------------------|-----:|---------:|-------------:|-------:|-------:|----:|-----:|
| rnaseq_novaseq   | 9.94 |    23.80 |    **29.16** |  25.27 |   9.59 | 9.01| 5.05 |
| mgi_mirna_22bp   |10.14 |    10.84 |    **15.14** |  12.22 |   6.38 | 7.57| 4.70 |
| rnaseq_fullrange | 7.60 |    10.20 |    **11.62** |  10.01 |   5.84 | 5.98| 3.71 |
| ecoli_miseq      | 4.92 |     7.41 |     **7.35** |   7.32 |   5.24 | 5.13| 2.94 |
| mgi_bgiseq_hard  | 3.80 |     3.80 |         3.89 |**4.05**|   3.14 | 3.26| 2.57 |

Every value is unchanged from 2026-07-22 (this cycle's codec work chunks
long-read quality only; short-read archives are byte-identical). SPRING still
wins the hardest short-read case (BGISEQ-500, 137 bp, 5.32 bits/base) at 4.05.
On **22 bp miRNA reads** fqxv wins by a wide margin (15.14 vs SPRING's 12.22):
at that length the read names dominate the archive and the tokenizer is the
whole game.

## Long-read, lossless

The best-ratio points are unchanged — **ONT still edges ahead of CoLoRd**
(3.06× vs 3.05× at `-l9`/`--max`) and **HiFi still leads on every stream**
(4.77×) — while the *default* level now trades ~0.2–1.5% of archive for
parallel decode (64 MiB Nanopore blocks + chunked quality). All rows round-trip
losslessly; the fqxv columns are the v0.7.0 (`102fa86`) build.

| dataset             |  fqxv | fqxv9   |    colord | zstd19 |   xz9 | gzip |
|---------------------|------:|--------:|----------:|-------:|------:|-----:|
| hifi_revio_amplicon | 21.79 |**21.80**|     19.88 |  15.12 | 14.45 | 8.96 |
| hifi_revio_wgs      | 16.96 |   16.96 | **18.77** |  12.97 | 12.75 | 9.25 |
| ecoli_hifi          |  4.76 |**4.77** |      4.44 |   3.83 |  3.85 | 2.27 |
| ecoli_ont           |  2.96 |**3.06** |      3.05 |   2.38 |  2.49 | 1.94 |

The per-stream split (fqxv best lossless point — the serial `-l9`/`--max`
layout, byte-identical to the previous release — from `fqxv info`) shows where
the levers landed:

| dataset             | seq (b/base) |  seq (bytes) | quality (bytes) | quality share |
|---------------------|-------------:|-------------:|----------------:|--------------:|
| ecoli_hifi          |        0.064 |   12,286,488 |     635,629,989 |       **98%** |
| hifi_revio_amplicon |        0.083 |    6,276,266 |      49,183,082 |           84% |
| hifi_revio_wgs      |    **0.683** |  115,329,481 |      44,786,803 |           28% |
| ecoli_ont (`-l9`)   |    **0.887** |   33,404,890 |     163,660,203 |           83% |

**The default level now buys parallel decode with those streams.** Chunked
quality (#279, mode 5/6) recodes each block's quality as a serially-coded warmup
prefix plus K−1 whole-read chunks that encode and decode in parallel: ONT
quality 164.5M → 165.2M (**+0.36%** of archive, matching the #277 sweep's
prediction to the third digit) and HiFi 635.6M → 636.9M (**+0.19%**), sequence
and names byte-identical. On the pre-binned Revio runs the cost is noise
(amplicon +0.07%, WGS +0.004%). What it buys: long-read full decode **~1.9×
faster at 16–32 threads** (ONT 17.7 → 9.1 s, HiFi 66.5 → 33.8 s on a dedicated
32-CPU node; the shared-node matrix cells corroborate at 16 threads with ONT
26.3 → 13.5 s), and HiFi *compress* **17% faster** — the chunks parallelize the
encoder too. A v0.6.x reader refuses a mode-5/6 archive loudly
(`unsupported quality context mode 5; this archive needs a newer fqxv`);
`--max` keeps the serial smallest-archive layout, byte-identical across the
releases.

**`ecoli_ont` default 3.01× → 2.96×** decomposes as: 64 MiB block budget (#275,
−1.2%: the benchmark file cuts 5 blocks instead of 2, sequence loses some
cross-read reach) plus chunked quality (#279, −0.36%). Decode of that same file
went from ~11 MB/s flat at every thread count to scaling through 32 threads —
see `docs/decode-scaling.md`. The `-l9`/`--max` point keeps maximal blocks and
serial quality: 3.06×, byte-identical, still ahead of CoLoRd's 3.05×.

**`ecoli_hifi` (Sequel II, ~300×) still leads CoLoRd on every stream** at
`--max` — sequence 12.3M vs 13.4M, quality 635.6M vs 684.3M, total 649.5M vs
697.7M (4.77× vs 4.44×) — and the chunked default gives up only +0.19% (4.76×).
This is why one HiFi dataset was never enough: `ecoli_hifi` is 98% quality by
bytes (300× coverage of a 4.6 Mb genome collapses the sequence stream), while
Revio WGS at ordinary coverage is 72% sequence. The two datasets exercise
opposite regimes.

**`hifi_revio_wgs`** stays a strong second behind CoLoRd (16.96× vs 18.77×,
ahead of zstd19 12.97× and xz9 12.75×) on the raw-LZMA sequence path; the
chunked-quality cost is invisible (+0.004%) because its quality stream is
pre-binned and only 28% of the archive.

## Long-read compress speed

Three waves of long-read work land here. First, the **byte-identical** speedups
of the previous cycles — gating off the always-discarded overlap-consensus
candidate on Nanopore (#223), de-packing the banded-DP traceback (#222), and
skipping the redundant shared-reference assembly on Nanopore (#211) — made
default ONT compress ~3× faster at identical output. Then the
anchor-restricted tiler coder (#231) made ONT another **1.54× faster at 16
threads** *and* smaller. This cycle, **chunked quality (#279) makes default
HiFi compress 17% faster** (317 s → 263 s at 16 threads on a dedicated node —
the quality encoder was the bottleneck, and the chunks encode in parallel);
ONT compress is tiler-bound, so its encode cost is flat (+1.8%). The deepest
sequence lever (best-of-4 tiling references, band 768) stays gated to
`-l9`/`--max`, which spends ~2.5× the default's compress time for the
best-ratio 3.06×.

## Lossy quality (fqxv binning)

Quality binning is the big lever, especially on long reads. `reorder-bin*` stacks
read-reordering on top of binning.

| dataset          | bin8 | bin4 | bin2 | reorder-bin2 |
|------------------|-----:|-----:|-----:|-------------:|
| rnaseq_novaseq   | 9.94 | 9.94 | 11.37|        33.98 |
| rnaseq_fullrange |10.84 |14.49 | 15.37|        31.95 |
| ecoli_ont (binont)| 7.05|   —  |   —  |            — |

The short-read rows are unchanged. `ecoli_ont --quality-bin ont` reads 7.17 →
**7.05** this cycle for the same reason the lossless default moved: the binned
ONT archive is sequence-dominated, so the 64 MiB block budget (#275) costs it
~1.7%, and chunking the (already small) binned quality stream adds +0.03%. In
exchange the lossy archive decodes with the same block/chunk parallelism as the
lossless one.

## fqxv archive vs native NCBI `.sra` (lossless `max` regime)

`fqxv/.sra` < 1 means the lossless `fqxv --max` archive is smaller than the native
`.sra` the run ships in; it wins on every platform (`sra_compare.sh`, both mates,
`.sra` sizes from `sracha info`).

| accession  | platform    | fqxv / .sra |
|------------|-------------|------------:|
| DRR174812  | NovaSeq6000 |       0.331 |
| SRR453566  | GAIIx       |       0.509 |
| SRR2627175 | MiSeq       |       0.538 |
| DRR205413  | ONT-MinION  |       0.721 |

fqxv is ~2× smaller than the `.sra` on average (geomean 0.51). NovaSeq wins most
(its quality is pre-binned, so the lossless point is already compact); ONT wins
least, since the long-read sequence stream dominates that archive. Measured
2026-07-23 with `--max` archives, which are byte-identical on this release, so
the table stands. A grouped size chart of this table is on the
[Benchmarks](../docs/benchmarks.md) page.
