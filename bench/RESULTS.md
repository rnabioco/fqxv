# Benchmark results — 2026-07-22

Snapshot from the unified parallel harness (`submit_parallel.sh`) on the Bodhi
`rna` partition: 9 datasets × platform-appropriate tools = 150 cells, each fanned
out one-tool-per-node, all COMPLETED. Covers short-read Illumina (MiSeq,
NovaSeq6000, GAIIx) and MGI/BGISEQ, plus long-read ONT (MinION) and PacBio HiFi
on **both** Sequel II and Revio. Ratios are deterministic (thread-count
independent); `rt=yes` means the round-trip content digest matched (lossless, or
lossy-expected for the binned points). Reproduce with `bash submit_parallel.sh`.

Re-run on `2468cb2`. Every short-read ratio, and the ONT and Revio-WGS long-read
ratios, reproduced **byte-for-byte** from the prior tables. What moved is **HiFi**:
the anchor-restricted HiFi consensus coder (#234) shaved the sequence stream and
the per-block quality quantizer (#235, trialled and kept only when smaller) shaved
the quality stream, lifting `ecoli_hifi` **4.73× → 4.77×** and `hifi_revio_amplicon`
21.79× → 21.80×. The same quality quantizer shrank the binned ONT stream, taking
`ecoli_ont --quality-bin ont` **6.68× → 7.17×**.

Two long-read sequence codecs landed since the 2026-07-20 run and move the
numbers below: a **raw-LZMA** sequence method for ordinary-coverage long reads
(#197/#205–208 — the PacBio Revio WGS lever) and a **multi-reference tiling**
codec for Nanopore (#212/#213 — the ONT lever, engaged at `-l9`/`--max`). Both
are coded per block and kept only when they beat the alternatives, so no
short-read or high-coverage result regresses. The robustness corpus (48 accessions
across ONT, HiFi, Illumina, and MGI, `corpus.sh`) round-trips losslessly and builds
deterministically on this build — no failures.

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

The short-read numbers are unchanged from 2026-07-20 (this cycle's codec work was
long-read only), with one improvement on **`mgi_bgiseq_hard`**: `fqxv-max` used to
be *worse* than plain `fqxv` (3.78 vs 3.80), because the level-9 reorder + hashed
tier cost more than it recovered on data with little cross-read redundancy. The
never-worse gate (#196/#202) now floors `--max`/`-l9` at the default cost, so both
sit at 3.80 and `--max` can no longer lose to the default. SPRING still wins this
hardest short-read case (BGISEQ-500, 137 bp, 5.32 bits/base) at 4.05. On **22 bp
miRNA reads** fqxv wins by a wide margin (15.14 vs SPRING's 12.22): at that length
the read names dominate the archive and the tokenizer is the whole game.

## Long-read, lossless

**HiFi Revio WGS — last cycle's single worst result — nearly doubled** (9.72× →
**16.95×**) via the raw-LZMA sequence path, and **ONT now edges ahead of CoLoRd**
(3.06× vs 3.05×) after multi-reference tiling plus the anchor-restricted tiler
coder (#231, which aligns only inter-anchor gaps). All rows round-trip losslessly.

| dataset             |     fqxv |  fqxv9 |    colord | zstd19 |   xz9 | gzip |
|---------------------|---------:|-------:|----------:|-------:|------:|-----:|
| hifi_revio_amplicon |**21.80** |  21.80 |     19.88 |  15.12 | 14.45 | 8.96 |
| hifi_revio_wgs      |    16.95 |  16.95 | **18.76** |  13.02 | 12.75 | 9.28 |
| ecoli_hifi          | **4.77** |   4.77 |      4.44 |   3.83 |  3.85 | 2.27 |
| ecoli_ont           |     3.01 |**3.06**|      3.05 |   2.38 |  2.49 | 1.94 |

The per-stream split (fqxv best lossless point, from `fqxv info`) shows where the
two levers landed:

| dataset             | seq (b/base) |  seq (bytes) | quality (bytes) | quality share |
|---------------------|-------------:|-------------:|----------------:|--------------:|
| ecoli_hifi          |        0.064 |   12,286,488 |     635,629,989 |       **98%** |
| hifi_revio_amplicon |        0.083 |    6,276,266 |      49,183,082 |           84% |
| hifi_revio_wgs      |    **0.683** |  115,329,481 |      44,786,803 |           28% |
| ecoli_ont (`-l9`)   |    **0.887** |   33,404,890 |     163,660,203 |           83% |

**`hifi_revio_wgs` was the worst result in the suite** — 1.391 b/base of sequence
(234.7M), an archive that barely cleared gzip and lost to zstd/xz. A real genome
at ordinary coverage carries **exact** cross-read matches that neither the
within-read order-k model nor a per-block voted consensus can exploit; the raw
large-window LZMA sequence method catches them and halves the stream to **0.683
b/base (115.3M)**. The archive goes from second-worst-lossless to a strong second
behind CoLoRd, now ahead of both zstd19 (13.02) and xz9 (12.75).

**`ecoli_ont` now beats CoLoRd.** The multi-reference tiling codec codes each read
against earlier *raw* reads (best-of-4 references at `-l9`/`--max`), and the
anchor-restricted coder (#231) then aligns only the inter-anchor gaps rather than
re-deriving the whole read — together cutting `-l9` ONT sequence to **0.887 b/base
(33.4M)** versus CoLoRd's 31.4M (now ~6% behind on sequence alone). fqxv's quality
coder makes up the rest: 163.7M versus CoLoRd's 166.5M (a 2.8M credit), so the ONT
total lands at **3.063×, ahead of CoLoRd's 3.052×** (188.0M vs 188.7M) — the whole
archive is now smaller. The anchor coder also lifts default `fqxv` (2.92× → 3.01×,
−2.9% archive) and makes ONT compress **1.54× faster at 16 threads**; the tiling
depth is effort-gated so the deepest references are spent only at `-l9`/`--max`.

**`ecoli_hifi` (Sequel II, ~300×) still leads CoLoRd on every stream** — sequence
12.3M vs 13.4M, quality 635.6M vs 684.3M, total 649.5M vs 697.7M — at 4.77×. This
is why one HiFi dataset was never enough: `ecoli_hifi` is 98% quality by bytes
(300× coverage of a 4.6 Mb genome collapses the sequence stream), while Revio WGS
at ordinary coverage is 72% sequence. The two datasets exercise opposite regimes,
and the sequence lever only shows up on the second.

## Long-read compress speed

Two waves of long-read work landed here. First, a set of **byte-identical**
speedups — gating off the always-discarded overlap-consensus candidate on Nanopore
(#223), de-packing the banded-DP traceback (#222), and skipping the redundant
shared-reference assembly on Nanopore (#211): on the same 16-thread `rna`-partition
cell, default-mode compress dropped from **438 s to 141 s (~3.1× faster) on the
noisy ONT run** and from **493 s to 450 s (~9%) on `ecoli_hifi`**, at identical
output. Then the **anchor-restricted tiler coder (#231)** — which aligns only
inter-anchor gaps rather than re-deriving the whole read — made ONT another **1.54×
faster at 16 threads** *and* smaller (the ratio gains in the table above). The
`--max` ONT point spends best-of-4 tiling references for its best-ratio **3.063×**,
now ahead of CoLoRd; the deep sequence lever is paid only when asked for.

## Lossy quality (fqxv binning)

Quality binning is the big lever, especially on long reads. `reorder-bin*` stacks
read-reordering on top of binning.

| dataset          | bin8 | bin4 | bin2 | reorder-bin2 |
|------------------|-----:|-----:|-----:|-------------:|
| rnaseq_novaseq   | 9.94 | 9.94 | 11.37|        33.98 |
| rnaseq_fullrange |10.84 |14.49 | 15.37|        31.95 |
| ecoli_ont (binont)| 7.17|   —  |   —  |            — |

`ecoli_ont --quality-bin ont` improved 5.83 → 6.68 (the tiling codec on its
sequence stream) → **7.17** this cycle, as the per-block quality quantizer (#235)
shrank the binned quality stream from ~53M to 47M. On the lossy long-read points
the sequence stream still dominates the archive, so the ONT sequence lever matters
more there than in the lossless totals above.

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
least, since the long-read sequence stream dominates that archive — though it
improved from 0.825 to **0.721** as the tiling codec shrank that stream. Refreshed
2026-07-23 at HEAD `b900125`; the short-read points are byte-stable across the
last several release cycles. A grouped size chart of this table is on the
[Benchmarks](../docs/benchmarks.md) page.
