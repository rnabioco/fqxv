# Benchmark results — 2026-07-20

Snapshot from the unified parallel harness (`submit_parallel.sh`) on the Bodhi
`rna` partition: 9 datasets × platform-appropriate tools = 150 cells, each fanned
out one-tool-per-node, all COMPLETED. Covers short-read Illumina (MiSeq,
NovaSeq6000, GAIIx) and MGI/BGISEQ, plus long-read ONT (MinION) and PacBio HiFi
on **both** Sequel II and Revio. Ratios are deterministic (thread-count
independent); `rt=yes` means the round-trip content digest matched (lossless, or
lossy-expected for the binned points). Reproduce with `bash submit_parallel.sh`.

> `fqzcomp5` recorded **no rows**: the binary is not installed, so `run_bench.sh`
> skipped it (`[miss] fqzcomp5`) rather than recording a failure. It was equally
> absent from the 2026-07-18 run, so this is a standing gap, not a regression —
> run `build_tools.sh` to populate it.

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
| mgi_bgiseq_hard  | 3.80 |     3.78 |         3.89 |**4.05**|   3.14 | 3.26| 2.57 |

Both MGI rows are new, and they pull in opposite directions. On **22 bp miRNA
reads** fqxv wins by a wide margin (15.14 vs SPRING's 12.22) — at that length the
read names dominate the archive and the tokenizer is the whole game. On
**`mgi_bgiseq_hard`** (BGISEQ-500, 137 bp, 36-symbol quality) SPRING wins, and
it is the one dataset where `fqxv-max` (3.78) is *worse* than plain `fqxv`
(3.80): the reorder pass costs more than it recovers on data with little
cross-read redundancy. This is the hardest short-read case in the corpus at 5.32
bits/base, and it was invisible before because no MGI run was benchmarked.

## Long-read, lossless

**The previous "fqxv leads on HiFi" claim does not survive on Revio.** It was
measured only on `ecoli_hifi`, and adding two Revio runs shows the lead was an
artifact of that dataset. All rows round-trip losslessly.

| dataset             |     fqxv | fqxv9 |    colord | zstd19 |   xz9 | gzip |
|---------------------|---------:|------:|----------:|-------:|------:|-----:|
| hifi_revio_amplicon |**21.79** | 21.79 |     19.88 |  15.12 | 14.45 | 8.96 |
| hifi_revio_wgs      |     9.72 |  9.72 | **18.76** |  13.02 | 12.75 | 9.28 |
| ecoli_hifi          | **4.73** |  4.73 |      4.44 |   3.83 |  3.85 | 2.27 |
| ecoli_ont           |     2.79 |  2.79 |  **3.05** |   2.38 |  2.49 | 1.94 |

**On `hifi_revio_wgs` fqxv is beaten by nearly 2× — and not only by CoLoRd.**
`zstd19` (13.02) and `xz9` (12.75) both beat it too; fqxv (9.72) barely clears
gzip (9.28), making it the second-worst tool in the lossless set on the most
mainstream modern PacBio case there is. This is the single worst result in the
suite and it had been invisible.

The per-stream split explains it exactly, and shows why one HiFi dataset was
never enough:

| dataset             | seq (b/base) |  seq (bytes) | quality (bytes) | quality share |
|---------------------|-------------:|-------------:|----------------:|--------------:|
| ecoli_hifi          |        0.065 |   12,643,309 |     641,788,411 |       **98%** |
| hifi_revio_amplicon |        0.083 |    6,299,312 |      49,187,324 |           86% |
| hifi_revio_wgs      |    **1.391** |  234,726,792 |      44,876,883 |       **16%** |
| ecoli_ont           |        1.283 |   48,292,740 |     163,660,203 |           41% |

`ecoli_hifi` is **98% quality by bytes**, so its result is essentially a
measurement of full-range quality coding — the one thing fqxv does better than
CoLoRd — and its 0.065 b/base sequence stream is an artifact of ~300× coverage of
a single 4.6 Mb genome, where cross-read redundancy is enormous. Revio flips
both halves: the 7-symbol quality alphabet (Phred 3–40) collapses quality to 16%
of the archive, and a real genome at ordinary coverage gives the overlap codec
almost nothing to exploit. The result is a **1.391 b/base** sequence stream —
*worse than ONT's 1.283*, on data with a fraction of ONT's error rate.

So the honest summary is: fqxv's quality coder is genuinely strong, and its
sequence codec only looks strong when coverage is high enough to make reads
redundant (300× E. coli, or amplicons). **The lever is the sequence stream on
ordinary-coverage long reads**, which is now measurable on `hifi_revio_wgs`
rather than hidden. `ecoli_ont` at 5.75 bits/base sits inside the corpus ONT-WGS
spread (5.53–6.21), so the ONT entry is representative; its deficit is likewise
entirely sequence, bounded by the quality of the assembled consensus reads are
coded against.

> **Known regression on ONT.** The ONT total moved the *wrong* way with the shared
> whole-file reference: 2.827 → 2.791. The reference shrinks the ONT edit streams
> by 1.58 MB but costs 4.37 MB to store, a net loss of ~2.8 MB. It is adopted
> anyway because the never-worse gate compares the reference layout against
> **order-k**, not against the per-block overlap layout it replaces — and order-k
> (~1.8 b/base) is a much weaker bar. HiFi is unaffected: there the reference costs
> 1.36 MB and saves 7.2 MB. Tracked as a follow-up; the fix is to gate against the
> better of the two layouts.

## Lossy quality (fqxv binning)

Quality binning is the big lever, especially on long reads. `reorder-bin*` stacks
read-reordering on top of binning.

| dataset          | bin8 | bin4 | bin2 | reorder-bin2 |
|------------------|-----:|-----:|-----:|-------------:|
| rnaseq_novaseq   | 9.94 | 9.94 | 11.37|        33.98 |
| rnaseq_fullrange |10.84 |14.49 | 15.37|        31.95 |
| ecoli_ont (binont)| 5.83|   —  |   —  |            — |

## fqxv archive vs native NCBI `.sra` (lossless `max` regime)

`fqxv/.sra` < 1 means fqxv beats the `.sra` archive; it wins on every platform.

| accession  | platform    | fqxv / .sra |
|------------|-------------|------------:|
| DRR174812  | NovaSeq6000 |       0.331 |
| SRR453566  | GAIIx       |       0.509 |
| SRR2627175 | MiSeq       |       0.538 |
| DRR205413  | ONT-MinION  |       0.825 |

NovaSeq wins most (its quality is pre-binned, so the lossless point is already
compact); ONT wins least, since the long-read sequence stream dominates that
archive. The ONT row predates the edit-stream context coding, the shared
whole-file reference, and syncmer seeding, so it should improve on a re-run.
