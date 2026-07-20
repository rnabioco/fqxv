# Benchmark results — 2026-07-19

Snapshot from the unified parallel harness (`submit_parallel.sh`) on the Bodhi
`rna` partition. One run covers all platforms: short-read Illumina (MiSeq,
NovaSeq6000, GAIIx) plus long-read ONT (MinION) and PacBio HiFi, each dataset
fanned out one-tool-per-node. Ratios are deterministic (thread-count independent);
`rt=yes` means the round-trip content digest matched (lossless, or lossy-expected
for the binned points). Reproduce with `bash submit_parallel.sh`.

## Short-read, lossless

`fqxv --order shuffle` (rename-and-renumber, lossless) beats SPRING on every
dataset; `fqxv-max` is the order-preserving lossless point.

| dataset          | fqxv | fqxv-max | **fqxv-shuffle** | spring | zstd19 | xz9 | gzip |
|------------------|-----:|---------:|-----------------:|-------:|-------:|----:|-----:|
| rnaseq_novaseq   | 9.94 |    23.80 |        **29.16** |  25.27 |   9.59 | 9.01| 5.05 |
| rnaseq_fullrange | 7.60 |    10.20 |        **11.62** |  10.01 |   5.84 | 5.98| 3.71 |
| ecoli_miseq      | 4.92 |     7.41 |         **7.35** |   7.32 |   5.24 | 5.13| 2.94 |

## Long-read, lossless

**fqxv leads on HiFi.** Its quality coder (binary-decomposition context mixing)
beats CoLoRd's lossless quality, and the shared whole-file reference has now taken
the HiFi *sequence* stream past CoLoRd too — **0.065 bits/base against CoLoRd's
0.068**, so fqxv is ahead on both streams. **CoLoRd still leads on ONT**, though
not on quality: fqxv's ONT quality stream is the smaller of the two. The ONT
deficit is entirely sequence, and it is no longer a seeding problem — closed
syncmers select on a k-mer's own bases, so the overlap codec beats the within-read
order-k baseline rather than falling back to it. What bounds it now is the quality
of the assembled consensus the reads are coded against. All rows round-trip
losslessly.

| dataset    |      fqxv | fqxv9 |   colord | zstd19 | xz9 | gzip |
|------------|----------:|------:|---------:|-------:|----:|-----:|
| ecoli_hifi | **4.73** |  4.73 |     4.44 |   3.83 |3.85 | 2.27 |
| ecoli_ont  |      2.79 |  2.79 | **3.05** |   2.38 |2.49 | 1.94 |

> **The HiFi lead is measured only on Sequel II.** `ecoli_hifi` has a 93-symbol
> full-range quality alphabet and ~640 MB of its 656 MB archive is quality — so
> the row above is largely a statement about full-range quality coding, which is
> exactly the lever that beats CoLoRd here. Every PacBio run in the robustness
> corpus is a **Revio** run with a **7-symbol** alphabet (Phred 3–40); that
> collapses the quality stream and leaves sequence dominant, which is the stream
> CoLoRd has historically been stronger on. `ecoli_hifi` is also subsampled from
> ~5000x to ~300x, and at 3.39 bits/base it is *harder* than all twelve corpus
> PacBio runs (0.22–2.55). Revio is the current PacBio instrument, so treat
> "fqxv leads on HiFi" as unverified on modern PacBio data until the
> `hifi_revio_wgs` / `hifi_revio_amplicon` rows land. By contrast `ecoli_ont` at
> 5.75 bits/base sits squarely inside the corpus ONT-WGS spread (5.53–6.21), so
> the ONT panel is representative and needs no such caveat.

Per-stream, over the run's own bases:

| dataset    | seq (b/base) | quality (bytes) | reference frame |
|------------|-------------:|----------------:|----------------:|
| ecoli_hifi |    **0.065** |     641,788,411 |       1,356,344 |
| ecoli_ont  |        1.283 |     163,660,203 |       4,369,101 |

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
