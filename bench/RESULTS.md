# Benchmark results — 2026-07-18

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

## Long-read, lossless — the open gap

CoLoRd still leads on lossless long-read: its cross-read sequence coding is the
edge (fqxv quality already matches; the sequence stream is the gap). fqxv numbers
match the standalone WFA-path runs exactly.

| dataset    | fqxv | fqxv9 | **colord** | zstd19 | xz9 | gzip |
|------------|-----:|------:|-----------:|-------:|----:|-----:|
| ecoli_hifi | 4.04 |  4.04 |   **4.44** |   3.83 |3.85 | 2.27 |
| ecoli_ont  | 2.68 |  2.68 |   **3.05** |   2.38 |2.49 | 1.94 |

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
compact); ONT wins least, since the long-read sequence stream — the CoLoRd gap
above — dominates that archive.
