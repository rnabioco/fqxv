# Benchmarks

All numbers below are from a real Illumina RNA-seq run (SRR453566, full-range
quality), measured on a single CPU node. The harness that produces them lives in
`bench/` — see the repository `bench/README.md`.

!!! note
    These are early numbers on one dataset. Ratios vary with data type (WGS vs
    RNA-seq), coverage, and quality regime (full-range vs binned). Run the
    harness on your own libraries before drawing conclusions.

## Per stream, vs the C references

`fqxv`'s clean-room codecs match or beat the reference C tools stream-for-stream
(1M reads × 101 bp):

| Stream | fqxv | C reference | |
| --- | --- | --- | --- |
| Quality | **1.737** bits/byte | `fqz_comp` 1.735 | matches |
| Sequence | **1.247** bits/base (order-11) | `fqz_comp` 1.474 | beats |
| Names | 9.2× (6.19 MB) | `fqz_comp` ~18× (3.06 MB) | behind |

The sequence model already beats `fqz_comp`'s on this duplicate-rich data; the
quality model matches it; the name tokenizer is a first version with a clear path
to close the gap (per-column adaptive models).

## Whole archive, vs the field

On the same 1M-read file (319.8 MB raw / ~85 MB gzipped), compared to
general-purpose and FASTQ-specific tools:

| Tool | × smaller than gzip | round-trips |
| --- | --- | --- |
| SPRING | 2.61× | yes (read reordering) |
| fqzcomp5 | 2.20× | yes |
| **fqxv** | **1.84×** | yes |
| fqz_comp | 1.95× | yes |
| xz -9 | 1.63× | yes |
| zstd -19 --long | 1.57× | yes |
| gzip | 1.00× | — (baseline) |

`fqxv` is competitive with `fqz_comp`, with sequence already ahead and quality
matched; the current gap is the name tokenizer plus the block-size/ratio
trade-off (smaller blocks parallelize better but train the sequence model on
fewer reads). Closing those, plus the reorder work, is the path toward the
`SPRING`/`fqzcomp5` tier.

## Throughput

Compression ran at ~230 MB/s across 16 threads (blocks fan out with `rayon`);
decompression at a similar rate. The rANS entropy coder reaches ~200 MB/s
single-threaded on the scalar path.
