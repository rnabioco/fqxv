# Benchmarks

Field comparison against general-purpose and FASTQ-specific compressors on two
real Illumina RNA-seq runs spanning both quality regimes. Produced by the
harness in `bench/`, 4M-read subsets, 44 threads per tool, single node.
[Long reads](#long-reads-ont-pacbio) are measured separately, against CoLoRd.

!!! note
    Ratios and throughput vary with data type (WGS vs RNA-seq), coverage, and
    quality regime (full-range vs binned). Run the `bench/` harness on your own
    libraries before drawing conclusions.

`fqxv` has three operating points, all **lossless** on sequence and quality and,
uniquely among the tools here, **deterministic** (byte-identical output
regardless of thread count):

- **`fqxv`** (default) — fast, block-local, preserves everything.
- **`--max`** — deepest context + read reordering, still preserves read order and
  names exactly (a permutation is stored to restore order).
- **`--order shuffle`** — the reorder-lossy best-ratio point: reads come back as
  the same sequence+quality multiset but **renumbered** (original names and order
  discarded), so no permutation is stored. This is the same trade SPRING makes by
  default, and it is where `fqxv` is smallest.

## NovaSeq 6000, binned quality (DRR174812) — 4M reads, 572 Mbase

| Tool | Ratio | Compress | Decompress | Lossless | Deterministic |
| --- | ---: | ---: | ---: | :---: | :---: |
| **`fqxv --order shuffle`** | **23.9×** | 31 MB/s | 169 MB/s | seq+qual (renumbered) | **yes** |
| SPRING | 21.9× | 70 MB/s | 153 MB/s | reordered | — |
| **`fqxv --max`** | 20.2× | 37 MB/s | 134 MB/s | **yes (order-preserving)** | **yes** |
| `fqxv` (default) | 9.8× | **129 MB/s** | 103 MB/s | **yes** | **yes** |
| fqz_comp | 9.6× | 109 MB/s | 74 MB/s | **no** | — |
| zstd -19 | 9.4× | 11 MB/s | 411 MB/s | yes | — |
| xz -9 | 8.9× | 9 MB/s | 326 MB/s | yes | — |
| gzip | 5.0× | 274 MB/s | 184 MB/s | yes | — |

## GAIIx, full-range quality (SRR453566) — 4M reads, 385 Mbase

| Tool | Ratio | Compress | Decompress | Lossless | Deterministic |
| --- | ---: | ---: | ---: | :---: | :---: |
| **`fqxv --order shuffle`** | **10.27×** | 55 MB/s | 150 MB/s | seq+qual (renumbered) | **yes** |
| SPRING | 9.81× | 60 MB/s | 131 MB/s | reordered | — |
| **`fqxv --max`** | 9.73× | 55 MB/s | 141 MB/s | **yes (order-preserving)** | **yes** |
| `fqxv` (default) | 7.12× | **139 MB/s** | 112 MB/s | **yes** | **yes** |
| fqz_comp | 7.05× | 92 MB/s | 71 MB/s | **no** | — |
| xz -9 | 6.07× | 6 MB/s | 267 MB/s | yes | — |
| zstd -19 | 5.87× | 11 MB/s | 459 MB/s | yes | — |
| gzip | 3.52× | 272 MB/s | 163 MB/s | yes | — |

## Reading the numbers

- **`fqxv --order shuffle` is the smallest lossless compressor here — it beats
  SPRING on both datasets** (23.9× vs 21.9× on NovaSeq, 10.27× vs 9.81× on
  full-range), under the same rules SPRING plays by (reorder + renumber). SPRING
  keeps the original read *names*; `fqxv --order shuffle` renumbers, but the name
  stream is a rounding error either way (~a few KB on a 60 MB archive), so this is
  an apples-to-apples ratio win from `fqxv`'s sequence and quality coding.
- **`fqxv --max` is the best-ratio *fully* lossless option** — it additionally
  preserves original read order and names, which SPRING does not. That guarantee
  costs a read-order permutation: on NovaSeq the `--max` archive is
  names 0.1 MB + sequence 42.7 MB (of which ~11 MB is the permutation) +
  quality 30.7 MB. Dropping just that permutation (via `--order shuffle`) is what
  takes the sequence stream from 42.7 MB to 31.3 MB and puts `fqxv` ahead of
  SPRING.
- **Determinism and verified losslessness are `fqxv`'s alone here.** Every `fqxv`
  archive round-trips its exact content (sequence+quality multiset, and names in
  the order-preserving modes) and is byte-identical across thread counts;
  `fqz_comp` fails the content round-trip in this harness on both sets.
- **The default mode trades ratio for speed cleanly**: still ahead of fqz_comp /
  zstd / xz on ratio, at >120 MB/s — an order of magnitude faster to compress
  than zstd -19 or xz -9.
- **Lossy quality is a separate tier.** `fqxv --quality-bin bin8|bin4|bin2` and
  SPRING's `ill_bin` / `binary` modes quantize quality for more ratio; compare
  those against each other, not against the lossless rows above.

## Long reads (ONT / PacBio)

Long-read FASTQ compresses correctly today — read lengths are `u32` end to end,
blocks are cut by a raw-byte budget so ragged reads still parallelize, and the
short-read reorder codec auto-disables. The tables below measure each stream
per-read against CoLoRd `-q org` (its lossless quality mode); the `fqxv -l9`
rows' seq column is the **within-read** sequence baseline. The cross-read overlap codec that
now ships for long reads (see [The sequence lever](#the-sequence-lever-wired))
replaces that baseline and codes the HiFi sequence stream ~6× smaller.

fqxv's columns come from `fqxv info`, which reports its real streams. CoLoRd has
no such split, so its columns are taken by difference — `-q none` discards
quality, and quality is `org - none` — which makes the rows additive by
construction. `-q none` still carries names and container overhead, so CoLoRd's
non-quality column is an upper bound on its sequence stream. See
[Long-read support](design/longread.md) for the method and `bench/scripts/colord_split.sh`
to re-run it.

All sizes are `MiB` (base-1024, matching `fqxv info` and the harness `report.py`).

**`ecoli_ont`** (DRR205413, 301 Mbase, mean Q≈11.5 — noisy older basecaller):

| Tool | Total | Non-quality | Quality | Non-quality bits/base |
| --- | ---: | ---: | ---: | ---: |
| CoLoRd `-q org` | **197.9M** | **31.4M** | **166.5M** | 0.88 |
| `fqxv` (binmix qual) | 222.6M | 58.8M (seq only) | 163.7M | 1.64 |

**`ecoli_hifi`** (SRR11434954 subset, 1.55G bases, mean Q≈27, ~300×):

| Tool | Total | Non-quality | Quality | Non-quality bits/base |
| --- | ---: | ---: | ---: | ---: |
| CoLoRd `-q org` | **697.7M** | **13.4M** | **684.3M** | 0.069 |
| `fqxv` (binmix qual) | 768.3M | 126.3M (seq only) | 641.8M | 0.685 |

Two facts hold on both platforms:

- **Quality now beats CoLoRd on both platforms.** fqxv's binary-decomposition
  context-mixing quality coder codes the HiFi quality stream to **641.8M vs
  CoLoRd's 684.3M** (~6% smaller) and ONT to **163.7M vs 166.5M** (~2% smaller).
  This quality win carries the **HiFi lossless total ahead of CoLoRd** (see the
  repository's `bench/RESULTS.md`; fqxv 4.68× vs CoLoRd 4.44× on the full file);
  on ONT the much larger sequence deficit below keeps CoLoRd ahead on total.
- **The remaining lossless gap is the sequence stream**, and it widens with
  coverage: on HiFi it is now the *only* deficit (quality is a credit), so closing
  it widens fqxv's lead; on ONT it is the larger of the two. At ~300× the same locus
  is read hundreds of times — CoLoRd codes each read against a similar earlier read,
  while fqxv's per-block overlap reference is re-stored every block (HiFi) or the
  exact seed anchors don't survive the error rate (ONT).

### Lossy quality

`--quality-bin ont` cuts the ONT quality stream from 163.7M to 49.2M (3.4×) at
mean |Δ| 3.35. The tables work — but binning removes the stream fqxv is good at
and leaves the one it is not, so sequence becomes 57–62% of the lossy archive and
the gap above dominates the total. `colord-lossy` remains smaller overall
(ONT 73.2M vs 114.6M for `fqxv --quality-bin ont`). Match the table to the
platform: see [Lossy quality binning](cli/compress.md#lossy-quality-binning).

### The sequence lever (wired)

The cross-read overlap codec `fqxv-lroverlap` (minimizers → overlaps → layout →
consensus → per-read banded edit script → rANS) is now the container's sequence
codec for long reads. It is selected automatically when a block's mean read
length exceeds 500 bp; both it and the within-read order-k model are coded for
every long-read block and the **smaller is kept**, so it never regresses a block
(low coverage, a large genome, or too few reads all fall back to order-k). Output
is byte-identical across thread counts, and every archive round-trips (enforced
by the per-block content digest and `compress --verify`).

**It operates per block.** Each 256 MiB block self-assembles its own reference,
which keeps the container's blocked parallelism, per-block random access, and
thread-determinism intact — but caps within-block coverage at about
`256 MiB / genome`, not the whole file's depth. Measured through a real,
round-trip-verified archive of the whole `ecoli_hifi` file (120k reads, 1.55
Gbase, 6 blocks at ~52× each), the sequence stream drops from **0.653 to 0.107
bits/base — 6.1× smaller** than the within-read model (total archive 4.04×).

That 0.107 does **not** reach CoLoRd's 0.068, and the difference is precisely the
per-block coverage cap: CoLoRd assembles across the whole 300× file, a per-block
reference sees ~52×. The measurement harness
(`crates/fqxv-lroverlap/examples/encode.rs`) codes the whole file as one
reference and reaches **0.067 bits/base — parity with CoLoRd's 0.068** (strides
4–14 span 0.067–0.072, sample noise; the oracle with true-reference placement is
0.040). So the ceiling is real and demonstrated; closing the container's gap to
it means giving the codec more coverage per reference — larger long-read blocks
or a reference shared across blocks — without breaking the per-block invariants.

See [Long-read support](design/longread.md) for the full analysis.

## Round-trip fidelity (alignment level)

A compression ratio only counts if the reads survive intact. Beyond `fqxv`'s
internal round-trip digest, the `bench/scripts/bam_identity.sh` harness proves fidelity
*through a real analysis*: it aligns the original reads and the `fqxv`
round-tripped reads with `bwa mem` and compares the BAMs with an order-independent
multiset digest (`bench/tools/bamcmp.rs`).

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
    does by default and in its `--max` mode on real data. This is an aligner
    property, not a `fqxv` effect: the reads going in are byte-identical.

## Codecs

Every codec is clean-room, pure Rust — there is **no external/C compressor
dependency**. The reordered sequence path assembles a shared reference over the
clustered reads and codes each read as a position on it; the reference itself is
2-bit-packed and coded with a clean-room LZMA (or the order-*k* `fqxv-seq` model,
whichever is smaller), split into fixed blocks compressed in parallel. Quality is
a fqzcomp-class context model coded over only the quality levels that occur.
Assembly, the overlap-merge, the reference coder, and the per-block codecs all fan
out with `rayon`, and every stage uses fixed, thread-independent boundaries so
output never depends on the thread count.
