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
| **`fqxv --order shuffle`** | **23.9×** | 13 MB/s | 148 MB/s | seq+qual (renumbered) | **yes** |
| SPRING | 21.9× | 70 MB/s | 136 MB/s | reordered | — |
| **`fqxv --max`** | 20.2× | 19 MB/s | 124 MB/s | **yes (order-preserving)** | **yes** |
| `fqxv` (default) | 9.8× | **123 MB/s** | 116 MB/s | **yes** | **yes** |
| fqz_comp | 9.6× | 105 MB/s | 73 MB/s | **no** | — |
| zstd -19 | 9.4× | 11 MB/s | 415 MB/s | yes | — |
| xz -9 | 8.9× | 10 MB/s | 321 MB/s | yes | — |
| gzip | 5.0× | 374 MB/s | 194 MB/s | yes | — |

## GAIIx, full-range quality (SRR453566) — 4M reads, 385 Mbase

| Tool | Ratio | Compress | Decompress | Lossless | Deterministic |
| --- | ---: | ---: | ---: | :---: | :---: |
| **`fqxv --order shuffle`** | **10.25×** | 18 MB/s | 147 MB/s | seq+qual (renumbered) | **yes** |
| SPRING | 9.81× | 58 MB/s | 131 MB/s | reordered | — |
| **`fqxv --max`** | 9.71× | 15 MB/s | 125 MB/s | **yes (order-preserving)** | **yes** |
| `fqxv` (default) | 7.12× | **124 MB/s** | 116 MB/s | **yes** | **yes** |
| fqz_comp | 7.05× | 92 MB/s | 70 MB/s | **no** | — |
| xz -9 | 6.07× | 6 MB/s | 288 MB/s | yes | — |
| zstd -19 | 5.87× | 11 MB/s | 390 MB/s | yes | — |
| gzip | 3.52× | 363 MB/s | 164 MB/s | yes | — |

## Reading the numbers

- **`fqxv --order shuffle` is the smallest lossless compressor here — it beats
  SPRING on both datasets** (23.9× vs 21.9× on NovaSeq, 10.25× vs 9.81× on
  full-range), under the same rules SPRING plays by (reorder + renumber). SPRING
  keeps the original read *names*; `fqxv --order shuffle` renumbers, but the name
  stream is a rounding error either way (~a few KB on a 60 MB archive), so this is
  an apples-to-apples ratio win from `fqxv`'s sequence and quality coding.
- **`fqxv --max` is the best-ratio *fully* lossless option** — it additionally
  preserves original read order and names, which SPRING does not. That guarantee
  costs a read-order permutation: on NovaSeq the `--max` archive is
  names 0.1 MB + sequence 40.7 MB (of which ~11 MB is the permutation) +
  quality 29.3 MB. Dropping just that permutation (via `--order shuffle`) is what
  takes the sequence stream from 40.7 MB to 29.8 MB and puts `fqxv` ahead of
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
short-read reorder codec auto-disables. Long-read blocks pick a sequence codec per
block from several candidates and keep the smallest (see
[the sequence codecs](#the-sequence-codecs) below), so no block ever regresses
below the within-read order-k model.

The tables below measure each stream against CoLoRd `-q org` (its lossless quality
mode). fqxv's columns come from `fqxv info`, which reports its real streams; CoLoRd
has no such split, so its columns are taken by difference — `-q none` discards
quality, and quality is `org - none`. `-q none` still carries names and container
overhead, so CoLoRd's non-quality column is an *upper bound* on its sequence
stream, while fqxv's is the sequence stream alone. Sizes are decimal MB
(bytes ÷ 10⁶); bits/base is over the true base count.

**`ecoli_ont`** (DRR205413, 301 Mbase, 21,140 reads, mean Q≈11.5 — noisy older
basecaller):

| Tool | Total | Non-quality | Quality | Non-quality bits/base |
| --- | ---: | ---: | ---: | ---: |
| CoLoRd `-q org` | 197.9M | **31.4M** | 166.5M | 0.84 |
| `fqxv -l9` | **197.1M** | 33.4M (seq) | **163.7M** | 0.89 |

**`ecoli_hifi`** (SRR11434954 subset, 1.55G bases, mean Q≈27, ~300×):

| Tool | Total | Non-quality | Quality | Non-quality bits/base |
| --- | ---: | ---: | ---: | ---: |
| CoLoRd `-q org` | 697.7M | 13.4M | 684.3M | 0.069 |
| `fqxv` | **656.0M** | **12.6M (seq)** | **641.8M** | **0.065** |

Two facts hold on both platforms:

- **Quality beats CoLoRd on both platforms.** fqxv's binary-decomposition,
  context-mixing quality coder codes the HiFi quality stream to **641.8M vs
  CoLoRd's 684.3M** (~6% smaller) and ONT to **163.7M vs 166.5M** (~2% smaller).
  On HiFi that quality win carries the **lossless total ahead of CoLoRd** (656.0M
  vs 697.7M, **4.73× vs 4.44×** on the full file); on ONT fqxv now edges ahead on
  the **total** too — **197.1M vs 197.9M** (3.06× vs 3.05×).
- **The sequence stream is no longer a blanket deficit.** On HiFi Sequel II it is a
  *credit* (12.6M vs 13.4M) — the same locus is read hundreds of times at ~300×,
  and fqxv now assembles one whole-file consensus and codes each read against it. On
  ONT it is down to 33.4M vs 31.4M — ~6% behind on sequence (the anchor-restricted
  coder aligns only inter-anchor gaps), which the quality credit now more than
  closes. The one regime it still trails materially is
  ordinary-coverage HiFi (Revio WGS), which the raw-LZMA path below addresses.
  Each lever is described under [the sequence codecs](#the-sequence-codecs).

**Compress speed.** The long-read codecs are the compute-heavy part of `fqxv`, and
recent work cut that cost without touching the output: on a 16-thread node,
default-mode compress of the noisy ONT run is now **~3× faster** (skipping an
always-discarded consensus candidate and the redundant shared-reference assembly
on Nanopore, and de-packing the alignment traceback), and HiFi is ~9% faster —
same archives, same ratios. On top of that the anchor-restricted tiler coder (#231)
makes ONT another **1.54× faster at 16 threads** *and* smaller (it aligns only
inter-anchor gaps). The deepest sequence lever (best-of-N tiling references) stays
gated to `--max`, so the default stays fast and `--max` buys the extra ratio only
when you ask for it.

### Lossy quality

`--quality-bin ont` cuts the ONT quality stream by ~3.5× (mean |Δ| 3.35). Binning
removes the stream fqxv is best at and leaves sequence dominant, so the platform
match matters — see [Lossy quality binning](cli/compress.md#lossy-quality-binning).
`colord-lossy` is still smaller on the aggressive long-read lossy points.

### The sequence codecs

Long-read blocks are coded by up to four sequence methods and the **smallest is
kept per block**, so selection is automatic and a block can never regress below the
within-read order-k model. All are clean-room Rust, byte-identical across thread
counts, and every archive round-trips (per-block content digest plus
`compress --verify`). Three of them are the long-read levers, each tuned to a
different coverage/error regime:

- **Shared whole-file overlap reference** — the high-coverage HiFi lever
  (`fqxv-lroverlap`: minimizers → overlaps → consensus → per-read banded edit
  script → rANS). The container assembles **one** consensus over the whole file,
  stores it once in a framed region before the first block, and codes every block's
  reads against that frozen frame — so a read codes identically no matter which
  block holds it, and blocks stay independently decodable. On `ecoli_hifi` (120k
  reads, 1.55 Gbase, ~300×) this puts the sequence stream at **0.065 bits/base**,
  past CoLoRd's 0.069 — a ~10× shrink over the within-read model. Storing the
  reference once instead of per block was the whole gap.
- **Raw large-window LZMA** — the ordinary-coverage lever, and the biggest change
  this cycle. A real genome at typical (not 300×) coverage carries **exact**
  cross-read matches that neither the within-read model nor a per-block voted
  consensus can reach. On the modern **Revio WGS** run (`hifi_revio_wgs`, 1.35
  Gbase) it cuts the sequence stream from **1.391 to 0.683 bits/base**, taking the
  whole archive from **9.7× to 17×** — from the worst lossless result in the suite
  to a strong second behind CoLoRd (18.8×), now ahead of both zstd -19 and xz -9.
- **Multi-reference tiling + anchor-restricted coding** — the Nanopore lever. Each
  read is coded against earlier *raw* reads with best-of-N reference selection
  (engaged at `-l9`/`--max`), which pays where a read sits closer to another single
  read than to a divergent consensus — the ~10% ONT-error regime. Rather than
  re-aligning the whole read, the coder seeds minimizers, chains them, and aligns
  only the short inter-anchor gaps (CoLoRd-style) — both faster and tighter.
  Together they cut `-l9` ONT sequence to **0.887 bits/base**, bringing the ONT
  total **ahead of CoLoRd** (3.06× vs 3.05×).

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
