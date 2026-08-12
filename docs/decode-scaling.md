# Decode scaling

How fast does the archive turn back into FASTQ when you throw cores at it?

The BINSEQ paper ([Teyssier & Dobin 2026, *PLOS Comput Biol*,
doi:10.1371/journal.pcbi.1014181](https://doi.org/10.1371/journal.pcbi.1014181))
makes the case that a sequencing format is an *analysis substrate*, not just an
archive: consumers reading gzip-compressed FASTQ saturate at 2–4 threads —
DEFLATE inflation is inherently serial, so allocating more cores buys nothing —
while block-based binary formats keep scaling near-linearly. `fqxv` has the
structural properties BINSEQ builds on (independent, parallel-decodable blocks;
a footer index for random access; paired mates in one file) *plus*
FASTQ-specialized compression, but the [benchmarks](benchmarks.md) page only
reports ratio and single-configuration throughput. This page measures the
scaling story directly: full-decompression throughput versus thread count, for
`fqxv` against `gzip`, `pigz`, and `zstd` decoding the same FASTQ.

## Results

Throughput is MB/s of decompressed FASTQ (decimal MB, median of 3 runs). Every
cell decoded to the byte-identical output — the piped byte count is checked for
every tool at every thread count. "ratio" is compressed-size ÷ plain FASTQ,
measured on these exact inputs (it differs slightly from `bench/RESULTS.md`,
which divides by the *raw* SRA dump before `+`-line normalization).

### NovaSeq 6000, paired 2×151 bp (DRR174812) — 10.09 GB, 29.4 M reads

| tool | size | ratio | 1 thr | 2 | 4 | 8 | 16 | 32 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `fqxv` (default) | 1.168 GB | 8.6× | 29 | 48 | 76 | 112 | 179 | **238** |
| `fqxv --max` | **0.496 GB** | **20.3×** | 65 | 101 | 139 | 172 | 194 | 194 |
| `gzip -dc` | 2.137 GB | 4.7× | 172 | — | — | — | — | — |
| `pigz -dc -p N` | 2.137 GB | 4.7× | 261 | 338 | 338 | 339 | 338 | 337 |
| `zstd -dc -T N` | 1.133 GB | 8.9× | 777 | 780 | 778 | 780 | 778 | 776 |

`fqxv` decode scales **8.1×** from 1→32 threads (29 → 238 MB/s), passing the
serial `gzip` baseline between 8 and 16 threads. `--max` — the archive you
would actually keep, **2.3× smaller than `zstd -19 --long`** — starts *faster*
than the default point (the reordered sequence stream leaves the range coders
far less work per read) and matches serial `gzip` at 8 threads, passing it by
16. `pigz` gains once from its separate read/write/CRC threads (1.3× over one
thread) and is then flat — its own documentation is clear that DEFLATE
decompression cannot be parallelized. `zstd` ignores `-T` on decode entirely
(single-threaded by design), it is just a very fast serial decoder.

### MiSeq, paired 2×301 bp, full-range quality (SRR2627175) — 2.34 GB, 4.4 M reads

| tool | size | ratio | 1 thr | 2 | 4 | 8 | 16 | 32 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `fqxv` (default) | 568 MB | 4.1× | 21 | 34 | 56 | 87 | 89 | 88 |
| `fqxv --max` | **393 MB** | **6.0×** | 37 | 63 | 96 | 128 | 154 | 168 |
| `gzip -dc` | 881 MB | 2.7× | 125 | — | — | — | — | — |
| `pigz -dc -p N` | 881 MB | 2.7× | 165 | 198 | 198 | 198 | 198 | 197 |
| `zstd -dc -T N` | 536 MB | 4.4× | 570 | 577 | 574 | 572 | 577 | 574 |

The default point's curve stops dead at 8 threads — this 2.34 GB archive has
only **5 blocks** (see below), so there is nothing for threads 9–32 to do. The
`--max` archive reorders into 17 smaller blocks and keeps climbing to
168 MB/s — 1.3× faster than serial `gzip` *and* 1.4× smaller than `zstd -19`.
Full-range quality also makes each decoded byte more expensive than NovaSeq's
4-symbol binned quality.

### ONT MinION, ~14 kb reads (DRR205413) — 0.60 GB, 21,140 reads

| tool | size | ratio | 1 thr | 2 | 4 | 8 | 16 | 32 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `fqxv` (default) | 204 MB | 3.0× | 8.3 | 16 | 27 | 45 | 62 | **67** |
| `fqxv --max` | **197 MB** | **3.1×** | 9.5 | 11 | 11 | 11 | 11 | 11 |
| `gzip -dc` | 311 MB | 1.9× | 105 | — | — | — | — | — |
| `pigz -dc -p N` | 311 MB | 1.9× | 152 | 180 | 180 | 180 | 180 | 180 |
| `zstd -dc -T N` | 254 MB | 2.4× | 353 | 353 | 351 | 352 | 351 | 353 |

Here `fqxv` wins the ratio comparison outright (20% smaller than
`zstd -19 --long` already at the default level; `--max` adds 3% more), and the
default archive now *scales*: **8.1×** from 1→32 threads (73.1 s → 9.0 s full
decode), where this same file used to be flat at ~11 MB/s at every thread
count. Two changes compose to unlock that. The 64 MiB Nanopore block budget
(#275) cuts this file
into **5 blocks** instead of 2, and chunked quality coding (#279, context mode
5, v0.7.0+) splits each block's quality stream — 92% of the decode work, and
until now one serial adaptive-coder pass per block — into a warmup prefix plus
7 chunks that decode in parallel once the bases exist. Together they cost
+1.6% of archive size versus the old serial default. The curve
bends at 16 threads not for lack of tasks (5 blocks × ~8 quality tasks) but
because concurrent chunk decoders random-walk ~336 MiB context models —
memory parallelism, not cores, is the ceiling. `--max` keeps the
smallest-archive contract: maximal blocks, serial quality, byte-identical to
what v0.6.x wrote, and flat at ~11 MB/s; `--block-reads` still trades block
granularity explicitly. A mode-5 archive needs a v0.7.0+ reader (older
binaries refuse it loudly by name).

## Why the curves look like this

**The unit of `fqxv` decode parallelism is the block.** Blocks are independent
by construction (the same property that makes compression deterministic and
random access possible), so decode fans out block-per-task — but it cannot fan
out further than the archive has blocks. At the default level 5, a block is
1 Mi reads:

| dataset | blocks (default) | saturates at | blocks (`--max`) | saturates at |
| --- | ---: | --- | ---: | --- |
| DRR174812 | 29 | ~32 threads | 113 | ~16 threads |
| SRR2627175 | 5 | ~8 threads | 17 | ~32 threads |
| DRR205413 | 5 | ~16 threads (chunked) | 2 | ~2 threads |

The `--max` short-read path (read reordering) cuts its own smaller blocks
(~260 k reads), which is what un-sticks SRR2627175's 5-block ceiling. With
blocks to spare (DRR174812's 113), the `--max` curve still bends at ~16
threads: restoring the original read order from the stored permutation and
reassembling interleaved output is work the block fan-out cannot hide, and
the memory-bandwidth ceiling (below) applies sooner because each of its
threads decodes faster. The Nanopore default block budget is 64 MiB (#275) —
5 blocks on this file — while ONT `--max` keeps maximal blocks (2).

Block count used to be the *whole* lever for long reads, because long-read
blocks could not be un-serialized from the inside: quality is coded
*conditioned on the decoded bases* (`fqzcomp` MODE_SEQ), so a block's sequence
must finish decoding before its quality — 92% of the work — can start, and
that quality stream was one adaptive coding pass over the whole block. The
on-disk format change that fixes the second half of that shipped in v0.7.0
(quality mode 5, issue #278/#279; design and measured costs in
`docs/design/parallel-decode.md`): long-read quality is now coded as a serial
warmup plus K−1 chunks that decode in parallel once the bases exist, so
within-block quality fan-out composes with the block-level parallelism — the
ONT table above measures exactly that composition. The sequence stream itself
(the tiler/overlap read-chains) still decodes serially per block. Short-read
archives never chunk (mode 5/6 is platform-gated to long reads), so the
Illumina tables measure archives that are byte-identical either side of the
change.

**`gzip`'s flat line is a property of the format, not a slow implementation.**
A `.gz` member is one DEFLATE stream whose back-references chain each byte to
the bytes before it; nothing can be inflated before what precedes it. `pigz`
parallelizes *compression*, but on decode it can only offload read/write/CRC —
the ~1.3× step from 1→2 threads above, then flat forever. This is exactly the
2–4-thread FASTQ saturation BINSEQ measures, reproduced here at the
decompression layer itself.

**`zstd` decode is serial too — just fast.** `-T` affects compression only;
decode of a standard single-frame `.zst` is one thread, and its ~0.35–0.8 GB/s
makes it the raw-decode-speed champion in this field. If pure decode
throughput on a few threads is the only thing that matters, `zstd` is the
honest recommendation — what it gives up is the ratio (`fqxv --max` is 2.3×
smaller on NovaSeq, 1.4× on MiSeq, 1.3× on ONT), record-level random access,
and the container's verify/determinism guarantees.

**`fqxv` buys its ratio with per-byte decode work.** Sequence and quality come
back through adaptive context-model range coders — inherently serial per block,
several times more work per byte than DEFLATE. One `fqxv` thread is therefore
slower than one `gzip` thread, and parallel blocks are how the archive earns it
back: with enough blocks (DRR174812) it passes the serial-gzip ceiling at
~16 threads and is still climbing at 32. At high thread counts the range-coder
working set makes decode increasingly memory-bandwidth-bound, so the curve
bends (179 → 238 MB/s from 16 → 32 threads, not 2×) — more cores past that
point buy little.

## Method

Produced by `bench/scripts/decode_scaling.sh` (submit
`bench/slurm/decode_scaling.sbatch`); raw numbers with full run metadata in
[`docs/charts/decode_scaling.tsv`](charts/decode_scaling.tsv).

- **Datasets** (staged with `sracha`, plain FASTQ): see tables above; paired
  runs are the concatenated R1+R2 bytes. The `+` separator lines are
  normalized (`+SRR… length=…` → `+`) *before* compressing the byte-stream
  tools' inputs, because `fqxv` performs that normalization by design (its one
  documented deviation from byte-losslessness, shared with SPRING/fqz_comp) —
  without it no byte-count comparison could succeed.
- **Compressed inputs**, synthesized once from the same bytes: `.fqxv` /
  `.max.fqxv` from this release (paired runs as one spot-interleaved archive),
  `.fastq.gz` from `pigz -6` (the harness's gzip baseline of record; standard
  DEFLATE, decodable by `gzip`), `.fastq.zst` from `zstd -19 --long=27` (the
  `zstd19` point in `bench/scripts/toolsets.sh`).
- **Measurement**: each cell is a full decompression with output piped to
  `wc -c` — a real pipe, never `/dev/null` (tools can fast-path a null sink)
  and never the filesystem. hyperfine, median of 3 timed runs after 1 warm-up.
  The byte count is asserted identical across every tool and thread count
  (for paired archives `fqxv -Z` emits the same records interleaved, so the
  byte *count* matches the R1+R2 concatenation).
- **Commands** at thread count N: `fqxv decompress in.fqxv -Z --quiet
  --threads N`, `gzip -dc` (no thread knob), `pigz -dc -p N`,
  `zstd -dc --long=27 -T N`.
- **Hardware/versions**: Bodhi `rna`-partition nodes (Intel Xeon Gold 6240R),
  dedicated 32-CPU / 64 GB Slurm allocations (not node-exclusive), one job per
  table. The Illumina tables were measured with build `v0.6.2-8-g4d6b434`
  (whose short-read archives are byte-identical to v0.7.0 — mode 5 never
  applies to short reads); the ONT table with `v0.6.2-9-g102fa86`, the v0.7.0
  chunked-quality codebase. gzip 1.12, pigz 2.8, zstd 1.5.7,
  hyperfine 1.20.0.
