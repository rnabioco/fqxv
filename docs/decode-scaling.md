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
| `fqxv` (default) | 1.168 GB | 8.6× | 29 | 47 | 76 | 111 | 181 | **236** |
| `fqxv --max` | **0.496 GB** | **20.3×** | 65 | 100 | 139 | 173 | 193 | 196 |
| `gzip -dc` | 2.137 GB | 4.7× | 171 | — | — | — | — | — |
| `pigz -dc -p N` | 2.137 GB | 4.7× | 259 | 336 | 334 | 335 | 334 | 334 |
| `zstd -dc -T N` | 1.133 GB | 8.9× | 778 | 778 | 778 | 774 | 777 | 776 |

`fqxv` decode scales **8.1×** from 1→32 threads (29 → 236 MB/s), passing the
serial `gzip` baseline between 8 and 16 threads. `--max` — the archive you
would actually keep, **2.3× smaller than `zstd -19 --long`** — starts *faster*
than the default point (the reordered sequence stream leaves the range coders
far less work per read) and beats serial `gzip` from 8 threads. `pigz` gains
once from its separate read/write/CRC threads (1.3× over one thread) and is
then flat — its own documentation is clear that DEFLATE decompression cannot
be parallelized. `zstd` ignores `-T` on decode entirely (single-threaded by
design), it is just a very fast serial decoder.

### MiSeq, paired 2×301 bp, full-range quality (SRR2627175) — 2.34 GB, 4.4 M reads

| tool | size | ratio | 1 thr | 2 | 4 | 8 | 16 | 32 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `fqxv` (default) | 568 MB | 4.1× | 21 | 34 | 57 | 88 | 88 | 88 |
| `fqxv --max` | **393 MB** | **6.0×** | 38 | 64 | 97 | 130 | 156 | 170 |
| `gzip -dc` | 881 MB | 2.7× | 123 | — | — | — | — | — |
| `pigz -dc -p N` | 881 MB | 2.7× | 162 | 193 | 193 | 193 | 193 | 194 |
| `zstd -dc -T N` | 536 MB | 4.4× | 572 | 575 | 575 | 573 | 563 | 573 |

The default point's curve stops dead at 8 threads — this 2.34 GB archive has
only **5 blocks** (see below), so there is nothing for threads 9–32 to do. The
`--max` archive reorders into 17 smaller blocks and keeps climbing to
170 MB/s — 1.4× faster than serial `gzip` *and* 1.4× smaller than `zstd -19`.
Full-range quality also makes each decoded byte more expensive than NovaSeq's
4-symbol binned quality.

### ONT MinION, ~14 kb reads (DRR205413) — 0.60 GB, 21,140 reads

| tool | size | ratio | 1 thr | 2 | 4 | 8 | 16 | 32 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `fqxv` (default) | 201 MB | 3.0× | 9.5 | 11 | 11 | 11 | 11 | 11 |
| `fqxv --max` | 197 MB | 3.1× | 9.5 | 11 | 11 | 11 | 11 | 11 |
| `gzip -dc` | 311 MB | 1.9× | 105 | — | — | — | — | — |
| `pigz -dc -p N` | 311 MB | 1.9× | 150 | 177 | 177 | 177 | 176 | 178 |
| `zstd -dc -T N` | 254 MB | 2.4× | 352 | 351 | 353 | 352 | 352 | 350 |

Here `fqxv` wins the ratio comparison outright (21% smaller than
`zstd -19 --long` already at the default level; `--max` adds 2% more) but
decode is flat at ~11 MB/s: the file holds only **2 blocks**, and long-read
decode does real work per base (the cross-read overlap codec reconstructs
each read from a consensus reference plus an edit script). Small long-read
files expose the serial cost; scaling returns with file size, not thread
count.

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
| DRR205413 | 2 | ~2 threads | 2 | ~2 threads |

The `--max` short-read path (read reordering) cuts its own smaller blocks
(~260 k reads), which is what un-sticks SRR2627175's 5-block ceiling. With
blocks to spare (DRR174812's 113), the `--max` curve still bends at ~16
threads: restoring the original read order from the stored permutation and
reassembling interleaved output is work the block fan-out cannot hide, and
the memory-bandwidth ceiling (below) applies sooner because each of its
threads decodes faster. The ONT tiling path keeps the same 2 blocks either
way.

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
bends (181 → 236 MB/s from 16 → 32 threads, not 2×) — more cores past that
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
- **Hardware/versions**: one Bodhi `rna`-partition node (Intel Xeon Gold
  6240R), dedicated 32-CPU Slurm allocations (not node-exclusive) — the main
  matrix and the `fqxv-max` cells ran as two back-to-back jobs pinned to the
  same node (64 GB and 128 GB; the `--max` reorder *compress* of the 10 GB
  input peaks at ~39 GB RSS — decode itself stays low). fqxv 0.6.2 (commit
  4940bf2), gzip 1.12, pigz 2.8, zstd 1.5.7, hyperfine 1.20.0.
