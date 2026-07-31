# Read Reordering

`fqxv-reorder` clusters similar reads together so a downstream sequence coder
sees runs of near-identical reads — the cross-read redundancy lever that
per-read context modeling can't reach on its own (the idea behind SPRING and
PgRC).

## How it clusters

Each read is reduced to its **minimum canonical k-mer**: over every k-mer in the
read, take the smaller of the k-mer and its reverse complement, and keep the
minimum. Reads are then sorted by that key (and by oriented sequence as a
tie-break). Exact duplicates, reverse-complement duplicates, and near-duplicates
that share a minimizer all land next to each other. A per-read *flip* flag marks
reads stored reverse-complemented so a read and its RC-duplicate become
byte-identical after orientation.

`plan()` returns the emission order and the flip flags; the caller reorders the
name/sequence/quality streams accordingly and stores a permutation to restore
the original order.

## The pairing matters

Feeding reordered reads straight into the order-k context model barely helps —
the model already captures much of the redundancy, so clustering double-counts
it (a naïve reorder gave only ~9% on RNA-seq, negated by the permutation cost).

The win comes from **explicit differential coding** on the clustered reads —
SPRING's actual mechanism. Reads are assembled into contigs, each a growing
plurality-consensus reference, and every read is coded as one of three ops:

- **MATCH** — byte-identical to the previous read (nearly free),
- **CONTIG** — the read overlaps the current contig at the shift its shared
  minimizer implies: store the offset (delta-coded), the read length, the
  overlap's mismatch positions + substituted bases, and the novel tail, which
  then extends the reference,
- **LITERAL** — no placement: seeds a new contig, coded with the `fqxv-seq`
  context model.

Coding against the *consensus of every read placed so far* rather than against
the immediate predecessor is what captures the shifted overlaps of deep coverage.
Duplicates collapse to a single op; the unique reads still get the context model.
Each sequence block leads with a **version byte**: 2 for the single-contig codec
above, 3 for the literal-rescue variant that attaches reads v2 would strand. One or
the other is coded, not both — v3 generalizes v2 and otherwise degenerates to the
same coding, so it is the block-local floor on its own, and `--no-rescue` picks the
faster v2 path.

**Global reference.** Beyond the block-local codecs, the container can also
assemble **one frozen whole-file reference** (SPRING-style) over every clustered
read and code reads as `(contig, offset, mismatches)` positions on it — a
**version-4** sequence block — so the cross-block overlaps the block-local codecs
strand as literals collapse to a cheap back-reference. The assembly runs the greedy
fold over a fixed set of windows **in parallel** (a fixed count, never derived from
`--threads`, so the reference is byte-identical regardless of thread count); the
reference is then **overlap-merged** — contigs whose suffix overlaps another's
prefix are chained into fewer, longer super-contigs, reclaiming the cross-window
deduplication the split gave up — and stored once.

The reference frame itself is coded two ways and the smaller kept: SPRING's
**2-bit-pack + LZMA** on the packed consensus (usually the winner — the packing is
a hard 2 bits/base floor and the byte-domain LZ then reaches the long-range
near-duplicate-contig repeats a context model cannot see) and a **block-parallel
order-k `fqxv-seq`** pass over a fixed 64 contig blocks. Both are in-tree and
clean-room: there is no external/C compressor in this path — the earlier xz
(`liblzma`) reference coder was removed, and the LZMA here is `fqxv-seq`'s own.

This is the adaptive `rescue` path (default under `--order any`; `--no-rescue`
turns it off). It is adopted only when the reference frame plus the per-block
winners (v4 where it is smaller, the block-local codec otherwise) beat the
block-local total, so it can only ever shrink the archive. See [Container
Format → Reordered archives](container.md#reordered-archives) for the on-disk frame
layout and the reference-frame method byte.

Measured on the sequence stream:

| dataset | `fqxv-seq` order-11 | reorder + delta + ctx-literals | gain |
| --- | --- | --- | --- |
| E. coli, ~119× coverage | 1.344 bits/base | **0.737** | **−45%** |
| RNA-seq, shallow | 1.247 bits/base | **0.949** | **−24%** |

These are *idealized* numbers: fixed read length, sequence stream only, order
not preserved.

## End to end, and the real-world caveats

`fqxv compress --order any` (or `--max`) turns on read reordering. On single-end
input the reads emerge clustered (order not preserved). On grouped (paired /
single-cell) input a permutation is stored so the mate interleaving is
reconstructed on decompress and `--split` — grouped reorder is therefore always
order-preserving. For single-end input, whether the original order is kept is
picked adaptively (a stored permutation wins for counter-style names that
delta-code to almost nothing in original order); the Advanced `--keep-order` flag
forces it on, and `--order shuffle` opts into discarding order entirely by
regenerating purely positional names from a template (reorder-lossy — reads are
renumbered, sequence and quality preserved exactly). All modes round-trip
exactly (or exactly as a set for single-end `any`). On a **full, real** deep
dataset (E. coli, 2.19 M variable-length reads) the whole-archive gain is modest:

| mode | size | vs plain |
| --- | --- | --- |
| plain (`--order preserve`) | 255.8 MB | — |
| `--order any` (order not kept) | 247.4 MB | −3.3% |
| reorder + stored permutation | 253.9 MB | −0.7% |

This table predates the whole-file global reference and the overlap-merge that
followed it; for current whole-archive numbers see
[Benchmarks](../benchmarks.md).

Three things erode the idealized gain on real data:

1. **Variable read lengths.** `MATCH` needs byte identity, so trimmed reads
   (249/250/251 bp) do not collapse to a single op — each pays its own offset,
   length, mismatch set, and novel tail as a `CONTIG` placement, or falls to
   `LITERAL` when it will not place at all. The 45% was measured on a fixed-251bp
   subset.
2. **Reordering scrambles read names**, which destroys the tokenizer's
   match/delta structure — a cost that partly offsets the sequence gain.
3. **The permutation** (order-preserving / grouped reorder) is expensive at scale.

Two of those three have since been addressed: contig placement removed the
equal-length requirement, and the whole-file global reference collects the
cross-block overlaps the block-local codecs stranded. Name scrambling is the one
that stands — which is why `--order shuffle`, the mode that stops paying for names
at all by regenerating them from a counter template, is the mode that reaches the
SPRING/PgRC tier (and, on both benchmark datasets, past SPRING; see
[Benchmarks](../benchmarks.md)). The remaining lever is the permutation itself:
keeping the original order costs bytes that no reordering gain currently offsets.
