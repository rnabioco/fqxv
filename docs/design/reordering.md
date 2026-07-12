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
SPRING's actual mechanism. Each read is coded relative to the previous read in
the reordered stream:

- **MATCH** — identical to the previous read (nearly free),
- **DELTA** — same length, a few mismatches: store the mismatch positions + bases,
- **LITERAL** — everything else, coded with the `fqxv-seq` context model.

Duplicates collapse to a single op; the unique reads still get the context
model. Measured on the sequence stream:

| dataset | `fqxv-seq` order-11 | reorder + delta + ctx-literals | gain |
| --- | --- | --- | --- |
| E. coli, ~119× coverage | 1.344 bits/base | **0.737** | **−45%** |
| RNA-seq, shallow | 1.247 bits/base | **0.949** | **−24%** |

These are *idealized* numbers: fixed read length, sequence stream only, order
not preserved.

## End to end, and the real-world caveats

`fqxv compress --reorder` ships as a single-end container mode (`--keep-order`
stores a permutation for byte-exact restore; without it, reads emerge clustered).
It round-trips exactly. But on a **full, real** deep dataset (E. coli, 2.19 M
variable-length reads) the whole-archive gain is modest:

| mode | size | vs plain |
| --- | --- | --- |
| plain | 255.8 MB | — |
| `--reorder` (order not kept) | 247.4 MB | −3.3% |
| `--reorder --keep-order` | 253.9 MB | −0.7% |

Three things erode the idealized gain on real data:

1. **Variable read lengths.** `MATCH`/`DELTA` require equal length, so trimmed
   reads (249/250/251 bp) mostly fall to `LITERAL`. The 45% was measured on a
   fixed-251bp subset.
2. **Reordering scrambles read names**, which destroys the tokenizer's
   match/delta structure — a cost that partly offsets the sequence gain.
3. **The permutation** for `--keep-order` is expensive at scale.

So the honest state: the mechanism is real and validated, but realizing it on
everyday data needs variable-length-aware differential coding (align, or allow
length-changing deltas) and a way to keep names in original order while the
sequence is reordered. That's the remaining work toward the SPRING/PgRC tier.
