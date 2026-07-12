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

The gain scales with coverage depth (how many reads are matchable) but is large
even on shallow data. It is big enough to stay net-positive after storing the
permutation needed to restore the original read order.

Productionizing this as a container mode — dedup coding, context-model literals,
and a compact permutation for order-preserving output — is the path to the
SPRING/PgRC tier.
