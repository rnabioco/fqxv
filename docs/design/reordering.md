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

## An honest finding

On real duplicate-rich RNA-seq (1M × 101 bp), reordering ahead of the order-11
`fqxv-seq` model gave:

| variant | bits/base |
| --- | --- |
| baseline (no reorder) | 1.247 |
| reordered (order not preserved) | 1.139 (**−8.7%**) |
| reordered + permutation (order preserved) | 1.770 (net **loss**) |

The reordering works, but the gain is modest and the permutation needed to
restore the original order costs more than it saves. The reason: `fqxv`'s
**high-order context model already captures most cross-read redundancy**, so
adding reordering double-counts it. SPRING and PgRC get large reordering gains
by pairing clustering with a *cheap* coder that can't see that redundancy on its
own.

So reordering ships as a library primitive rather than a default. The promising
directions are (1) an opt-in reorder-free mode for workflows that don't need the
original order, and (2) pairing the clustering with a lightweight
differential/assembly coder — future work.
