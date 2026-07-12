# fqxv

The container format and compression library for the
[`fqxv`](https://github.com/rnabioco/fqxv) reference-free FASTQ archiver. It
composes the per-stream codecs — `fqxv-tokenizer` (names), `fqxv-seq`
(sequence), `fqxv-fqzcomp` (quality) — into a block-based, `rayon`-parallel
`.fqxv` container.

Supports single-end, paired-end, and N-way single-cell (R1/R2/I1/I2) inputs
(interleaved into one archive and split back out), interleaved streaming for
aligners, and an opt-in reorder mode.

```rust
use fqxv::{compress, decompress, Params};

let fastq = b"@r1\nACGT\n+\nIIII\n";
let mut archive = Vec::new();
compress(&fastq[..], &mut archive, Params::default()).unwrap();

let mut out = Vec::new();
decompress(&archive[..], &mut out, 1).unwrap();
assert_eq!(out, fastq); // name/seq/qual preserved; '+' normalized
```

For the command-line tool, see
[`fqxv-cli`](https://crates.io/crates/fqxv-cli).

License: MIT OR Apache-2.0.
