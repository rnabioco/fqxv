# fqxv-lroverlap

Cross-read overlap coding for long reads (Oxford Nanopore / PacBio). `fqxv-seq`
models each read independently with an order-k context model, which at long-read
coverage leaves the dominant redundancy on the table — the same locus is read
hundreds of times and each copy is coded from scratch. This crate closes that gap
by modeling reads against each other: minimizers → overlaps → layout → voted
consensus → a per-read banded edit script → rANS.

`encode`/`decode` are the container's sequence path for long-read blocks. The
container auto-selects the codec per block via the sequence stream's method byte
and keeps this one only when it beats the order-k model, so it can never enlarge
the archive. See [`docs/design/longread.md`](../../docs/design/longread.md) for
the measurements (parity with CoLoRd at depth on HiFi).

Part of the [`fqxv`](https://github.com/rnabioco/fqxv) FASTQ-archiver workspace.

```rust
use fqxv_lroverlap::{encode, decode, EncodeOpts};

let lens = [8u32, 8];
let seq = b"ACGTACGTACGTACGT";
let enc = encode(&lens, seq, &EncodeOpts::default()).unwrap();
let (out_lens, out_seq) = decode(&enc).unwrap();
assert_eq!((out_lens, out_seq), (lens.to_vec(), seq.to_vec()));
```

License: MIT OR Apache-2.0.
