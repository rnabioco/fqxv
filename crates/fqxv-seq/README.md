# fqxv-seq

Nucleotide sequence coding via an order-k adaptive context model over a 2-bit
A/C/G/T alphabet, range-coded. Non-ACGT bytes (`N`, IUPAC codes, lowercase) go
to a delta-coded exception list, so it is byte-exact.

Part of the [`fqxv`](https://github.com/rnabioco/fqxv) FASTQ-archiver workspace.

```rust
use fqxv_seq::{encode, decode};

let lens = [4u32, 5];
let seq = b"ACGTACGTN";
let enc = encode(&lens, seq, 6).unwrap();
let (out_lens, out_seq) = decode(&enc).unwrap();
assert_eq!((out_lens, out_seq), (lens.to_vec(), seq.to_vec()));
```

License: MIT OR Apache-2.0.
