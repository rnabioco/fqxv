# fqxv-fqzcomp

A fqzcomp-style quality-score context model: each symbol is range-coded under a
context of the two previous qualities and position, one adaptive model per
context, reset at read boundaries. Opt-in lossy Illumina 2/4/8-level binning.

On real full-range Illumina quality it matches the C `fqz_comp`
(≈1.74 bits/byte). Part of the
[`fqxv`](https://github.com/rnabioco/fqxv) FASTQ-archiver workspace.

```rust
use fqxv_fqzcomp::{encode, decode, QualityBinning};

let lens = [5u32, 3];
let quals = b"IIIII##F"; // two reads
let enc = encode(&lens, quals, QualityBinning::Lossless).unwrap();
let (out_lens, out_quals) = decode(&enc).unwrap();
assert_eq!((out_lens, out_quals), (lens.to_vec(), quals.to_vec()));
```

License: MIT OR Apache-2.0.
