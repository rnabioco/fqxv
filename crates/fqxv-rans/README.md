# fqxv-rans

A clean-room [rANS](https://en.wikipedia.org/wiki/Asymmetric_numeral_systems)
Nx16 entropy coder — 32 interleaved states, 16-bit renormalization, order-0 and
order-1 models — with scalar and AVX2 decode backends selected at runtime.
Implemented from the CRAM 3.1 codecs specification.

Part of the [`fqxv`](https://github.com/rnabioco/fqxv) FASTQ-archiver workspace.

```rust
use fqxv_rans::{encode, decode, Order};

let data = b"the quick brown fox";
let compressed = encode(data, Order::One).unwrap();
assert_eq!(decode(&compressed).unwrap(), data);
```

> Note: on AMD Zen 3 the AVX2 gather path measured slower than the autovectorized
> scalar path, so `decode` defaults to scalar; select AVX2 explicitly with
> `decode_with(_, Backend::Avx2)`.

License: MIT OR Apache-2.0.
