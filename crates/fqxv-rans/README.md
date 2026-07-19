# fqxv-rans

A clean-room [rANS](https://en.wikipedia.org/wiki/Asymmetric_numeral_systems)
Nx16 entropy coder — 32 interleaved states, 16-bit renormalization, order-0 and
order-1 models — implemented from the CRAM 3.1 codecs specification.

Order-0 encode **and** decode have scalar, AVX2, and AVX-512 backends chosen at
runtime; the widest path the CPU supports wins (AVX-512, else AVX2, else scalar).
Order-1 and every path below AVX2 run the scalar reference. There is no SSE4.2
vector path (gather requires AVX2), so `Backend::Sse42` runs scalar. **Every
backend produces byte-identical output**, so a stream encoded on one host decodes
bit-for-bit the same on any other.

Part of the [`fqxv`](https://github.com/rnabioco/fqxv) FASTQ-archiver workspace.

```rust
use fqxv_rans::{encode, decode, Order};

let data = b"the quick brown fox";
let compressed = encode(data, Order::One).unwrap();
assert_eq!(decode(&compressed).unwrap(), data);
```

`encode` and `decode` pick the backend automatically; `decode_with(_, backend)`
forces a specific `Backend` for testing or benchmarking.

License: MIT OR Apache-2.0.
