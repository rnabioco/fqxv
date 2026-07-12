# fqxv-range

A Subbotin carryless range coder plus an adaptive frequency model
(`SimpleModel<N>`) — the serial entropy backend for context models such as the
[`fqxv-fqzcomp`](https://crates.io/crates/fqxv-fqzcomp) quality coder.

Part of the [`fqxv`](https://github.com/rnabioco/fqxv) FASTQ-archiver workspace.

```rust
use fqxv_range::{Encoder, Decoder, SimpleModel};

let data = [3usize, 1, 4, 1, 5, 9, 2, 6];
let mut enc = Encoder::new();
let mut m = SimpleModel::<10>::new();
for &s in &data { m.encode(&mut enc, s); }
let bytes = enc.finish();

let mut dec = Decoder::new(&bytes);
let mut m = SimpleModel::<10>::new();
let out: Vec<usize> = (0..data.len()).map(|_| m.decode(&mut dec)).collect();
assert_eq!(out, data);
```

License: MIT OR Apache-2.0.
