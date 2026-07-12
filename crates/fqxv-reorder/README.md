# fqxv-reorder

Reference-free read reordering by canonical-minimizer clustering, plus a
byte-exact clustered differential codec (`encode_clustered` / `decode_clustered`)
that codes each read relative to the previous one (match / few-mismatch delta /
literal) — the SPRING/PgRC-style lever for cross-read duplicate redundancy.

Part of the [`fqxv`](https://github.com/rnabioco/fqxv) FASTQ-archiver workspace.

```rust
use fqxv_reorder::{plan, revcomp};

let a = b"ACGTTTGACCGATT";
let ra = revcomp(a);
let mut seq = a.to_vec();
seq.extend_from_slice(&ra);
let lens = [a.len() as u32, ra.len() as u32];
let p = plan(&lens, &seq, 7);
// the reverse-complement pair clusters; exactly one is flipped to match.
assert_ne!(p.flip[0], p.flip[1]);
```

License: MIT OR Apache-2.0.
