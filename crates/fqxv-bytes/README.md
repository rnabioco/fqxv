# fqxv-bytes

Shared on-disk byte-serialization primitives for the `fqxv` codec crates:
unsigned LEB128 varints, zig-zag transforms, the read-length array codec
(`write_lens`/`read_lens`), and `Reader` — a bounds-checked cursor generic over
each crate's error type (via the `ReaderError` trait). These encodings were once
copy-pasted, byte-identical, into each codec; this leaf crate is the single
source of truth, so the on-disk byte layout can't drift between them.

Part of the [`fqxv`](https://github.com/rnabioco/fqxv) FASTQ-archiver workspace.

```rust
use fqxv_bytes::{write_varint, read_varint};

let mut buf = Vec::new();
write_varint(&mut buf, 300);
let mut pos = 0;
assert_eq!(read_varint(&buf, &mut pos), Some(300));
```

License: MIT OR Apache-2.0.
