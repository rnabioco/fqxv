# fqxv-tokenizer

A positional read-name tokenizer for FASTQ headers. Names are split into
digit/non-digit runs, each token modeled against the previous record's token at
the same position (match / numeric delta / literal), and the role streams are
rANS-coded — so constant instrument prefixes and incrementing tile/x/y
coordinates collapse.

Part of the [`fqxv`](https://github.com/rnabioco/fqxv) FASTQ-archiver workspace.

```rust
use fqxv_tokenizer::{encode, decode};

let names: Vec<&[u8]> = vec![
    b"INST:1:FC:1:1101:1000:2000",
    b"INST:1:FC:1:1101:1005:2050",
];
let enc = encode(&names).unwrap();
assert_eq!(decode(&enc).unwrap(), vec![names[0].to_vec(), names[1].to_vec()]);
```

License: MIT OR Apache-2.0.
