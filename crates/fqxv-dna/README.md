# fqxv-dna

Shared nucleotide primitives for the `fqxv` sequence codecs: the 2-bit `ACGT`
lookup (`BASE_LUT`, `SYM2BASE`, `code_strict` case-sensitive / `code_fold`
case-insensitive) and reverse complement (`revcomp`/`revcomp_into` complement
both cases; `revcomp_acgt` complements only uppercase and passes everything else
through). The 2-bit encodings and the reverse-complement mapping were once
copy-pasted — not always identically — into `fqxv-seq`, `fqxv-reorder`, and
`fqxv-lroverlap`; this leaf crate is the single source of truth so they can't
drift apart silently.

Part of the [`fqxv`](https://github.com/rnabioco/fqxv) FASTQ-archiver workspace.

```rust
use fqxv_dna::{code_strict, revcomp, NON_ACGT};

assert_eq!(code_strict(b'G'), 2);
assert_eq!(code_strict(b'N'), NON_ACGT); // non-ACGT sentinel
assert_eq!(revcomp(b"ACGTn"), b"nACGT");
```

License: MIT OR Apache-2.0.
