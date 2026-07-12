# fqxv-cli

The `fqxv` command-line tool: a reference-free FASTQ archiver for short-read
data. Installs the `fqxv` binary.

Part of the [`fqxv`](https://github.com/rnabioco/fqxv) workspace; see the
[documentation](https://rnabioco.github.io/fqxv/) for the full reference.

```bash
cargo install fqxv-cli

# single-end (gzip input auto-detected)
fqxv compress reads.fastq.gz -o reads.fqxv
fqxv decompress reads.fqxv -o reads.fastq

# paired-end / single-cell into one archive, split back or stream to an aligner
fqxv compress R1.fq R2.fq -o sample.fqxv
fqxv decompress sample.fqxv --split out
fqxv decompress sample.fqxv | bwa mem -p ref.fa -

fqxv info sample.fqxv
```

License: MIT OR Apache-2.0.
