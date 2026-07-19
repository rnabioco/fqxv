# Python API

`fqxv` ships a **read-only** Python package for consuming `.fqxv` archives
directly from Python — no subprocess, no intermediate FASTQ file. Compression
stays in the [CLI](../cli/index.md); the Python side reads.

```bash
uv pip install fqxv
```

The wheels are `abi3` (one per platform, CPython ≥ 3.9) and carry the native
codecs, so there is no separate Rust toolchain to install. To build from a
checkout instead:

```bash
uv pip install maturin
maturin develop            # from crates/fqxv-python/
```

## Iterating records

`fqxv.open()` streams records in original file order. It works on **every**
archive layout, including globally-reordered (`--order shuffle`) archives, and
decodes on a background thread so memory stays bounded no matter the archive
size.

```python
import fqxv

for rec in fqxv.open("reads.fqxv"):
    print(rec.name, rec.sequence, rec.quality)   # all bytes
    print(len(rec))                              # sequence length
```

Each `Record` exposes `name` (the header with the leading `@` stripped),
`sequence`, and `quality` as `bytes`. Bytes — not `str` — because reads are
stored losslessly and read names are not guaranteed to be UTF-8.

In-memory input works anywhere a path does — pass `bytes` instead of a filename:

```python
data = open("reads.fqxv", "rb").read()
n = sum(1 for _ in fqxv.open(data))
```

`open()` also works as a context manager, and breaking out early is safe (the
decode thread is stopped cleanly):

```python
with fqxv.open("reads.fqxv") as reader:
    for rec in reader:
        if rec.sequence.startswith(b"ACGT"):
            break        # no hang; the decoder is torn down on exit
```

Pass `threads=` to control the decode pool (default `0` = a sensible default).

## Whole-archive convenience

```python
# Write interleaved FASTQ to a file; returns the read count.
n = fqxv.decompress_to_path("reads.fqxv", "reads.fastq")

# ...or get the FASTQ bytes directly.
raw = fqxv.decompress_to_bytes("reads.fqxv")

# Metadata only — no payload decode.
info = fqxv.inspect("reads.fqxv")
print(info.reads, info.blocks, info.format_version, info.platform)
```

`inspect()` returns an `Info` with `reads`, `blocks`, `group_size`,
`reordered`, `keep_order`, `regenerated_names`, `plus_normalized`,
`format_version`, `seq_order`, `quality_binning`, `names_bytes`,
`sequence_bytes`, `quality_bytes`, `platform`, and `whole_file_crc`.

## Column projection & random access

An `.fqxv` archive carries a footer index with per-stream byte offsets and CRCs.
That makes it **Parquet-shaped**: you can fetch and decode a
single column, or a single row group, without touching the rest. Read names are
typically under 1% of the archive, so an ID-only pass reads ~100× less.

```python
# Just the read names (IDs) — a fraction of the archive.
ids = fqxv.read_names("reads.fqxv")

# Just the sequences, skipping quality entirely.
seqs = fqxv.read_sequences("reads.fqxv")     # list[bytes], one per read

# Restrict to specific row groups.
first_group = fqxv.read_sequences("reads.fqxv", groups=[0])
```

`read_names`, `read_sequences`, and `read_qualities` each return a flat
`list[bytes]` across the requested groups (all groups when `groups` is omitted),
verifying every fetched stream against its stored CRC before decoding.

Inspect the index and decode a whole row group as records:

```python
idx = fqxv.open_index("reads.fqxv")
print(idx.total_reads, idx.num_groups)
for g in idx.groups():
    print(g.block_offset, g.read_count)

block0 = fqxv.read_block("reads.fqxv", 0)     # list[Record]
```

!!! note "Reordered archives"

    Projection and `open_index` require the plain (per-block) layout. A
    globally-reordered archive (`--order shuffle`) has no footer index — its
    streams are mutually dependent — so these raise `fqxv.FqxvError`. Use
    `fqxv.open()` to iterate those. Check `fqxv.inspect(path).reordered` if you
    need to branch. Compressing with a smaller `--block-reads` makes projection
    finer-grained on the plain layout.

## Errors

Decode and I/O failures raise exceptions: missing/unreadable files raise
`OSError`; a corrupt archive, an unsupported format version, or a projection on
a reordered archive raise `fqxv.FqxvError`.

```python
try:
    list(fqxv.open("truncated.fqxv"))
except fqxv.FqxvError as e:
    print("archive problem:", e)
```

## API reference

| Function | Returns | Notes |
| --- | --- | --- |
| `open(source, *, threads=0)` | `Reader` | Iterator of `Record`; every layout |
| `decompress_to_path(source, dest, *, threads=0)` | `int` | Read count; writes interleaved FASTQ |
| `decompress_to_bytes(source, *, threads=0)` | `bytes` | Interleaved FASTQ |
| `inspect(source)` | `Info` | Header + footer metadata |
| `open_index(source)` | `Index` | Footer row-group index (plain layout) |
| `read_names(source, groups=None)` | `list[bytes]` | Names for the groups (or all) |
| `read_sequences(source, groups=None)` | `list[bytes]` | Sequences for the groups (or all) |
| `read_qualities(source, groups=None)` | `list[bytes]` | Qualities for the groups (or all) |
| `read_block(source, group)` | `list[Record]` | Decode one whole row group |

`source` is a path (`str` / `os.PathLike`) or `bytes`. `Record` has
`.name` / `.sequence` / `.quality` (bytes) and `len()`; `Index` has
`.total_reads` / `.num_groups` / `.whole_file_crc` and `.groups()`.
