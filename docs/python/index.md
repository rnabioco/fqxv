# Python API

`fqxv` ships a **read-only** Python package for consuming `.fqxv` archives
directly from Python — no subprocess, no intermediate FASTQ file. Compression
stays in the [CLI](../cli/index.md); the Python side reads, and with `estimate()`
measures without writing anything.

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

`fqxv.open()` streams records in the order the archive stores them — the original
file order for every order-preserving archive, which is the default and includes
every paired archive. It works on **every** layout, including the
globally-reordered ones (`--order any`, `--max`, `--order shuffle`), where
single-end reads come back clustered rather than in file order; check
`fqxv.inspect(path).keep_order` if that matters. Decoding runs on a background
thread behind a bounded channel, so a plain archive streams in constant memory;
a globally-reordered archive has mutually dependent streams and is decoded
whole-file, so it stays resident while you iterate.

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

So does any file-like object with `.read()` — an HTTP response or an AWS SDK
body streams straight in, with no staging file and no HTTP code in fqxv:

```python
import boto3
body = boto3.client("s3").get_object(Bucket=bucket, Key=key)["Body"]
for rec in fqxv.open(body):
    ...
```

File-like sources work for the forward-only entry points — `open()`,
`decompress_to_path()`, `decompress_to_bytes()`. Everything else (`inspect`,
`open_index`, the projections, `read_block`, `verify`) seeks, so it takes a path
or `bytes`.

`open()` also works as a context manager, and breaking out early is safe (the
decode thread is stopped cleanly):

```python
with fqxv.open("reads.fqxv") as reader:
    for rec in reader:
        if rec.sequence.startswith(b"ACGT"):
            break        # no hang; the decoder is torn down on exit
```

Pass `threads=` to control the decode pool (default `0` = all cores).

## Decoding only some streams

Many consumers never read all three streams — k-mer counters and classifiers
want only sequences, an ID audit only names. `open()` takes a `streams=`
selection: one stream name or an iterable of names out of `"names"`,
`"sequence"`, `"quality"` (the short forms `"name"`/`"seq"`/`"qual"` also
work). Deselected fields come back as empty `bytes`, and their coded streams
are **skipped, never entropy-decoded** — the decode cost of everything you did
not ask for disappears.

```python
# Sequence-only pass: k-mer counting / classification input, straight from the archive.
for rec in fqxv.open("reads.fqxv", streams=("seq",)):
    count(rec.sequence)            # rec.name == rec.quality == b""

ids = [r.name for r in fqxv.open("reads.fqxv", streams="names")]

# streams=() decodes nothing: count reads at framing speed.
n = sum(1 for _ in fqxv.open("reads.fqxv", streams=()))
```

How much skipping buys depends on where the archive's decode compute sits:
sequence-only decode measured **~12× faster** on an ONT long-read archive (the
sequence-conditioned quality model dominates decode there) and ~1.3–1.6×
single-threaded on short-read Illumina; at high thread counts the short-read
win shows up as freed CPU time rather than wall time. Two caveats:

- **Long-read quality needs the sequence.** On long-read archives quality is
  coded against the bases, so selecting quality still decodes the sequence
  stream internally (it is returned empty unless selected). Names and sequence
  have no such dependency.
- **Skipped streams are not verified.** Every stream that *is* decoded keeps
  its integrity checks, but a skipped stream's content digest cannot be checked
  (its content is never reconstructed). Run `fqxv.verify()` when you need the
  whole archive vouched for.

Selection works on **every** layout, including globally-reordered archives —
unlike the [column projections](#column-projection-random-access) below, which
need the plain layout's footer index but can also skip *reading* the deselected
bytes (the streaming selection reads the whole file and skips the decode).

### FASTA output

`decompress_to_path()` and `decompress_to_bytes()` take `fasta=True` — the CLI's
`decompress --fasta` as an API: single-line FASTA (`>name` + sequence), with the
quality stream skipped entirely, at the same speedups as above.

```python
fqxv.decompress_to_path("reads.fqxv", "reads.fasta", fasta=True)
fa = fqxv.decompress_to_bytes("reads.fqxv", fasta=True)
```

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
`sequence_bytes`, `quality_bytes`, `platform`, `whole_file_crc`,
`required_features`, and `member_labels`.

`format_version` is the **container format** version packed as
`(major << 8) | minor` — `256` for format 1.0 — and is independent of the
`fqxv` package version. `member_labels` holds the original per-slot labels
(`["R1", "R2"]`, `["2", "4"]`, …) when the archive recorded them, and is empty
otherwise. `platform` is the human label (`"Illumina"`, `"Oxford Nanopore"`,
`"PacBio"`, `"MGI/BGI"`, `"unknown"`); `seq_order`, `quality_binning`, and
`required_features` are the raw header codes (`quality_binning` `0` = lossless);
`whole_file_crc` is `None` for the footer-less globally-reordered layout.

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
verifying every fetched stream against its stored CRC before decoding. Long-read
quality is coded against the sequence, so `read_qualities` also reads the
sequence column for those groups; names and sequences stand alone.

Inspect the index and decode a whole row group as records:

```python
idx = fqxv.open_index("reads.fqxv")
print(idx.total_reads, idx.num_groups)
for g in idx.groups():
    print(g.block_offset, g.read_count)

block0 = fqxv.read_block("reads.fqxv", 0)     # list[Record]
```

`groups()` returns `GroupLoc` objects with `block_offset` and `read_count`.
`Index` also carries the per-stream primitives a remote client drives:
`stream_range(group, stream)` → `(start, end)` (half-open),
`stream_crc(group, stream)`, and `verify_stream(group, stream, coded)`, where
`stream` is `"names"`, `"sequence"`, or `"quality"`.

!!! note "Reordered archives"

    Projection and `open_index` require the plain (per-block) layout. A
    globally-reordered archive (`--order any`, `--max`, `--order shuffle`) has no
    footer index — its streams are mutually dependent — so these raise
    `fqxv.FqxvError`. Use `fqxv.open()` to iterate those (with `streams=` if you
    only need a subset). Check
    `fqxv.inspect(path).reordered` if you need to branch. Compressing with a
    smaller `--block-reads` makes projection finer-grained on the plain layout.

## Estimating a compression ratio

`fqxv.estimate()` projects an archive's size from a **FASTQ** input (not an
archive) by coding a bounded leading sample with the real codecs — the library
behind `fqxv compress --estimate`. gzip/BGZF input is decoded transparently, and
nothing is written.

```python
est = fqxv.estimate("reads.fastq.gz")            # level=5 by default
print(est.ratio, est.archive_bytes, est.exhausted)

# Mates that would compress into one archive: pass a list (or tuple).
est = fqxv.estimate(["R1.fastq.gz", "R2.fastq.gz"], level=9)
```

`Estimate` carries `sample_reads`, `sample_bases`, `raw_bytes`, `names_bytes`,
`sequence_bytes`, `quality_bytes`, `archive_bytes`, `exhausted`, `platform`, and
`ratio`. The byte counts describe the sample; `ratio` (raw ÷ archive) is
scale-invariant and holds for the whole file. `exhausted` is `True` when the whole
input fit in the sample, so the numbers are exact rather than an extrapolation
base. `platform` is the lowercase token detected from the sample (`"illumina"`,
`"nanopore"`, `"pacbio"`, `"mgi"`, `"unknown"`) — not the same spelling as
`Info.platform`, which is the human label.

Keyword arguments: `level` (1-9, default `5`), `quality_binning` (`"lossless"`,
`"bin8"`, `"bin4"`, `"bin2"`, `"ont"`, `"hifi"`), `sample_reads` (default
`1048576`), and `threads`. A grouped estimate splits `sample_reads` across its
inputs and rejects a cross-platform mix. Reordering and interleaving are not
modelled, so the real archive comes out this size or smaller.

## Verifying an archive

```python
fqxv.verify("reads.fqxv")     # returns None; raises if the archive is bad
```

Checks the header, footer, and stored whole-file CRC in one linear pass (the
globally-reordered layout, which has no footer, is checked by a full decode).
This is the tool for "is this archive intact": a *streaming* read of a truncated
archive can end early and silently, since the block region carries no running
read count.

## Reading over HTTP

`fqxv.remote` reads an archive over HTTP using the standard library only (no
dependencies). Streaming feeds the response straight into the decoder; projection
fetches the footer index from the archive tail and then only the columns you ask
for, with `Range` requests.

```python
import fqxv.remote as remote

for rec in remote.stream("https://host/reads.fqxv"):   # or a presigned S3 URL
    ...                                                # streams; no full download

arc = remote.open_index("https://host/reads.fqxv")     # one tail GET
names = arc.names()                                    # ~1% of the file, CRC-checked
print(arc.bytes_fetched, "of", arc.size)

n = remote.download("https://host/reads.fqxv", "reads.fastq")   # read count
```

`RemoteArchive` exposes `.names()`, `.sequences()`, `.qualities()`, and
`.records()` — each taking an optional `groups=` — plus `.index`, `.size`, and
`.bytes_fetched`. The module-level shortcuts are `open_index`, `read_names`,
`read_sequences`, `read_qualities`, `stream`, and `download`; all accept
`headers=` for an `Authorization` header on a private object. `stream()` also
takes [`streams=`](#decoding-only-some-streams) (decode only some streams while
streaming — the archive body is still transferred) and `download()` takes
`fasta=True`.

To drive the same thing from a different client — an async `httpx`/`aiohttp`
session issuing range fetches concurrently — call the IO-free primitives
directly:

```python
index, need = fqxv.parse_index_suffix(tail, file_size)   # tail = bytes=-65536 body
start, end = index.stream_range(0, "names")              # then GET bytes=start-(end-1)
index.verify_stream(0, "names", coded)
names = fqxv.decode_names_bytes(coded)
```

`parse_index_suffix(suffix, file_len)` returns `(Index, None)` when the fetched
tail reached the footer, or `(None, need_at_least)` when it fell short — refetch
that many trailing bytes and call again; the second length is exact, so it costs
at most one extra round trip. Long-read quality is coded against the sequence:
`quality_needs_sequence_bytes(coded)` reports that, and
`decode_qualities_bytes(coded, seq)` takes the group's concatenated decoded bases.

## Errors

Decode and I/O failures raise exceptions: a missing or unreadable file raises
`OSError`; a corrupt or truncated archive, an unsupported format version or
feature, and any projection on a reordered archive raise `fqxv.FqxvError`. A
`source` of the wrong type raises `TypeError`, and a bad argument value (an
unknown `quality_binning` or stream name, an empty source list) raises
`ValueError`.

```python
try:
    fqxv.verify("suspect.fqxv")
except fqxv.FqxvError as e:
    print("archive problem:", e)
except OSError as e:
    print("cannot read it:", e)
```

## API reference

| Function | Returns | Notes |
| --- | --- | --- |
| `open(source, *, threads=0, streams=None)` | `Reader` | Iterator of `Record`; every layout. `streams=` decodes a subset — deselected fields are `b""`, their streams skipped |
| `decompress_to_path(source, dest, *, threads=0, fasta=False)` | `int` | Read count; interleaved FASTQ, or single-line FASTA with `fasta=True` (quality skipped) |
| `decompress_to_bytes(source, *, threads=0, fasta=False)` | `bytes` | Interleaved FASTQ, or single-line FASTA with `fasta=True` |
| `inspect(source)` | `Info` | Header + footer metadata |
| `open_index(source)` | `Index` | Footer row-group index (plain layout) |
| `read_names(source, groups=None)` | `list[bytes]` | Names for the groups (or all) |
| `read_sequences(source, groups=None)` | `list[bytes]` | Sequences for the groups (or all) |
| `read_qualities(source, groups=None)` | `list[bytes]` | Qualities for the groups (or all) |
| `read_block(source, group)` | `list[Record]` | Decode one whole row group |
| `estimate(source, *, level=5, quality_binning="lossless", sample_reads=1048576, threads=0)` | `Estimate` | Projected size/ratio from a FASTQ; `source` may be a list of mates |
| `verify(source, *, threads=0)` | `None` | Integrity check; raises if the archive is bad |

Random-access primitives — no I/O of their own, for driving a custom (e.g. async)
remote client; `fqxv.remote` is built on them:

| Function | Returns | Notes |
| --- | --- | --- |
| `parse_index_suffix(suffix, file_len)` | `(Index \| None, int \| None)` | Footer index from a fetched tail, else the bytes needed |
| `decode_names_bytes(coded)` | `list[bytes]` | Decode one fetched names column |
| `decode_sequences_bytes(coded)` | `list[bytes]` | Decode one fetched sequence column |
| `decode_qualities_bytes(coded, seq=None)` | `list[bytes]` | Decode one fetched quality column |
| `quality_needs_sequence_bytes(coded)` | `bool` | Whether that column needs `seq` |

`source` is a path (`str` / `os.PathLike`) or `bytes`; `open`,
`decompress_to_path`, and `decompress_to_bytes` also accept a file-like object
with `.read()`. `Record` has `.name` / `.sequence` / `.quality` (bytes) and
`len()`; `Index` has `.total_reads` / `.num_groups` / `.whole_file_crc`,
`.groups()`, and `.stream_range()` / `.stream_crc()` / `.verify_stream()`;
`GroupLoc` has `.block_offset` / `.read_count`.
