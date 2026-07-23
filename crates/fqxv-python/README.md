# fqxv (Python)

Read-only Python bindings for [`fqxv`](https://github.com/rnabioco/fqxv), a
Rust toolkit for lossless FASTQ archiving.

```python
import fqxv

# Stream records (works on every layout, including reordered archives)
for rec in fqxv.open("reads.fqxv"):
    print(rec.name, rec.sequence, rec.quality)  # all bytes

# In-memory input works too
data = open("reads.fqxv", "rb").read()
n = sum(1 for _ in fqxv.open(data))

# Whole-archive convenience
fqxv.decompress_to_path("reads.fqxv", "reads.fastq")
raw = fqxv.decompress_to_bytes("reads.fqxv")
info = fqxv.inspect("reads.fqxv")
print(info.reads, info.format_version, info.platform)

# Column projection / random access (plain layout only)
idx = fqxv.open_index("reads.fqxv")
seqs = fqxv.read_sequences("reads.fqxv")        # list[bytes], skips quality
ids = fqxv.read_names("reads.fqxv", groups=[0]) # just the first row group
block0 = fqxv.read_block("reads.fqxv", 0)       # list[Record]

# Read over the network (fqxv.remote, standard-library HTTP). Streaming decodes the
# whole archive on the fly; projection fetches only the column you ask for via HTTP
# byte-range requests.
import fqxv.remote as remote
for rec in remote.stream("https://host/reads.fqxv"):   # or a presigned S3 URL
    ...                                                # streams; no full download
arc = remote.open_index("https://host/reads.fqxv")     # 1 tail GET → footer index
names = arc.names()                                    # ~1% of the file, CRC-checked
print(arc.bytes_fetched, "of", arc.size)

# Any file-like works, so an AWS SDK response streams straight in — no fqxv HTTP:
import boto3
body = boto3.client("s3").get_object(Bucket=b, Key=k)["Body"]
for rec in fqxv.open(body):
    ...
# For concurrent async range fetches, drive fqxv.parse_index_suffix /
# Index.stream_range / fqxv.decode_*_bytes with your own httpx/aiohttp session
# (see fqxv.remote's module docstring).

# Integrity check — raises fqxv.FqxvError on a corrupt archive
fqxv.verify("reads.fqxv")

# Project the archive size/ratio from a FASTQ *without* compressing (gzip/BGZF ok)
est = fqxv.estimate("reads.fastq.gz", level=5)  # also: quality_binning=, sample_reads=
print(est.ratio, est.archive_bytes, est.exhausted)
```

Projection and `open_index` are unavailable for globally-reordered archives
(`--order shuffle`), whose streams are mutually dependent; use `fqxv.open()` to
iterate those. Everything here is read-only: `verify` and `estimate` only
*measure* — neither writes an archive — and full compression stays in the CLI.

## Build from source

```bash
uv pip install maturin
maturin develop            # from crates/fqxv-python/
```
