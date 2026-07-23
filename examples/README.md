# Examples

## `remote_smoke.sh` — live remote-read validation

An opt-in, end-to-end smoke test of reading a `.fqxv` over the network. It is a
*network* test against a real Range-capable host, so it is **not** wired into
`cargo nextest` / CI — run it by hand against a URL you control. It exercises both
halves of the remote-read feature and asserts five things:

1. the server supports HTTP `Range` (a `206` byte-range probe returns the `FQXV` magic),
2. streaming decode over the pipe succeeds and yields reads (`curl | fqxv decompress -`),
3. the streamed decode is **byte-identical** to a downloaded-file decode,
4. a **truncated** transfer is rejected (premature EOF before the terminator frame), and
5. the Python `fqxv.remote.RemoteArchive` column projection reads the names over
   `Range` requests, matches the full decode, and transfers a fraction of the file.

### Run it against any URL

```bash
# defaults to a shared fixture; override with an argument or $FQXV_REMOTE_URL
examples/remote_smoke.sh https://my-host/reads.fqxv

# point it at a specific fqxv binary / Python build
FQXV_BIN=./target/release/fqxv examples/remote_smoke.sh https://my-host/reads.fqxv
```

The script resolves the CLI via `FQXV_BIN`, then `PATH`, then `target/{release,debug}`,
and adds the in-repo Python package to `PYTHONPATH` if it was built with
`maturin develop` (the `_fqxv` extension is present). Step 5 prints `SKIP` if the
`fqxv` Python module isn't importable.

---

## Validating against AWS S3

The remote-read design is a **stdin pipe** — the transfer tool owns the download and
fqxv just decodes the stream, so there is no HTTP/AWS code in fqxv itself. Any tool
that can write the object to stdout works (`aws`, `gsutil`, `curl` on a presigned URL).

**Prerequisites:** the AWS CLI v2 and working credentials. Verify with:

```bash
aws sts get-caller-identity     # should print your account / user ARN
```

### Setup — a small test bucket

Bucket names are globally unique; leave the defaults on (Block Public Access stays
*enabled* — access is via your credentials or a presigned URL, never public). Pick a
nearby region (`us-west-2` shown):

```bash
aws s3 mb s3://YOURNAME-fqxv-test --region us-west-2

# upload some archives (regenerate old ones first — see note below)
aws s3 cp ./my-archives/ s3://YOURNAME-fqxv-test/ \
  --recursive --exclude "*" --include "*.fqxv"

aws s3 ls s3://YOURNAME-fqxv-test/
```

### Read — two ways

```bash
# A) native pipe — uses your AWS credential chain, no presigning, no public access
aws s3 cp s3://YOURNAME-fqxv-test/reads.fqxv - | fqxv decompress - -Z | bwa mem -p ref.fa -

# B) presigned https URL — needed for curl and for the Python range projection,
#    which uses urllib and cannot resolve an s3:// scheme
URL=$(aws s3 presign s3://YOURNAME-fqxv-test/reads.fqxv --expires-in 3600)
examples/remote_smoke.sh "$URL"
```

### Teardown

```bash
aws s3 rb s3://YOURNAME-fqxv-test --force     # deletes all objects + the bucket
```

### Notes

- **Cost is negligible** for testing — S3 Standard storage is ~$0.023/GB·month
  (≈$0.002/mo for 100 MB), requests are fractions of a cent, and egress has a
  100 GB/month free allowance. Free-Tier accounts pay nothing.
- **Use S3 Standard**, not Intelligent-Tiering / IA / Glacier (their
  minimum-duration and retrieval fees make transient test data *more* expensive),
  and **skip CloudFront** — a CDN/proxy in front is the one place `Range: bytes=-N`
  (suffix ranges, which the projection's tail-index fetch relies on) can be
  mishandled. Plain S3 honors suffix ranges.
- **Presigned URLs and `Range`:** SigV4 does not sign the `Range` header, so
  suffix/byte-range requests pass through a presigned URL fine.
- **Regenerate stale archives** before uploading: an archive written by an older
  build may lack the footer index the projection needs — a good one ends with the
  `FQXF` trailer (`tail -c 4 file.fqxv | xxd` → `4651 5846`). A globally-reordered
  (`--order shuffle`) archive has no index and supports streaming only, not projection.
- **Credentials on a shared machine:** `chmod 600 ~/.aws/credentials` and remove them
  when finished; long-lived keys shouldn't linger on shared storage.

GCS is the same shape (`gsutil cat gs://… | fqxv decompress -`); presign with
`gcloud storage sign-url` for the projection path.
