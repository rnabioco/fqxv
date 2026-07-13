# fqxv verify

Check the integrity of a `.fqxv` archive without decoding it or writing any
output. `verify` validates the CRC-32C checksums the container carries — the
whole-file checksum, the footer checksum, and every block frame — and reports
whether the archive is intact.

## Usage

```bash
fqxv verify <INPUT>
```

## Arguments

| Argument | Description |
| --- | --- |
| `<INPUT>` | Input `.fqxv` file. |

## Output and exit status

On success it prints to stdout and exits `0`:

```text
sample.fqxv: OK
```

On corruption it prints the reason to stderr and exits **non-zero** (`1`), so it
drops cleanly into scripts and CI:

```text
sample.fqxv: CORRUPT — <description of the failed check>
```

## Examples

```bash
# quick integrity check
fqxv verify sample.fqxv

# gate a pipeline on integrity (only decompress if the archive is sound)
fqxv verify sample.fqxv && fqxv decompress sample.fqxv -o sample.fastq

# verify a directory of archives, reporting any that fail
for f in *.fqxv; do fqxv verify "$f" || echo "FAILED: $f"; done
```

## Notes

- `verify` is a fast, one-pass integrity check: it reads the archive and
  recomputes checksums but never runs the codecs, so it is far cheaper than a
  full `decompress`.
- The checks are layered — a whole-file CRC-32C over the archive, a footer CRC
  guarding the row-group index, and a per-block payload CRC. See
  [Container Format](../design/container.md#integrity) for where each checksum
  lives on disk.
- To *recover data* from an archive that fails verification rather than just
  detect the damage, use [`decompress --recover`](decompress.md#recovering-a-corrupted-archive),
  which skips the corrupt blocks and decodes the rest.
