# fqxv verify

Check the integrity of a `.fqxv` archive without decoding it or writing any
output. `verify` validates the CRC-32C checksums the container carries — the
whole-file checksum, the footer checksum, and every block frame — and reports
each check in a table.

## Usage

```bash
fqxv verify [--quick] [--tsv | --json] <INPUT>
```

## Arguments

| Argument | Description |
| --- | --- |
| `<INPUT>` | Input `.fqxv` file. |

## Options

| Option | Description |
| --- | --- |
| `--quick` | Faster, weaker check: validate each block's stored CRC via the footer index (parallel positioned reads) instead of the whole-file digest. |
| `--tsv` | Emit tab-separated per-check rows instead of the table. |
| `--json` | Emit a JSON object instead of the table. |
| `--threads <N>` | Worker threads (0 = all cores) for the parallelized CRC pass. |

## Output and exit status

The default output is a table of per-check results followed by an overall
verdict line. It exits `0` when intact and **non-zero (`1`) when corrupt** in
every output format, so it drops cleanly into scripts and CI.

```text
sample.fqxv
╭────────────────┬────────┬───────────────────────────╮
│ check          │ result │ detail                    │
├────────────────┼────────┼───────────────────────────┤
│ header         │ ok     │ format v1.0, plain layout │
│ footer         │ ok     │ 128 blocks, 33200000 …    │
│ block CRCs     │ ok     │ 128/128 intact            │
│ whole-file CRC │ ok     │                           │
╰────────────────┴────────┴───────────────────────────╯
sample.fqxv: OK
```

When a check fails, its row shows `FAIL` and the reason. A whole-file CRC failure
triggers a block-by-block scan that names the damaged blocks:

```text
│ block CRCs     │ FAIL   │ 126/128 intact; failed: 5, 91 │
│ whole-file CRC │ FAIL   │ digest mismatch               │
...
sample.fqxv: CORRUPT
```

`--quick` runs only the per-block payload CRCs, so its table omits the
`whole-file CRC` row (see the note below on what it trades away).

### Machine-readable output

`--tsv` emits a header line plus one row per check (`result` is `ok`/`fail`):

```bash
fqxv verify sample.fqxv --tsv
```

```text
check	result	detail
header	ok	format v1.0, plain layout
footer	ok	128 blocks, 33200000 reads
block CRCs	ok	128/128 intact
whole-file CRC	ok	
```

`--json` emits an object with the overall `passed` verdict, the per-check array,
and any `failed_blocks` (`--tsv` and `--json` are mutually exclusive):

```json
{
  "file": "sample.fqxv",
  "passed": true,
  "checks": [
    { "name": "header", "ok": true, "detail": "format v1.0, plain layout" },
    { "name": "footer", "ok": true, "detail": "128 blocks, 33200000 reads" },
    { "name": "block CRCs", "ok": true, "detail": "128/128 intact" },
    { "name": "whole-file CRC", "ok": true, "detail": "" }
  ]
}
```

## Examples

```bash
# quick integrity check
fqxv verify sample.fqxv

# gate a pipeline on integrity (only decompress if the archive is sound)
fqxv verify sample.fqxv >/dev/null && fqxv decompress sample.fqxv -o sample.fastq

# verify a directory of archives, reporting any that fail
for f in *.fqxv; do fqxv verify "$f" >/dev/null || echo "FAILED: $f"; done
```

## Notes

- `verify` never runs the codecs — the default check reads the archive once and
  recomputes the whole-file CRC-32C (parallelized across cores), which covers
  **every** byte: the header, the footer index, the inter-block framing, and
  every payload. Only when that CRC fails does it make a second pass to localize
  the damaged blocks.
- `--quick` instead checks only the per-block payload CRCs, reading blocks
  concurrently from their footer offsets. On a large multi-block archive this is
  typically a few times faster and still localizes any failure to specific
  blocks, but it is a **weaker** check: corruption confined to the header,
  footer, or framing bytes is not detected, and a single-block archive gains no
  parallelism. The globally-reordered layout has no per-block index, so both
  modes verify it by decoding into a sink (one `streams (decode)` check).
- To *recover data* from an archive that fails verification rather than just
  detect the damage, use [`decompress --recover`](decompress.md#recovering-a-corrupted-archive),
  which skips the corrupt blocks and decodes the rest.
