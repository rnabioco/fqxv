# fqxv info

Print container metadata and per-stream compressed sizes for a `.fqxv` archive.

## Usage

```bash
fqxv info <INPUT> [--tsv | --json]
```

By default `info` prints a human-readable report of two tables (metadata and
per-stream sizes). Pass `--tsv` or `--json` for machine-readable output; the two
are mutually exclusive.

## Example

```bash
fqxv info sample.fqxv
```

```text
sample.fqxv
╭────────────────┬──────────────────────────────╮
│ property       │ value                        │
├────────────────┼──────────────────────────────┤
│ layout         │ paired (group size 2)        │
│ reads          │ 20,000                       │
│ spots          │ 10,000                       │
│ blocks         │ 1                            │
│ platform       │ Illumina                     │
│ sequence order │ 11                           │
│ quality        │ lossless                     │
│ reordered      │ no                           │
│ plus line      │ normalized                   │
│ file size      │ 1.01 MB (1,058,497 bytes)    │
╰────────────────┴──────────────────────────────╯
╭──────────┬───────────┬────────╮
│ stream   │     bytes │  share │
├──────────┼───────────┼────────┤
│ names    │    76,416 │   7.2% │
│ sequence │   455,180 │  43.0% │
│ quality  │   526,867 │  49.8% │
│ total    │ 1,058,463 │ 100.0% │
╰──────────┴───────────┴────────╯
52.92 bytes/read
```

## Machine-readable output

`--tsv` emits a fixed header line plus one data line. The columns are stable
(new fields are only ever appended), so it is safe to parse in scripts and the
benchmark harness:

```bash
fqxv info sample.fqxv --tsv
```

```text
file_size	reads	blocks	group_size	seq_order	quality_binning	reordered	names_bytes	seq_bytes	qual_bytes	platform
1058497	20000	1	2	11	0	0	76416	455180	526867	illumina
```

`--json` emits a single object with the same facts plus derived fields
(percentage shares, `bytes_per_read`, and human labels). `spots` and
`bytes_per_read` are omitted for single-end or empty archives:

```bash
fqxv info sample.fqxv --json
```

```json
{
  "file": "sample.fqxv",
  "file_size": 1058497,
  "file_size_human": "1.01 MB",
  "platform": "illumina",
  "reads": 20000,
  "spots": 10000,
  "blocks": 1,
  "layout": "paired",
  "group_size": 2,
  "sequence_order": 11,
  "quality": "lossless",
  "quality_binning": 0,
  "reordered": false,
  "read_order_preserved": true,
  "plus_normalized": true,
  "streams": {
    "names": { "bytes": 76416, "pct": 7.2 },
    "sequence": { "bytes": 455180, "pct": 43.0 },
    "quality": { "bytes": 526867, "pct": 49.8 },
    "total": { "bytes": 1058463, "pct": 100.0 }
  },
  "bytes_per_read": 52.92
}
```

## Fields

| Field | Meaning |
| --- | --- |
| `layout` | `single-end`, `paired`, or `grouped xG (single-cell)`. |
| `reads` / `spots` | Total reads; spots = reads / group size (paired/grouped only). |
| `blocks` | Number of independently-coded blocks. |
| `platform` | Sequencing platform (`Illumina`, `Oxford Nanopore`, `PacBio`, `MGI/BGI`, or `unknown`), detected from read names at compress time. Override with `compress --platform`. |
| `sequence order` | Order of the sequence context model. |
| `quality` | `lossless` or the lossy binning level used. |
| `reordered` | Whether reads were clustered/reordered for a better ratio. |
| `read order` | For reordered archives, whether the original order is restored. |
| `plus line` | Whether the `+` line was normalized. |
| `names` / `sequence` / `quality` | Compressed bytes per stream, with share of the three-stream total. |
