# fqxv info

Print container metadata and per-stream compressed sizes for a `.fqxv` archive.

## Usage

```bash
fqxv info <INPUT>
```

## Example

```bash
fqxv info sample.fqxv
```

```text
sample.fqxv
  layout         paired (group size 2)
  reads          20000
  spots          10000
  blocks         1
  sequence order 11
  quality        lossless
  plus line      normalized
  file size      1058497 bytes
  names     76416 bytes (7.2%)
  seq      455180 bytes (43.0%)
  qual     526867 bytes (49.8%)
  streams total 1058463 bytes (52.92 bytes/read)
```

## Fields

| Field | Meaning |
| --- | --- |
| `layout` | `single-end`, `paired`, or `grouped xG (single-cell)`. |
| `reads` / `spots` | Total reads; spots = reads / group size. |
| `blocks` | Number of independently-coded blocks. |
| `sequence order` | Order of the sequence context model. |
| `quality` | `lossless` or the lossy binning level used. |
| `plus line` | Whether the `+` line was normalized. |
| `names` / `seq` / `qual` | Compressed bytes per stream, with share of total. |
