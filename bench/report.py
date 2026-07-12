#!/usr/bin/env python3
"""Render the fqxv baseline benchmark as a per-dataset table.

Reads results.tsv + meta.tsv from $FQXV_RESULTS_DIR (default
$SCRATCH/fqxv/results) and prints, per dataset, each tool's compression ratio,
size relative to the gzip baseline, bits/base, throughput, and peak RSS.

The gzip row is the "baseline of record": `vs_gzip` is how many times smaller a
tool is than plain `.fastq.gz` (higher is better).
"""
from __future__ import annotations

import csv
import os
from collections import defaultdict
from pathlib import Path

RESULTS_DIR = Path(
    os.environ.get("FQXV_RESULTS_DIR")
    or f"{os.environ.get('SCRATCH', Path.home() / 'scratch')}/fqxv/results"
)


def load_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as fh:
        return list(csv.DictReader((r for r in fh if not r.startswith("#")), delimiter="\t"))


def main() -> None:
    results = load_tsv(RESULTS_DIR / "results.tsv")
    meta = {m["dataset"]: m for m in load_tsv(RESULTS_DIR / "meta.tsv")}

    by_ds: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in results:
        by_ds[row["dataset"]].append(row)

    hdr = f"{'tool':<10} {'comp_size':>10} {'ratio':>7} {'vs_gzip':>8} {'bits/base':>10} {'c_MB/s':>8} {'d_MB/s':>8} {'c_RSS':>8} {'rt':>4}"
    for ds, rows in by_ds.items():
        m = meta.get(ds, {})
        nbases = int(m.get("n_bases", 0) or 0)
        orig = int(m.get("orig_bytes", 0) or 0)
        gzip_bytes = next((int(r["comp_bytes"]) for r in rows if r["tool"] == "gzip"), 0)

        print(f"\n### {ds}  ({fmt_bytes(orig)}, {int(m.get('n_records', 0) or 0):,} reads, {fmt_bytes(nbases)} bases)")
        print(hdr)
        print("-" * len(hdr))
        for r in sorted(rows, key=lambda x: float(x["comp_bytes"]) or 1e18):
            cb = int(r["comp_bytes"])
            vs_gzip = gzip_bytes / cb if cb and gzip_bytes else 0.0
            bpb = (cb * 8) / nbases if nbases else 0.0
            c_mb = (orig / 1e6) / float(r["c_secs"]) if float(r["c_secs"] or 0) else 0.0
            d_mb = (orig / 1e6) / float(r["d_secs"]) if float(r["d_secs"] or 0) else 0.0
            rss = int(r["c_rss_kb"])
            print(
                f"{r['tool']:<10} {fmt_bytes(cb):>10} {r['ratio']:>7} "
                f"{vs_gzip:>7.2f}x {bpb:>10.3f} {c_mb:>8.1f} {d_mb:>8.1f} "
                f"{(fmt_bytes(rss*1024) if rss>=0 else 'n/a'):>8} {r['rt_ok']:>4}"
            )


def fmt_bytes(n: int) -> str:
    x = float(n)
    for unit in ("B", "K", "M", "G", "T"):
        if x < 1024 or unit == "T":
            return f"{x:.1f}{unit}" if unit != "B" else f"{int(x)}B"
        x /= 1024
    return f"{x:.1f}T"


if __name__ == "__main__":
    main()
