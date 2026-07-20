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
import math
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

    hdr = (
        f"{'tool':<13} {'comp_size':>10} {'ratio':>7} {'vs_gzip':>8} "
        f"{'bits/base':>10} {'c_MB/s':>8} {'d_MB/s':>8} {'c_RSS':>8} "
        f"{'rt':>4} {'det':>4}"
    )
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
            # Guard on > 0, not merely non-zero: -1 is the "not measured"
            # sentinel (a tool whose binary is absent, rt=miss), and dividing by
            # it printed a negative throughput as if it were a measurement.
            c_mb = (orig / 1e6) / float(r["c_secs"]) if float(r["c_secs"] or 0) > 0 else 0.0
            d_mb = (orig / 1e6) / float(r["d_secs"]) if float(r["d_secs"] or 0) > 0 else 0.0
            rss = int(r["c_rss_kb"])
            det = r.get("deterministic", "n/a")
            print(
                f"{r['tool']:<13} {fmt_bytes(cb):>10} {r['ratio']:>7} "
                f"{vs_gzip:>7.2f}x {bpb:>10.3f} {c_mb:>8.1f} {d_mb:>8.1f} "
                f"{(fmt_bytes(rss*1024) if rss>=0 else 'n/a'):>8} "
                f"{r['rt_ok']:>4} {det:>4}"
            )
            # Per-stream breakdown (fqxv rows carry names/seq/qual bytes; other
            # tools store -1). Shows where the bits go — the lever for tuning.
            nb, sb, qb = (int(r.get(k, -1) or -1) for k in ("names_bytes", "seq_bytes", "qual_bytes"))
            if nb >= 0 and sb >= 0 and qb >= 0 and (nb + sb + qb) > 0:
                tot = nb + sb + qb
                print(
                    f"{'':<13} └─ names {fmt_bytes(nb)} ({100*nb/tot:.0f}%)  "
                    f"seq {fmt_bytes(sb)} ({100*sb/tot:.0f}%)  "
                    f"qual {fmt_bytes(qb)} ({100*qb/tot:.0f}%)"
                )
            # Quality distortion (lossy rows only; -1 = lossless / not measured).
            # The fidelity half of the lossy tradeoff: how far the reconstructed
            # qualities drift from the originals, alongside the ratio above.
            mae = float(r.get("qual_mae", -1) or -1)
            if mae >= 0:
                rmse = float(r.get("qual_rmse", -1) or -1)
                pct = float(r.get("qual_pct_changed", -1) or -1)
                print(
                    f"{'':<13} └─ Δqual  mae {mae:.2f}  rmse {rmse:.2f}  "
                    f"changed {pct:.1f}% of bases"
                )

    # Optional section: fqxv vs the native .sra (from sra_compare.sh). Only shown
    # when that harness has been run and its table is present.
    sra_path = RESULTS_DIR / "sra_compare.tsv"
    if sra_path.exists():
        render_sra_compare(load_tsv(sra_path))


def render_sra_compare(rows: list[dict[str, str]]) -> None:
    """fqxv archive size vs the native .sra the run shipped in.

    sra_compare.tsv is long-format: one row per (accession, fqxv point). We print
    the shared .sra/.fastq sizes once per run, then each fqxv point's size and its
    ratio to the .sra (< 1.00 = fqxv is smaller). A closing line geomeans the
    lossless (`max`) fqxv/.sra across the panel — the headline "is it worth it".
    """
    by_acc: dict[str, list[dict[str, str]]] = defaultdict(list)
    for r in rows:
        by_acc[r["accession"]].append(r)

    print("\n\n## fqxv vs native .sra  (the archive the data actually shipped in)")
    hdr = (
        f"{'accession':<12} {'platform':<14} {'regime':<11} {'.sra':>9} "
        f"{'.fastq':>9} {'point':<6} {'fqxv':>9} {'fqxv/.sra':>10} {'fqxv/.fastq':>12}"
    )
    print(hdr)
    print("-" * len(hdr))

    lossless_ratios: list[float] = []
    for acc, prows in by_acc.items():
        # Stable point order (lossless first, then increasingly lossy).
        order = {"max": 0, "bin8": 1, "bin4": 2, "bin2": 3}
        prows = sorted(prows, key=lambda r: order.get(r["point"], 9))
        for i, r in enumerate(prows):
            sra_b = int(r.get("sra_bytes", -1) or -1)
            fq_b = int(r.get("fastq_bytes", 0) or 0)
            fx_b = int(r.get("fqxv_bytes", 0) or 0)
            over_sra = r.get("fqxv_over_sra", "NA")
            over_fq = r.get("fqxv_over_fastq", "NA")
            if r["point"] == "max" and over_sra not in ("NA", ""):
                lossless_ratios.append(float(over_sra))
            first = i == 0
            print(
                f"{(acc if first else ''):<12} {(r['platform'] if first else ''):<14} "
                f"{(r['regime'] if first else ''):<11} "
                f"{(fmt_bytes(sra_b) if first and sra_b >= 0 else ''):>9} "
                f"{(fmt_bytes(fq_b) if first else ''):>9} "
                f"{r['point']:<6} {fmt_bytes(fx_b):>9} {over_sra:>10} {over_fq:>12}"
            )
        print("-" * len(hdr))

    if lossless_ratios:
        geo = math.exp(sum(math.log(x) for x in lossless_ratios) / len(lossless_ratios))
        verdict = "smaller than" if geo < 1 else "larger than"
        print(
            f"lossless fqxv (--max) is {geo:.2f}x the .sra size on average "
            f"(geomean over {len(lossless_ratios)} runs) — {1/geo:.2f}x {verdict} .sra."
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
