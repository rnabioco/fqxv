#!/usr/bin/env python3
"""Render the fqxv-vs-native-.sra size comparison as static PNGs.

Reads ``docs/charts/sra.tsv`` (one row per SRA run) and writes a light and a
dark PNG into ``docs/images/`` (``sra_vs_native.light.png`` / ``.dark.png``).
Grouped horizontal bars: one group per run, three bars per group — the raw
``.fastq`` a user downloads, the native ``.sra`` NCBI ships, and the lossless
``fqxv --max`` archive. Every bar is direct-labelled with its size, so the
message ("fqxv is the smallest of the three, on every platform") reads without
axes.

Same rationale as ``generate.py``: PNG (not ``<img>``-embedded SVG) so long
axis labels don't reflow and clip; ``vl-convert-python`` compiles the Vega-Lite
spec with no browser. Run via ``pixi run charts`` (renders this and the topline
charts). Data refreshed by ``bench/scripts/sra_compare.sh``.

Palette: the categorical slots blue / orange / aqua from the validated default
(``dataviz`` skill) — colorblind-safe as a trio in both modes. fqxv takes aqua
(nearest the brand teal, the hero series). The light-mode aqua sits just under
3:1 on the surface, so the direct value labels are load-bearing (the relief
rule), not decoration.
"""

from __future__ import annotations

import csv
import json
import pathlib

import vl_convert as vlc

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent.parent
DATA = HERE / "sra.tsv"
OUT = REPO / "docs" / "images"

FONT = "-apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif"

# Series order (top→bottom within a group) and the format label each bar wears.
SERIES = [
    ("fastq", "FASTQ"),
    ("sra", ".sra (NCBI)"),
    ("fqxv", "fqxv --max"),
]

# Per-theme palette. Categorical slots blue/orange/aqua, stepped per surface;
# text tokens mirror docs/stylesheets/fqxv.css and the dataviz reference.
THEMES = {
    "light": {
        "surface": "#fcfcfb",
        "text": "#41564f",
        "colors": {"fastq": "#2a78d6", "sra": "#eb6834", "fqxv": "#1baf7a"},
    },
    "dark": {
        "surface": None,  # transparent; the site card supplies the dark surface
        "text": "#c3d6d1",
        "colors": {"fastq": "#3987e5", "sra": "#d95926", "fqxv": "#199e70"},
    },
}


def fmt_bytes(n: int) -> str:
    """Decimal size, matching the CLI summary (GB at/above 1e9, else MB)."""
    return f"{n / 1e9:.2f} GB" if n >= 1e9 else f"{n / 1e6:.0f} MB"


def load_rows() -> list[dict]:
    with DATA.open(newline="") as fh:
        reader = csv.DictReader(
            (line for line in fh if not line.startswith("#")), delimiter="\t"
        )
        return list(reader)


def spec_for(rows: list[dict], theme: dict) -> dict:
    order = [f"{r['platform']}\n{r['accession']}" for r in rows]  # group order = TSV order
    fmt_domain = [name for _, name in SERIES]
    color_range = [theme["colors"][key] for key, _ in SERIES]

    values = []
    vmax = 0.0
    for r in rows:
        group = f"{r['platform']}\n{r['accession']}"
        for key, name in SERIES:
            b = int(r[f"{key}_bytes"])
            gb = b / 1e9
            vmax = max(vmax, gb)
            values.append(
                {"group": group, "format": name, "gb": gb, "vlabel": fmt_bytes(b)}
            )

    return {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "background": theme["surface"],
        "width": 500,
        "height": {"step": 26},
        "padding": 6,
        "data": {"values": values},
        "encoding": {
            "y": {
                "field": "group",
                "type": "nominal",
                "sort": order,
                # paddingInner separates the accession clusters; without it the
                # three offset sub-bars space evenly and the grouping is lost.
                "scale": {"paddingInner": 0.45, "paddingOuter": 0.25},
                "axis": {
                    "title": None,
                    "labelColor": theme["text"],
                    "labelFontSize": 13,
                    "labelFontWeight": "bold",
                    "labelPadding": 10,
                    "labelExpr": "split(datum.label, '\\n')",
                    "domain": False,
                    "ticks": False,
                },
            }
        },
        "layer": [
            {
                "mark": {"type": "bar", "height": 16, "cornerRadiusEnd": 3},
                "encoding": {
                    "x": {
                        "field": "gb",
                        "type": "quantitative",
                        "scale": {"domain": [0, vmax * 1.2]},
                        "axis": None,
                    },
                    "yOffset": {"field": "format", "sort": fmt_domain},
                    "color": {
                        "field": "format",
                        "type": "nominal",
                        "scale": {"domain": fmt_domain, "range": color_range},
                        "legend": {
                            "orient": "top",
                            "direction": "horizontal",
                            "title": None,
                            "labelColor": theme["text"],
                            "labelFontSize": 12,
                            "symbolType": "square",
                            "symbolSize": 130,
                        },
                    },
                },
            },
            {
                "mark": {
                    "type": "text",
                    "align": "left",
                    "baseline": "middle",
                    "dx": 5,
                    "color": theme["text"],
                    "fontSize": 11,
                },
                "encoding": {
                    "x": {"field": "gb", "type": "quantitative", "scale": {"domain": [0, vmax * 1.2]}},
                    "yOffset": {"field": "format", "sort": fmt_domain},
                    "text": {"field": "vlabel"},
                },
            },
        ],
        "config": {"font": FONT, "view": {"stroke": None}},
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    for name, theme in THEMES.items():
        png = vlc.vegalite_to_png(json.dumps(spec_for(rows, theme)), scale=2)
        dest = OUT / f"sra_vs_native.{name}.png"
        dest.write_bytes(png)
        print(f"wrote {dest.relative_to(REPO)}")


if __name__ == "__main__":
    main()
