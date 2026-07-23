#!/usr/bin/env python3
"""Render the topline benchmark bar charts as static PNGs.

Reads the committed source of truth ``docs/charts/topline.tsv`` and writes one
light and one dark PNG per chart into ``docs/images/`` (e.g.
``illumina_ratio_novaseq.light.png`` / ``.dark.png``). The Markdown references
both variants and ``docs/stylesheets/fqxv.css`` shows the one matching the
active ``data-md-color-scheme``, so the charts track the site's light/dark
toggle.

PNG rather than SVG on purpose: an ``<img>``-embedded SVG's ``<text>`` is
re-rendered by the browser with its own font, whose metrics differ from the
renderer's, so long axis labels overflow the left margin and clip. Rasterizing
bakes text to pixels — WYSIWYG.

Aesthetic: the ruff benchmark look — simple horizontal bars, sorted, one accent
(fqxv teal) with the rest of the field in a flat muted gray, the value drawn at
the end of each bar, and no gridlines/axes/title chrome.

Pure Python: stdlib ``csv`` for the data and ``vl-convert-python`` to compile a
Vega-Lite spec to PNG (no browser, no page JS). Run via ``pixi run charts``.
"""

from __future__ import annotations

import csv
import json
import pathlib

import vl_convert as vlc

HERE = pathlib.Path(__file__).resolve().parent
REPO = HERE.parent.parent
DATA = HERE / "topline.tsv"
OUT = REPO / "docs" / "images"

# A neutral sans stack; vl-convert falls back cleanly if a face is unavailable.
FONT = "-apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif"

# Per-theme palette. Only the text/field colors differ between light and dark;
# the fqxv accent teal reads on both. Values mirror docs/stylesheets/fqxv.css.
THEMES = {
    "dark": {"accent": "#3fd6c8", "field": "#41564f", "text": "#c3d6d1"},
    "light": {"accent": "#0e7c72", "field": "#c9d6d2", "text": "#41564f"},
}


def load_charts() -> dict[str, list[dict]]:
    """Group the TSV rows by their ``chart`` column, order preserved."""
    charts: dict[str, list[dict]] = {}
    with DATA.open(newline="") as fh:
        reader = csv.DictReader(
            (line for line in fh if not line.startswith("#")), delimiter="\t"
        )
        for row in reader:
            charts.setdefault(row["chart"], []).append(row)
    return charts


def spec_for(rows: list[dict], theme: dict[str, str]) -> dict:
    """Build a Vega-Lite spec for one chart under one theme palette."""
    values = [
        {
            "label": r["label"],
            "group": r["group"],
            "value": float(r["value"]),
            "vlabel": r["vlabel"],
        }
        for r in rows
    ]
    # Sort bars by value (largest at top) with an explicit category order — a
    # channel-based `sort` is silently dropped in a layered spec.
    values.sort(key=lambda v: v["value"], reverse=True)
    order = [v["label"] for v in values]
    # Headroom so the value text drawn past the bar end isn't clipped.
    vmax = max(v["value"] for v in values) * 1.14

    return {
        "$schema": "https://vega.github.io/schema/vega-lite/v5.json",
        "background": None,
        "width": 460,
        "height": {"step": 28},
        "padding": 4,
        "data": {"values": values},
        "encoding": {
            "y": {
                "field": "label",
                "type": "nominal",
                "sort": order,
                "axis": {
                    "title": None,
                    "labelColor": theme["text"],
                    "labelFontSize": 13,
                    "labelPadding": 8,
                    "domain": False,
                    "ticks": False,
                },
            },
            "x": {
                "field": "value",
                "type": "quantitative",
                "scale": {"domain": [0, vmax]},
                "axis": None,
            },
        },
        "layer": [
            {
                "mark": {"type": "bar", "height": 15, "cornerRadiusEnd": 3},
                "encoding": {
                    "color": {
                        "field": "group",
                        "type": "nominal",
                        "scale": {
                            "domain": ["fqxv", "other"],
                            "range": [theme["accent"], theme["field"]],
                        },
                        "legend": None,
                    }
                },
            },
            {
                "mark": {
                    "type": "text",
                    "align": "left",
                    "baseline": "middle",
                    "dx": 5,
                    "color": theme["text"],
                    "fontSize": 12,
                },
                "encoding": {"text": {"field": "vlabel"}},
            },
        ],
        "config": {
            "font": FONT,
            "view": {"stroke": None},
        },
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    charts = load_charts()
    for chart, rows in charts.items():
        for name, theme in THEMES.items():
            # PNG, not SVG: the axis labels are laid out with vl-convert's font
            # metrics, but a browser re-renders an <img>-embedded SVG's <text>
            # with its own (wider) font, overflowing the left margin and clipping
            # long labels. Rasterizing bakes text to pixels — WYSIWYG. scale=2
            # keeps it crisp on hi-dpi displays; CSS caps width at the column.
            png = vlc.vegalite_to_png(json.dumps(spec_for(rows, theme)), scale=2)
            dest = OUT / f"{chart}.{name}.png"
            dest.write_bytes(png)
            print(f"wrote {dest.relative_to(REPO)}")


if __name__ == "__main__":
    main()
