#!/usr/bin/env python3
"""Sync the Bioconda / BioContainers version pins in README.md and docs/.

The docs pin a concrete BioContainers tag (`fqxv:<version>--<build>`) because
quay.io publishes no `latest` tag for Bioconda-derived images — a tag has to
name a version *and* a conda build hash that actually exist. That makes the
pins rot every release, so this script looks up what is current and rewrites
them in place.

Two independent sources, because they can legitimately disagree:

  * the container tag comes from quay.io, so we only ever write a tag that has
    really been published (BioContainers lags the Bioconda package by hours to
    days);
  * the version pins (`bioconda::fqxv=X.Y.Z`, `pixi add`, `environment.yml`)
    come from anaconda.org, which is the package `conda`/`pixi` would resolve.

Usage:
    sync-bioconda-version.py            # rewrite the files in place
    sync-bioconda-version.py --check    # exit 1 if anything is out of date
"""

from __future__ import annotations

import argparse
import json
import re
import sys
import urllib.error
import urllib.request
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
TARGETS = ["README.md", "docs/index.md", "docs/getting-started/installation.md"]

QUAY_TAGS = (
    "https://quay.io/api/v1/repository/biocontainers/fqxv/tag/"
    "?limit=100&onlyActiveTags=true"
)
ANACONDA_PKG = "https://api.anaconda.org/package/bioconda/fqxv"

# `quay.io/biocontainers/fqxv:0.7.0--h1234567_0` and the matching Singularity
# URL on depot.galaxyproject.org. The tag is matched loosely (everything up to
# the next quote or whitespace) so the `<version>--<build>` placeholder the
# docs ship with — written before BioContainers has published anything — gets
# replaced by the first successful run just like a stale real tag would.
CONTAINER_TAG_RE = re.compile(r"(?<=/fqxv:)[^\s'\"`]+")
# `conda 'bioconda::fqxv=0.7.0'` and `pixi add "bioconda::fqxv==0.7.0"`
CONDA_PIN_RE = re.compile(r"(?<=bioconda::fqxv==)(\d[\d.]*)|(?<=bioconda::fqxv=)(\d[\d.]*)")
# `fqxv = { version = "==0.7.0", channel = "bioconda" }` and `- fqxv=0.7.0`
TOML_PIN_RE = re.compile(r'(?<=fqxv = \{ version = "==)(\d[\d.]*)')
YAML_PIN_RE = re.compile(r"(?<=^  - fqxv=)(\d[\d.]*)", re.MULTILINE)


def fetch_json(url: str) -> dict:
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        return json.load(resp)


def version_key(version: str) -> tuple[int, ...]:
    return tuple(int(part) for part in version.split("."))


def latest_container_tag() -> str | None:
    """Newest `<version>--<build>` tag published on quay.io, if any.

    Returns None when the repository does not exist yet: quay answers 401 for
    an unknown public repo, and a brand-new Bioconda package has no image for
    a day or two. That is not a failure — leave the placeholder in place and
    try again next week.
    """
    try:
        payload = fetch_json(QUAY_TAGS)
    except urllib.error.HTTPError as err:
        if err.code in (401, 404):
            print("no BioContainers image published yet; leaving container pins alone")
            return None
        raise
    tags = []
    for tag in payload.get("tags", []):
        match = re.fullmatch(r"(\d[\d.]*)--([A-Za-z0-9_]+)", tag.get("name", ""))
        if match:
            tags.append((version_key(match.group(1)), match.group(0)))
    if not tags:
        print("no versioned tags on quay.io/biocontainers/fqxv yet")
        return None
    return max(tags)[1]


def latest_conda_version() -> str:
    """Newest version of the `fqxv` package on anaconda.org/bioconda."""
    version = fetch_json(ANACONDA_PKG).get("latest_version")
    if not version:
        raise SystemExit("anaconda.org returned no latest_version for bioconda::fqxv")
    return version


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="report drift and exit 1 instead of rewriting files",
    )
    args = parser.parse_args()

    try:
        container_tag = latest_container_tag()
        conda_version = latest_conda_version()
    except (urllib.error.URLError, TimeoutError) as err:
        print(f"error: upstream lookup failed: {err}", file=sys.stderr)
        return 2

    print(f"biocontainers tag: {container_tag or '(none published)'}")
    print(f"bioconda version:  {conda_version}")

    stale = []
    for relative in TARGETS:
        path = REPO_ROOT / relative
        original = path.read_text()
        updated = original
        if container_tag:
            updated = CONTAINER_TAG_RE.sub(container_tag, updated)
        updated = CONDA_PIN_RE.sub(conda_version, updated)
        updated = TOML_PIN_RE.sub(conda_version, updated)
        updated = YAML_PIN_RE.sub(conda_version, updated)
        if updated == original:
            continue
        stale.append(relative)
        if not args.check:
            path.write_text(updated)

    if not stale:
        print("pins are up to date")
        return 0

    if args.check:
        print(f"out of date: {', '.join(stale)}", file=sys.stderr)
        return 1

    print(f"updated: {', '.join(stale)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
