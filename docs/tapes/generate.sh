#!/bin/bash
# Generate all VHS gifs for the docs. Run from anywhere; it cds to repo root.
#
#   pixi run tapes           # builds the release binary first, then renders
#   ./docs/tapes/generate.sh # if fqxv + vhs are already on PATH
#
# vhs writes the final GIFs directly — no post-processing. The tapes use the
# Source Code Pro monospace font shipped by the pixi env (VHS's headless
# chromium reads it via fontconfig); don't rename it to a font that isn't
# installed or glyphs render with broken metrics.

set -euo pipefail

cd "$(dirname "$0")/../.."

# Prefer the pixi-built release binary; fall back to a plain cargo build so the
# script also works outside pixi.
if [[ ! -x target/release/fqxv ]]; then
    cargo build --profile release
fi
export PATH="$PWD/target/release:$PATH"
# Silence the library's info-level diagnostics (e.g. the "detected layout" span)
# so the tapes show only the user-facing run summary. The summary is printed via
# eprintln! and is independent of the log level, so it still appears. VHS's shell
# inherits this env, just like PATH.
export RUST_LOG=warn

mkdir -p docs/images

# --- self-contained demo fixture -------------------------------------------
# A small synthetic FASTQ (30k reads x 100 bp) so the tapes render identically
# anywhere without downloading anything. Modeled loosely on real Illumina data —
# reads sampled from a short reference (cross-read redundancy) with occasional
# mismatches, and low-entropy quality skewed toward high scores — so the demo
# ratios are representative rather than worst-case random. Ignored by git
# (*.fastq.gz / *.fqxv) and cleaned up below.
echo "Building demo fixture…"
awk 'BEGIN {
    srand(42);
    b = "ACGT";
    # Quality alphabet weighted toward high Illumina scores (mostly F/J).
    qs = "FFFFFFFFFFFFFFFFFFFFFFJJJJJJJJJJ::::,,##";
    # A 4 kbp reference; reads are 100 bp windows into it, so the sequence
    # stream carries realistic cross-read redundancy.
    ref = "";
    for (k = 0; k < 4000; k++) ref = ref substr(b, int(rand() * 4) + 1, 1);
    reflen = length(ref);
    for (i = 0; i < 30000; i++) {
        start = int(rand() * (reflen - 100)) + 1;
        s = substr(ref, start, 100);
        # ~1% per-base sequencing error so reads are not identical copies.
        q = "";
        for (j = 1; j <= 100; j++) {
            if (rand() < 0.01) {
                s = substr(s, 1, j - 1) substr(b, int(rand() * 4) + 1, 1) substr(s, j + 1);
            }
            q = q substr(qs, int(rand() * length(qs)) + 1, 1);
        }
        printf "@READ%d\n%s\n+\n%s\n", i, s, q;
    }
}' | gzip > demo.fastq.gz

# Pre-build the archive so the info/decompress tapes have it before compress
# runs (glob order is alphabetical, but don't depend on it).
fqxv --quiet compress demo.fastq.gz -o demo.fqxv

# --- render -----------------------------------------------------------------
for tape in docs/tapes/*.tape; do
    echo "Generating: $tape"
    vhs "$tape"
done

# --- cleanup ----------------------------------------------------------------
rm -f demo.fastq.gz demo.fqxv restored.fastq

echo "Done — gifs are in docs/images/"
