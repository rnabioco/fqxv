#!/usr/bin/env bash
# Fetch benchmark FASTQ to $SCRATCH/fqxv/data using sracha (raw, uncompressed).
#
# We download uncompressed FASTQ so every compressor starts from the identical
# source bytes; the harness synthesizes the .fastq.gz "baseline of record"
# itself. Run from the bench/ dir under `pixi run` (or `pixi shell`):
#
#   pixi run bash fetch.sh                 # fetch all datasets in datasets.tsv
#   pixi run bash fetch.sh SRR2627175      # fetch just one accession
#
# Data staging is network/IO bound and fine on a login node; the compression
# *benchmark* itself must go through srun/sbatch (see bench.sbatch).
set -euo pipefail

DATA_DIR="${FQXV_DATA_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/data}"
THREADS="${FQXV_FETCH_THREADS:-8}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

mkdir -p "$DATA_DIR"

# Accessions: CLI args override the manifest.
if [[ $# -gt 0 ]]; then
  accessions=("$@")
else
  mapfile -t accessions < <(grep -v '^#' "$HERE/../panels/datasets.tsv" | awk 'NF{print $1}')
fi

echo "==> fetching ${#accessions[@]} accession(s) to $DATA_DIR"
for acc in "${accessions[@]}"; do
  if compgen -G "$DATA_DIR/${acc}*.fastq" > /dev/null; then
    echo "  [skip] $acc (already present)"
    continue
  fi
  echo "  [get ] $acc"
  # Per-accession failures (e.g. cSRA archives sracha can't reconstruct) must
  # not abort the whole batch.
  if ! sracha get "$acc" \
      --output-dir "$DATA_DIR" \
      --threads "$THREADS" \
      --split split-3 \
      --no-gzip \
      --no-progress; then
    echo "  [FAIL] $acc — skipping" >&2
    continue
  fi
done

echo "==> done. Contents:"
ls -lh "$DATA_DIR"/*.fastq 2>/dev/null || true
