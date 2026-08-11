#!/usr/bin/env bash
# Stage the PacBio HiFi benchmark dataset: fetch SRR11434954 (Wick E. coli
# E2348/69 CCS, BioProject PRJNA602597 — ~44 GB raw FASTQ) and keep the first
# 120k reads (~1.5 Gbp, ~300x) as ecoli_hifi.fastq. The prefix take
# deliberately reproduces the retired longread_hifi.sbatch recipe
# (`git show bf421e3:bench/longread_hifi.sbatch`) so new numbers stay
# comparable with the historical HiFi results.
#
# Fetching is network/IO bound and fine on a login node (same policy as
# fetch.sh); anything compute-heavy must go through srun/sbatch.
set -euo pipefail

DATA_DIR="${FQXV_DATA_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/data}"
THREADS="${FQXV_FETCH_THREADS:-8}"
READS="${HIFI_SUBSAMPLE_READS:-120000}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# `pixi run` has no manifest to find here; use the env's binaries directly.
# Worktrees don't carry the gitignored bench/.pixi env — point FQXV_PIXI_BIN
# at a checkout that has one when running from a worktree.
PIXI_BIN="${FQXV_PIXI_BIN:-$HERE/../.pixi/envs/default/bin}"
export PATH="$PIXI_BIN:$PATH"

mkdir -p "$DATA_DIR"

if [[ -s "$DATA_DIR/ecoli_hifi.fastq" ]]; then
  echo "[skip] ecoli_hifi.fastq already staged"
  exit 0
fi

if ! compgen -G "$DATA_DIR/SRR11434954*.fastq" > /dev/null; then
  avail_gb=$(df --output=avail -BG "$DATA_DIR" | tail -1 | tr -dc '0-9')
  if (( avail_gb < 60 )); then
    echo "error: need >=60 GB free in $DATA_DIR (have ${avail_gb}G)" >&2
    exit 1
  fi
  sracha get SRR11434954 \
    --output-dir "$DATA_DIR" \
    --threads "$THREADS" \
    --split split-3 \
    --no-gzip \
    --no-progress
fi

# sracha writes single-end runs as <acc>.fastq or <acc>_0.fastq.
raw="$DATA_DIR/SRR11434954.fastq"
[[ -s "$raw" ]] || raw="$DATA_DIR/SRR11434954_0.fastq"
[[ -s "$raw" ]] || { echo "error: fetch produced no FASTQ" >&2; exit 1; }

seqkit head -n "$READS" "$raw" -o "$DATA_DIR/ecoli_hifi.fastq"
seqkit stats -a "$DATA_DIR/ecoli_hifi.fastq"

# The 44 GB raw file is kept by default (scratch is roomy); set
# HIFI_KEEP_RAW=0 to drop it after subsampling.
if [[ "${HIFI_KEEP_RAW:-1}" == "0" ]]; then
  rm -f "$raw"
fi
