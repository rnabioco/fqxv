#!/usr/bin/env bash
# Submit the benchmark as a parallel Slurm pipeline:
#
#   prep (build + shared input digests)  ->  array of (dataset x tool) cells,
#   each on its own exclusive node  ->  merge (assemble results.tsv + report).
#
# One measured tool per node keeps throughput clean while the whole matrix fans
# out across the cluster. Run from the login node: `bash submit_parallel.sh`.
#
# Env knobs:
#   FQXV_TOOLS    tool subset (default: all)
#   FQXV_MAXPAR   optional cap on concurrent array tasks (default: none — the
#                 scheduler decides; set only to deliberately throttle)
#   FQXV_DATA_DIR / FQXV_RESULTS_DIR  as in run_bench.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${FQXV_RESULTS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/results}"
DATA_DIR="${FQXV_DATA_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/data}"
# Tool sets are platform-aware, so short- and long-read datasets share one matrix
# (this replaces the old standalone longread*.sbatch serial loops). Short-read
# (Illumina) gets the full field matrix; ONT/PacBio get the long-read set with
# CoLoRd and platform-appropriate quality bins.
TOOLS="${FQXV_TOOLS:-fqxv fqxv9 fqxv-reorder fqxv-bin4 gzip zstd19 xz9 fqz_comp fqzcomp5 spring}"
LR_TOOLS_ONT="${FQXV_LR_TOOLS_ONT:-fqxv fqxv9 fqxv-max fqxv-bin4 fqxv-bin2 fqxv-binont gzip zstd19 xz9 fqz_comp colord colord-lossy}"
LR_TOOLS_HIFI="${FQXV_LR_TOOLS_HIFI:-fqxv fqxv9 fqxv-max fqxv-binhifi fqxv-binont gzip zstd19 xz9 colord colord-lossy}"
MAXPAR="${FQXV_MAXPAR:-}"   # empty = no cap; the scheduler manages concurrency

mkdir -p "$RESULTS_DIR"
CELLS="$RESULTS_DIR/cells.tsv"
: > "$CELLS"
# datasets.tsv columns: accession label platform layout quality approx_gz reference
while read -r acc label platform layout quality approx _; do
  [[ -z "${acc:-}" || "$acc" == \#* ]] && continue
  # Tool set by platform.
  case "$platform" in
    MinION|GridION|PromethION|*[Nn]anopore*|ONT) sel="$LR_TOOLS_ONT" ;;
    SequelII|Sequel*|Revio|*[Hh]iFi*|PacBio*)    sel="$LR_TOOLS_HIFI" ;;
    *)                                           sel="$TOOLS" ;;
  esac
  # Require data now, EXCEPT datasets prep stages itself (approx_gz='subsampled',
  # e.g. HiFi) — prep runs before the array via the afterok dependency.
  if [[ "$approx" != "subsampled" \
        && ! -f "$DATA_DIR/${acc}_1.fastq" && ! -f "$DATA_DIR/${acc}.fastq" ]]; then
    echo "[skip] $label ($acc): no data in $DATA_DIR"
    continue
  fi
  for t in $sel; do
    printf '%s\t%s\n' "$label" "$t" >> "$CELLS"
  done
done < <(grep -v '^#' "$HERE/datasets.tsv")

N="$(wc -l < "$CELLS")"
[[ "$N" -gt 0 ]] || { echo "no cells generated (no datasets present?)"; exit 1; }
echo "==> $N cells ($(wc -l < <(grep -v '^#' "$HERE/datasets.tsv" | awk 'NF')) datasets x tools) -> $CELLS"

# Fresh run: clear previous parts so the merge only sees this submission.
rm -f "$RESULTS_DIR/parts"/results.*.tsv "$RESULTS_DIR/parts"/meta.*.tsv 2>/dev/null || true

jid_prep=$(sbatch --parsable "$HERE/prep.sbatch")
echo "prep   -> job $jid_prep"
# Default: no concurrency cap — let the scheduler decide how many cells run at
# once (each cell is one exclusive-ish node, so Slurm's partition limits already
# bound it). Set FQXV_MAXPAR only to *deliberately* throttle.
array="1-${N}"
[[ -n "$MAXPAR" ]] && array="${array}%${MAXPAR}"
jid_arr=$(sbatch --parsable --dependency=afterok:"$jid_prep" \
                 --array="$array" "$HERE/bench_cell.sbatch")
echo "cells  -> job $jid_arr  (array $array)"
jid_merge=$(sbatch --parsable --dependency=afterany:"$jid_arr" "$HERE/merge.sbatch")
echo "merge  -> job $jid_merge"
echo
echo "watch:   squeue -u $USER"
echo "results: $RESULTS_DIR/results.tsv  (rendered in fqxv-merge-${jid_merge}.out)"
