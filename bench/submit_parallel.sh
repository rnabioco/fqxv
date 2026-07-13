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
#   FQXV_MAXPAR   cap on concurrent array tasks (default 12)
#   FQXV_DATA_DIR / FQXV_RESULTS_DIR  as in run_bench.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULTS_DIR="${FQXV_RESULTS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/results}"
DATA_DIR="${FQXV_DATA_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/data}"
TOOLS="${FQXV_TOOLS:-fqxv fqxv9 fqxv-reorder fqxv-bin4 gzip zstd19 xz9 fqz_comp fqzcomp5 spring}"
MAXPAR="${FQXV_MAXPAR:-12}"

mkdir -p "$RESULTS_DIR"
CELLS="$RESULTS_DIR/cells.tsv"
: > "$CELLS"
while read -r acc label _; do
  [[ -z "${acc:-}" || "$acc" == \#* ]] && continue
  # Paired R1, or single-end `${acc}.fastq` (Nanopore etc.).
  if [[ ! -f "$DATA_DIR/${acc}_1.fastq" && ! -f "$DATA_DIR/${acc}.fastq" ]]; then
    echo "[skip] $label ($acc): no data in $DATA_DIR"
    continue
  fi
  for t in $TOOLS; do
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
jid_arr=$(sbatch --parsable --dependency=afterok:"$jid_prep" \
                 --array="1-${N}%${MAXPAR}" "$HERE/bench_cell.sbatch")
echo "cells  -> job $jid_arr  (array 1-$N, up to $MAXPAR concurrent)"
jid_merge=$(sbatch --parsable --dependency=afterany:"$jid_arr" "$HERE/merge.sbatch")
echo "merge  -> job $jid_merge"
echo
echo "watch:   squeue -u $USER"
echo "results: $RESULTS_DIR/results.tsv  (rendered in fqxv-merge-${jid_merge}.out)"
