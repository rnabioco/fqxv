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
# CoLoRd and platform-appropriate quality bins. The sets themselves live in
# toolsets.sh, shared with run_bench.sh so the two drivers cannot drift apart.
# shellcheck source=./toolsets.sh
. "$HERE/toolsets.sh"
# Partition/QoS/account for whichever cluster this is (Bodhi `rna` by default,
# Alpine `amilan` when detected or FQXV_CLUSTER=alpine). Passed explicitly to
# every sbatch below so the submission does not depend on each .sbatch file's
# own default — those had drifted and a plain run on Bodhi queued against
# `amilan`, which does not exist there.
# shellcheck source=./cluster.sh
. "$HERE/cluster.sh"
# FQXV_TOOLS overrides every platform's set with one explicit list (no platform
# filtering) — for running a specific tool everywhere. FQXV_LR_TOOLS_ONT /
# FQXV_LR_TOOLS_HIFI override just that platform.
TOOLS="${FQXV_TOOLS:-}"
LR_TOOLS_ONT="${FQXV_LR_TOOLS_ONT:-}"
LR_TOOLS_HIFI="${FQXV_LR_TOOLS_HIFI:-}"
MAXPAR="${FQXV_MAXPAR:-}"   # empty = no cap; the scheduler manages concurrency

mkdir -p "$RESULTS_DIR"
# Slurm writes job output relative to the submission dir, so the #SBATCH --output
# directives would litter the repo root. Directives cannot expand variables, so
# the submitter overrides them here and keeps every log with its results.
LOG_DIR="$RESULTS_DIR/logs"
mkdir -p "$LOG_DIR"
CELLS="$RESULTS_DIR/cells.tsv"
: > "$CELLS"
# datasets.tsv columns: accession label platform layout quality approx_gz reference
while read -r acc label platform layout quality approx _; do
  [[ -z "${acc:-}" || "$acc" == \#* ]] && continue
  # Optional read-class filter: FQXV_ONLY=long|short restricts the matrix to one
  # class (e.g. `FQXV_ONLY=long` for a CoLoRd head-to-head without the short-read
  # field). Empty runs everything.
  [[ -n "${FQXV_ONLY:-}" && "$(fqxv_read_class "$platform")" != "$FQXV_ONLY" ]] && continue
  # Tool set by platform, unless explicitly overridden.
  if [[ -n "$TOOLS" ]]; then
    sel="$TOOLS"
  else
    case "$platform" in
      MinION|GridION|PromethION|*[Nn]anopore*|ONT)
        sel="${LR_TOOLS_ONT:-$FQXV_TOOLSET_ONT}" ;;
      SequelII|Sequel*|Revio|*[Hh]iFi*|PacBio*)
        sel="${LR_TOOLS_HIFI:-$FQXV_TOOLSET_HIFI}" ;;
      *)
        sel="$(fqxv_toolset_for_platform "$platform")" ;;
    esac
  fi
  # Require data now, EXCEPT datasets prep stages itself (approx_gz='subsampled',
  # e.g. HiFi) — prep runs before the array via the afterok dependency.
  # Single-end runs land as either `${acc}.fastq` or `${acc}_0.fastq` depending
  # on the run's member layout (sracha decides, not the --split flag), so test
  # both — matching run_bench.sh. Missing the `_0` form silently dropped whole
  # datasets from the matrix behind a single "[skip]" line.
  if [[ "$approx" != "subsampled" \
        && ! -f "$DATA_DIR/${acc}_1.fastq" && ! -f "$DATA_DIR/${acc}.fastq" \
        && ! -f "$DATA_DIR/${acc}_0.fastq" ]]; then
    echo "[skip] $label ($acc): no data in $DATA_DIR"
    continue
  fi
  for t in $sel; do
    printf '%s\t%s\n' "$label" "$t" >> "$CELLS"
  done
done < <(grep -v '^#' "$HERE/../panels/datasets.tsv")

N="$(wc -l < "$CELLS")"
[[ "$N" -gt 0 ]] || { echo "no cells generated (no datasets present?)"; exit 1; }
echo "==> $N cells ($(wc -l < <(grep -v '^#' "$HERE/../panels/datasets.tsv" | awk 'NF')) datasets x tools) -> $CELLS"
echo "==> cluster $FQXV_CLUSTER_RESOLVED: partition=$FQXV_PARTITION qos=$FQXV_QOS${FQXV_ACCOUNT:+ account=$FQXV_ACCOUNT}"

# Fresh run: clear previous parts so the merge only sees this submission.
rm -f "$RESULTS_DIR/parts"/results.*.tsv "$RESULTS_DIR/parts"/meta.*.tsv 2>/dev/null || true

jid_prep=$(sbatch --parsable "${FQXV_SBATCH_OPTS[@]}" \
                  --output="$LOG_DIR/%x-%j.out" "$HERE/../slurm/prep.sbatch")
echo "prep   -> job $jid_prep"
# Default: no concurrency cap — let the scheduler decide how many cells run at
# once (each cell is one exclusive-ish node, so Slurm's partition limits already
# bound it). Set FQXV_MAXPAR only to *deliberately* throttle.
array="1-${N}"
[[ -n "$MAXPAR" ]] && array="${array}%${MAXPAR}"
jid_arr=$(sbatch --parsable "${FQXV_SBATCH_OPTS[@]}" --dependency=afterok:"$jid_prep" \
                 --output="$LOG_DIR/%x-%A_%a.out" \
                 --array="$array" "$HERE/../slurm/bench_cell.sbatch")
echo "cells  -> job $jid_arr  (array $array)"
jid_merge=$(sbatch --parsable "${FQXV_SBATCH_OPTS[@]}" --dependency=afterany:"$jid_arr" \
                   --output="$LOG_DIR/%x-%j.out" "$HERE/../slurm/merge.sbatch")
echo "merge  -> job $jid_merge"
echo
echo "watch:   squeue -u $USER"
echo "results: $RESULTS_DIR/results.tsv  (rendered in fqxv-merge-${jid_merge}.out)"
