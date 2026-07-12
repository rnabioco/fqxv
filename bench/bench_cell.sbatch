#!/usr/bin/env bash
#SBATCH --job-name=fqxv-cell
#SBATCH --comment=fqxv-cell
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32       # 32 cores, not whole-node --exclusive, so cells
#SBATCH --mem=90G                # pack onto partially-free nodes and queue faster
#SBATCH --time=02:00:00
#SBATCH --output=fqxv-cell-%A_%a.out
## #SBATCH --account=amc-general
#
# Phase 2 of the parallel harness: one array task per (dataset, tool) cell. Each
# cell still runs one measured tool at a time on its 32 bound cores — the point
# is never to co-schedule tools within a cell (every tool already uses all its
# cores, so that would only oversubscribe and blur MB/s). Requesting 32 cores
# rather than a whole exclusive node schedules much faster; the cost is that a
# co-tenant on the other socket can add some memory-bandwidth noise to the
# throughput columns (ratios and rt/det correctness are unaffected). For a
# publication-clean throughput run, add `#SBATCH --exclusive` back.
# The cell index maps into $RESULTS_DIR/cells.tsv (submit_parallel.sh generated
# it); results append to a private part file.
set -euo pipefail

export PATH="$HOME/.pixi/bin:$HOME/.cargo/bin:$PATH"
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
RESULTS_DIR="${FQXV_RESULTS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/results}"
CELLS="$RESULTS_DIR/cells.tsv"

line="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$CELLS")"
[[ -n "$line" ]] || { echo "no cell at index ${SLURM_ARRAY_TASK_ID}"; exit 1; }
label="$(cut -f1 <<<"$line")"
tool="$(cut -f2 <<<"$line")"

export FQXV_THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
export RAYON_NUM_THREADS="$FQXV_THREADS"
export FQXV_ONLY_DATASET="$label"
export FQXV_TOOLS="$tool"
export FQXV_PART_TAG="${label}.${tool}"

echo "==> cell ${SLURM_ARRAY_TASK_ID}: ${label} / ${tool} on $(hostname), ${FQXV_THREADS} threads, $(date)"
pixi run --manifest-path "$HERE/pixi.toml" bash "$HERE/run_bench.sh"
echo "==> cell done $(date)"
