#!/usr/bin/env bash
#SBATCH --job-name=fqxv-cell
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --exclusive              # one measured tool per node: no co-tenant
#SBATCH --mem=0                  # memory-bandwidth noise, clean throughput
#SBATCH --time=02:00:00
#SBATCH --output=fqxv-cell-%A_%a.out
## #SBATCH --account=amc-general
#
# Phase 2 of the parallel harness: one array task per (dataset, tool) cell, each
# on its own exclusive node so the throughput numbers stay clean (running tools
# concurrently on a shared node would halve and blur MB/s — every tool already
# uses all the cores). The cell index maps into $RESULTS_DIR/cells.tsv, which
# submit_parallel.sh generated. Results append to a private part file.
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
