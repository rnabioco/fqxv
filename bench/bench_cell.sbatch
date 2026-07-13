#!/usr/bin/env bash
#SBATCH --job-name=fqxv-cell
#SBATCH --comment=fqxv-cell
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16       # Memory is the real floor: the worst cell (xz9
#SBATCH --mem=28G                # on NovaSeq) peaks ~22.5G, and amilan caps mem
#SBATCH --time=01:00:00          # at 3840M/core, so 28G needs >=8 cores; 16 adds
#SBATCH --output=fqxv-cell-%A_%a.out                # throughput headroom while
## #SBATCH --account=amc-general                     # still packing onto a node.
#
# Phase 2 of the parallel harness: one array task per (dataset, tool) cell. Each
# cell runs one measured tool at a time on its bound cores — the point is never
# to co-schedule tools within a cell (every tool already uses all its cores, so
# that would only oversubscribe and blur MB/s). The ask is sized to memory, not
# throughput: xz9/fqxv-reorder on the 5.5G NovaSeq set peak ~16-22.5G, and this
# partition hard-caps memory at 3840M/core, so 28G already pulls in ~8 cores.
# Fewer cores means every cell packs onto a partially-free node and schedules
# far faster; the cost is lower absolute MB/s in the throughput columns (ratios
# and rt/det correctness are core-count-independent). Tools use fewer per-thread
# buffers at 8 threads, so peak memory is *below* the 32-thread figures above.
# For a publication-clean throughput run, raise cpus-per-task and add
# `#SBATCH --exclusive`.
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
