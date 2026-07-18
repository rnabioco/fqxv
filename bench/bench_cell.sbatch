#!/usr/bin/env bash
#SBATCH --job-name=fqxv-cell
#SBATCH --comment=fqxv-cell
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16       # Memory is the real floor: short-read xz9 on
#SBATCH --mem=48G                # NovaSeq peaks ~22.5G, and the long-read cells
#SBATCH --time=02:00:00          # (HiFi fqxv ~11G, CoLoRd on 1.5 Gbp more) want
                                 # headroom — 48G covers both. On Bodhi's rna
                                 # partition (no per-core mem cap) override the
                                 # amilan defaults via SBATCH_PARTITION=rna.
#SBATCH --output=fqxv-cell-%A_%a.out    # throughput headroom while still packing.
## #SBATCH --account=amc-general
# Post-processing (round-trip digest + Δqual) is now all compiled fqdigest — the
# per-variant single-thread determinism recompress and the interpreted per-byte
# awk are both gone — so cells finish in roughly compress+decompress time. 2h is
# generous headroom (the slow HiFi/CoLoRd cells are the long pole, not verify).
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
