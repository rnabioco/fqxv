#!/usr/bin/env bash
#SBATCH --job-name=fqxv-longread
#SBATCH --comment=fqxv-longread
#SBATCH --partition=rna
#SBATCH --qos=normal
#SBATCH --account=rbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=44
#SBATCH --mem=96G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
#
# Long-read (ONT) benchmark on Bodhi: build the *worktree's* fqxv, then run the
# harness on the ecoli_ont dataset only, comparing fqxv against CoLoRd (the
# long-read SOTA) plus gzip/zstd/xz baselines. Ratios are thread-independent;
# throughput uses all 44 cores. Submitted from the worktree bench/ dir.
set -euo pipefail

export PATH="$HOME/.pixi/bin:$HOME/.cargo/bin:$PATH"
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
ROOT="$(cd "$HERE/.." && pwd)"

THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
export FQXV_THREADS="$THREADS"
export RAYON_NUM_THREADS="$THREADS"

echo "==> node $(hostname), $THREADS threads, worktree $ROOT, $(date)"
echo "==> git HEAD: $(git -C "$ROOT" rev-parse --short HEAD) ($(git -C "$ROOT" rev-parse --abbrev-ref HEAD))"

echo "==> build release fqxv in the worktree"
( cd "$ROOT" && cargo build --release -p fqxv-cli -j "$THREADS" )

# Benchmark the worktree binary (NOT a stale $SCRATCH build).
export FQXV_BIN="$ROOT/target/release/fqxv"
export FQXV_DATA_DIR="${FQXV_DATA_DIR:-$HOME/scratch/fqxv/data}"
export FQXV_RESULTS_DIR="${FQXV_RESULTS_DIR:-$HOME/scratch/fqxv/results-longread}"
export FQXV_ONLY_DATASET="ecoli_ont"
export FQXV_INPUT="r1"
# Focused set: fqxv points + CoLoRd (lossless -q org) + general baselines.
# fqz_comp/fqzcomp5 are included to confirm whether they parse long reads.
export FQXV_TOOLS="fqxv fqxv9 fqxv-max fqxv-bin4 fqxv-bin2 gzip zstd19 xz9 fqz_comp fqzcomp5 colord"

rm -rf "$FQXV_RESULTS_DIR"; mkdir -p "$FQXV_RESULTS_DIR"
echo "==> run_bench.sh (dataset=$FQXV_ONLY_DATASET, tools=$FQXV_TOOLS)"
pixi run --manifest-path "$HERE/pixi.toml" bash "$HERE/run_bench.sh"

echo "==> report"
pixi run --manifest-path "$HERE/pixi.toml" python "$HERE/report.py" | tee "$FQXV_RESULTS_DIR/report.txt"
echo "==> done $(date)"
