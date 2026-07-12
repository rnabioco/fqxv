#!/usr/bin/env bash
#SBATCH --job-name=fqxv-bench
#SBATCH --comment=fqxv-bench
#SBATCH --partition=amilan          # AMD EPYC 7713P; renamed 'acpu' after 2026-08-05
#SBATCH --qos=normal                # renamed 'cpu-normal' after 2026-08-05
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
## #SBATCH --account=amc-general    # uncomment if your default association is unset
#
# Build the release binary AND run the full baseline benchmark on one compute
# node (never the login node): compile with 32 cores, then benchmark with 32
# threads. Results land in $SCRATCH/fqxv/results/results.tsv.
set -euo pipefail

export PATH="$HOME/.pixi/bin:$HOME/.cargo/bin:$PATH"
# In a batch job ${BASH_SOURCE[0]} is Slurm's spooled copy, not the repo — use
# the submit dir (this script lives in and is submitted from bench/).
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
ROOT="$(cd "$HERE/.." && pwd)"

export FQXV_THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
export RAYON_NUM_THREADS="$FQXV_THREADS"

echo "==> node $(hostname), $FQXV_THREADS threads, $(date)"
echo "==> building release binary on the compute node"
( cd "$ROOT" && cargo build --release -j "$FQXV_THREADS" )

echo "==> running benchmark harness"
pixi run --manifest-path "$HERE/pixi.toml" bash "$HERE/run_bench.sh"

echo "==> rendering table"
pixi run --manifest-path "$HERE/pixi.toml" python "$HERE/report.py"
echo "==> done $(date)"
