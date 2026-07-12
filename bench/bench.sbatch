#!/usr/bin/env bash
#SBATCH --job-name=fqxv-bench
#SBATCH --partition=amilan          # AMD EPYC 7713P; renamed 'acpu' after 2026-08-05
#SBATCH --qos=normal                # renamed 'cpu-normal' after 2026-08-05
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --exclusive                 # whole node: no co-tenant memory-bandwidth noise
#SBATCH --mem=0                     # all node memory
#SBATCH --time=04:00:00
#SBATCH --output=%x-%j.out
## #SBATCH --account=amc-general    # uncomment if your default association is unset
#
# Run the full baseline benchmark on one amilan node:
#   sbatch bench.sbatch
# Results land in $SCRATCH/fqxv/results/results.tsv; render with report.py.
set -euo pipefail

export PATH="$HOME/.pixi/bin:$PATH"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

export FQXV_THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
export RAYON_NUM_THREADS="$FQXV_THREADS"

echo "==> node $(hostname), $FQXV_THREADS threads, $(date)"
pixi run --manifest-path "$HERE/pixi.toml" bash "$HERE/run_bench.sh"
echo "==> rendering table"
pixi run --manifest-path "$HERE/pixi.toml" python "$HERE/report.py"
