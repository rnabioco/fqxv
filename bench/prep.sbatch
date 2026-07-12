#!/usr/bin/env bash
#SBATCH --job-name=fqxv-prep
#SBATCH --comment=fqxv-prep
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=90G
#SBATCH --time=01:00:00
#SBATCH --output=%x-%j.out
## #SBATCH --account=amc-general
#
# Phase 1 of the parallel harness: build the release binary once, and compute
# the per-dataset meta (record/base counts) + the expensive input-content
# digests once, into $RESULTS_DIR/prep/<label>.env. The array cells then reuse
# these instead of each re-sorting multi-GB FASTQ. Submitted by submit_parallel.sh.
set -euo pipefail

export PATH="$HOME/.pixi/bin:$HOME/.cargo/bin:$PATH"
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
ROOT="$(cd "$HERE/.." && pwd)"
export FQXV_THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

echo "==> prep on $(hostname), $FQXV_THREADS threads, $(date)"
echo "==> building release binary"
( cd "$ROOT" && cargo build --release -j "$FQXV_THREADS" )

echo "==> computing per-dataset meta + input digests"
export FQXV_PREP_ONLY=1
pixi run --manifest-path "$HERE/pixi.toml" bash "$HERE/run_bench.sh"
echo "==> prep done $(date)"
