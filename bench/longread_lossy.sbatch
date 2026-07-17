#!/usr/bin/env bash
#SBATCH --job-name=fqxv-lrlossy
#SBATCH --comment=fqxv-lrlossy
#SBATCH --partition=rna
#SBATCH --qos=normal
#SBATCH --account=rbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=44
#SBATCH --mem=128G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
#
# Lossy long-read quality comparison on BOTH ONT and HiFi: fqxv platform bins
# (binont/binhifi) vs CoLoRd's default lossy mode, with lossless fqxv9/colord as
# anchors. Reuses the release binary built by the HiFi job (rebuilds if absent).
set -euo pipefail
export PATH="$HOME/.pixi/bin:$HOME/.cargo/bin:$PATH"
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
ROOT="$(cd "$HERE/.." && pwd)"
THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
export FQXV_THREADS="$THREADS" RAYON_NUM_THREADS="$THREADS"
export FQXV_BIN="$ROOT/target/release/fqxv"
[[ -x "$FQXV_BIN" ]] || ( cd "$ROOT" && cargo build --release -p fqxv-cli -j "$THREADS" )

export FQXV_DATA_DIR="$HOME/scratch/fqxv/data"
export FQXV_INPUT="r1"
export FQXV_TOOLS="fqxv9 colord fqxv-binont fqxv-binhifi colord-lossy"

for ds in ecoli_ont ecoli_hifi; do
  echo "############ $ds ############"
  export FQXV_ONLY_DATASET="$ds"
  export FQXV_RESULTS_DIR="$HOME/scratch/fqxv/results-lossy-$ds"
  rm -rf "$FQXV_RESULTS_DIR"; mkdir -p "$FQXV_RESULTS_DIR"
  pixi run --manifest-path "$HERE/pixi.toml" bash "$HERE/run_bench.sh"
  pixi run --manifest-path "$HERE/pixi.toml" python "$HERE/report.py" | tee "$FQXV_RESULTS_DIR/report.txt"
done
echo "==> done $(date)"
