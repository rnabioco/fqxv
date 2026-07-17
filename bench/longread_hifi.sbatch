#!/usr/bin/env bash
#SBATCH --job-name=fqxv-hifi
#SBATCH --comment=fqxv-hifi
#SBATCH --partition=rna
#SBATCH --qos=normal
#SBATCH --account=rbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=44
#SBATCH --mem=128G
#SBATCH --time=03:00:00
#SBATCH --output=%x-%j.out
#
# PacBio HiFi benchmark on Bodhi: rebuild the worktree fqxv (with the new
# ONT/HiFi bins), subsample the deep 44 GB HiFi run to ~120k reads, then run the
# harness comparing fqxv (incl. fqxv-binhifi) against CoLoRd lossless + lossy.
set -euo pipefail
export PATH="$HOME/.pixi/bin:$HOME/.cargo/bin:$PATH"
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
ROOT="$(cd "$HERE/.." && pwd)"
THREADS="${SLURM_CPUS_PER_TASK:-$(nproc)}"
export FQXV_THREADS="$THREADS" RAYON_NUM_THREADS="$THREADS"

DATA_DIR="${FQXV_DATA_DIR:-$HOME/scratch/fqxv/data}"
SRC="$DATA_DIR/SRR11434954.fastq"
SUB="$DATA_DIR/ecoli_hifi.fastq"
N_READS="${HIFI_SUBSAMPLE_READS:-120000}"

echo "==> node $(hostname), $THREADS threads, $(date)"
echo "==> build release fqxv (with ONT/HiFi bins), HEAD $(git -C "$ROOT" rev-parse --short HEAD)"
( cd "$ROOT" && cargo build --release -p fqxv-cli -j "$THREADS" )

if [[ ! -s "$SUB" ]]; then
  echo "==> subsample $SRC -> $SUB (first $N_READS reads)"
  pixi run --manifest-path "$HERE/pixi.toml" seqkit head -n "$N_READS" "$SRC" -o "$SUB"
fi
echo "==> HiFi subset stats:"
pixi run --manifest-path "$HERE/pixi.toml" seqkit stats -a "$SUB" | tail -2

export FQXV_BIN="$ROOT/target/release/fqxv"
export FQXV_DATA_DIR="$DATA_DIR"
export FQXV_RESULTS_DIR="${FQXV_RESULTS_DIR:-$HOME/scratch/fqxv/results-hifi}"
export FQXV_ONLY_DATASET="ecoli_hifi"
export FQXV_INPUT="r1"
export FQXV_TOOLS="fqxv fqxv9 fqxv-max fqxv-binhifi fqxv-binont gzip zstd19 xz9 colord colord-lossy"

rm -rf "$FQXV_RESULTS_DIR"; mkdir -p "$FQXV_RESULTS_DIR"
echo "==> run_bench.sh (dataset=$FQXV_ONLY_DATASET, tools=$FQXV_TOOLS)"
pixi run --manifest-path "$HERE/pixi.toml" bash "$HERE/run_bench.sh"
echo "==> report"
pixi run --manifest-path "$HERE/pixi.toml" python "$HERE/report.py" | tee "$FQXV_RESULTS_DIR/report.txt"
echo "==> done $(date)"
