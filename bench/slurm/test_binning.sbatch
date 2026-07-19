#!/usr/bin/env bash
#SBATCH --job-name=fqxv-testbin
#SBATCH --comment=fqxv-testbin
#SBATCH --partition=rna
#SBATCH --qos=normal
#SBATCH --account=rbi
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=00:30:00
#SBATCH --output=%x-%j.out
#
# Compile + test the ONT/HiFi quality-binning additions on a real compute node.
set -euo pipefail
export PATH="$HOME/.cargo/bin:$HOME/.pixi/bin:$PATH"
BENCH="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# Tolerate submission from the repo root, bench/, or bench/slurm/.
for _c in "$BENCH" "$BENCH/.." "$BENCH/bench"; do
  [ -d "$_c/scripts" ] && [ -d "$_c/panels" ] && { BENCH="$(cd "$_c" && pwd)"; break; }
done
ROOT="$(cd "$BENCH/.." && pwd)"
cd "$ROOT"
echo "==> node $(hostname), $(nproc) cores, HEAD $(git rev-parse --short HEAD)"
echo "==> cargo test -p fqxv-fqzcomp (binning tables)"
cargo test -p fqxv-fqzcomp 2>&1 | tail -30
echo "==> cargo test -p fqxv (container: binning_tag, apply_binning)"
cargo test -p fqxv --lib 2>&1 | tail -15
echo "==> done $(date)"
