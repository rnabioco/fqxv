#!/usr/bin/env bash
#SBATCH --job-name=fqxv-fmtclippy
#SBATCH --comment=fqxv-fmtclippy
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
# Match CI: apply rustfmt, then clippy with -Dwarnings (+ the bench feature CI
# passes). Run before committing the ONT/HiFi binning work.
set -euo pipefail
export PATH="$HOME/.cargo/bin:$PATH"
ROOT="${SLURM_SUBMIT_DIR:-$(pwd)}/.."
cd "$ROOT"
echo "==> node $(hostname), $(nproc) cores"
echo "==> cargo fmt --all (applying)"
cargo fmt --all
echo "==> fmt diff (files rustfmt changed):"
git diff --stat -- '*.rs' | tail -5
echo "==> cargo fmt --all --check"
cargo fmt --all --check && echo "FMT: clean"
echo "==> cargo clippy (CI flags)"
RUSTFLAGS=-Dwarnings cargo clippy --workspace --all-targets --features fqxv-rans/bench 2>&1 | tail -20
echo "CLIPPY exit: $?"
echo "==> done $(date)"
