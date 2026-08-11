#!/usr/bin/env bash
#SBATCH --job-name=fqxv-qchunk-diag
#SBATCH --comment=fqxv-qchunk-diag
#SBATCH --partition=rna
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=44
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --output=%x-%j.out
#
# Chunked-quality diagnostic sweep (see scripts/qchunk_diag.sh): builds the
# CLI, stages the HiFi dataset if missing, then runs the FQXV_DIAG_QCHUNK /
# FQXV_DIAG_TILECHUNK sweeps + the decode-share measurement. Results land in
# $SCRATCH/fqxv/results/qchunk_{diag,tilechunk,decode_share}.tsv. Re-runnable:
# archives and diag logs are mtime-checked per dataset (FQXV_REPREP=1 to force),
# so a partial run resumes where it stopped.
#
# The qchunk sweep costs ~21x the quality-encode work per block and the
# tilechunk sweep two extra tile encodes per ONT block — this belongs in its own
# allocation, never a shared build session.
#
# From a worktree (bench/.pixi is gitignored, so worktrees lack it), point
# FQXV_PIXI_BIN at a checkout that has the env, for stage_hifi's sracha/seqkit:
#   FQXV_PIXI_BIN=~/devel/rnabioco/fqxv/bench/.pixi/envs/default/bin \
#     sbatch bench/slurm/qchunk_diag.sbatch
set -euo pipefail

BENCH="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# Tolerate submission from the repo root, bench/, or bench/slurm/.
for _c in "$BENCH" "$BENCH/.." "$BENCH/bench"; do
  [ -d "$_c/scripts" ] && [ -d "$_c/slurm" ] && { BENCH="$(cd "$_c" && pwd)"; break; }
done
ROOT="$(cd "$BENCH/.." && pwd)"
DATA_DIR="${FQXV_DATA_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/data}"

echo "==> node $(hostname), ${SLURM_CPUS_PER_TASK:-$(nproc)} cpus, $(date)"
echo "==> repo $ROOT @ $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"

# Build the binary the sweep measures (skip with FQXV_BIN pointing elsewhere).
if [[ -z "${FQXV_BIN:-}" ]]; then
  (cd "$ROOT" && cargo build --release -p fqxv-cli)
fi

# Stage the HiFi dataset unless it is already there (stage_hifi.sh is
# idempotent; FQXV_PIXI_BIN passes through for worktree submissions).
if [[ ! -s "$DATA_DIR/ecoli_hifi.fastq" ]]; then
  bash "$BENCH/scripts/stage_hifi.sh"
fi

bash "$BENCH/scripts/qchunk_diag.sh"
