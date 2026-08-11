#!/usr/bin/env bash
#SBATCH --job-name=fqxv-decode-scaling
#SBATCH --comment=fqxv-decode-scaling
#SBATCH --partition=rna             # Bodhi (Intel, 88+ cores / 754G). On Alpine use
#SBATCH --qos=normal                # `SBATCH_PARTITION=amilan sbatch ...`, or source
#                                   # scripts/cluster.sh, which detects the cluster.
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=%x-%j.out
#
# Decode thread-scaling benchmark (see scripts/decode_scaling.sh):
#   sbatch slurm/decode_scaling.sbatch
# Results land in $SCRATCH/fqxv/results/decode_scaling.tsv. Re-runnable: the
# compressed inputs are synthesized once and reused (FQXV_REPREP=1 to rebuild);
# knobs (FQXV_DECODE_DATASETS, FQXV_THREAD_STEPS, FQXV_RESULTS_DIR, ...) pass
# through the environment. Deliberately NOT --exclusive: a dedicated 32-core
# allocation is the measurement environment the docs describe.
set -euo pipefail

export PATH="$HOME/.pixi/bin:$PATH"
BENCH="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)}"
# Tolerate submission from the repo root, bench/, or bench/slurm/.
for _c in "$BENCH" "$BENCH/.." "$BENCH/bench"; do
  [ -d "$_c/scripts" ] && [ -d "$_c/panels" ] && { BENCH="$(cd "$_c" && pwd)"; break; }
done
ROOT="$(cd "$BENCH/.." && pwd)"

echo "==> node $(hostname), ${SLURM_CPUS_PER_TASK:-$(nproc)} cpus, $(date)"
pixi run -e bench --manifest-path "$ROOT/pixi.toml" bash "$BENCH/scripts/decode_scaling.sh"
