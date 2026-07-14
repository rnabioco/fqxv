#!/usr/bin/env bash
#SBATCH --job-name=fqxv-bam-id
#SBATCH --comment=fqxv-bam-id
#SBATCH --partition=rna
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=%x-%j.out
#
# Alignment-level round-trip proof (see bam_identity.sh). Run from bench/:
#   sbatch bam_identity.sbatch                 # default dataset (ecoli_miseq)
#   sbatch bam_identity.sbatch ecoli_miseq     # explicit label
set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
export FQXV_THREADS="${SLURM_CPUS_PER_TASK:-32}"
exec pixi run bash bam_identity.sh "$@"
