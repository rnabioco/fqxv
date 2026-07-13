#!/usr/bin/env bash
#SBATCH --job-name=fqxv-merge
#SBATCH --comment=fqxv-merge
#SBATCH --partition=amilan
#SBATCH --qos=normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
#SBATCH --output=%x-%j.out
## #SBATCH --account=amc-general
#
# Phase 3 of the parallel harness: concatenate the per-cell part files into one
# results.tsv (with a header) and render the report. Runs afterany the array, so
# a single failed cell doesn't block the table. meta.tsv came from the prep phase.
set -euo pipefail

export PATH="$HOME/.pixi/bin:$PATH"
HERE="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
RESULTS_DIR="${FQXV_RESULTS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/results}"
R="$RESULTS_DIR/results.tsv"

echo -e "dataset\ttool\torig_bytes\tcomp_bytes\tratio\tc_secs\td_secs\tc_rss_kb\td_rss_kb\tnames_bytes\tseq_bytes\tqual_bytes\trt_ok\tdeterministic\tqual_mae\tqual_rmse\tqual_pct_changed" > "$R"
cat "$RESULTS_DIR/parts"/results.*.tsv >> "$R" 2>/dev/null || true
echo "==> merged $(( $(wc -l < "$R") - 1 )) result rows into $R"

echo "==> rendering table"
pixi run --manifest-path "$HERE/pixi.toml" python "$HERE/report.py"
echo "==> done $(date)"
