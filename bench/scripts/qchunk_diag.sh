#!/usr/bin/env bash
# Chunked-quality diagnostic sweep (design/parallel-decode, docs/design/parallel-decode.md).
#
# Measures, on real archives, what a parallel-first quality stream would COST:
# the encoder-side FQXV_DIAG_QCHUNK probe re-codes every quality block under
# {reset, warm-clone, warm-frozen} x K in {4,8,16,32} x warmup in {total/K, 8 MiB}
# and prints one machine-parsable stderr line per cell; FQXV_DIAG_TILECHUNK does
# the same for the ONT tiler's chunk-confined reference selection. Neither flag
# changes the emitted archive — this script proves that with cmp on every
# dataset (a hard repo invariant, so a probe regression fails the sweep).
#
# Alongside the ratio sweep it measures the DECODE share the chunks would
# parallelize: full decode vs --fasta (no quality stream) at 1 thread, 3 reps,
# median, output piped to `wc -c` (never /dev/null — tools fast-path a null
# sink), plus 8/16-thread full-decode anchors. quality_share = 1 - t_fasta/t_full.
#
# Datasets (all single-stream):
#   DRR205413    ONT MinION      staged by fetch.sh
#   ecoli_hifi   PacBio HiFi     staged by stage_hifi.sh (SRR11434954 head)
#   novaseq_4m   NovaSeq control first 4M reads of DRR174812_1 (MODE_POS path)
#
# Meant to run inside a dedicated sbatch (slurm/qchunk_diag.sbatch), never in a
# shared build allocation. Env knobs:
#   FQXV_DATA_DIR       staged FASTQ dir     (default $SCRATCH/fqxv/data)
#   FQXV_RESULTS_DIR    output dir           (default $SCRATCH/fqxv/results)
#   FQXV_BIN            fqxv binary          (default <repo>/target/release/fqxv)
#   FQXV_QCHUNK_DATASETS  subset of "DRR205413 ecoli_hifi novaseq_4m"
#   FQXV_REPS           timed repetitions    (default 3; median reported)
#   FQXV_REPREP=1       force re-compress of the archives
#
# Outputs in $FQXV_RESULTS_DIR:
#   qchunk_diag.tsv        per-block cells + per-dataset TOTAL rows (payload- and
#                          archive-relative deltas)
#   qchunk_tilechunk.tsv   ONT tiler confinement cells + TOTAL rows
#   qchunk_decode_share.tsv decode timings + quality_share rows
# Per-dataset raw logs and archives stay under $FQXV_RESULTS_DIR/qchunk_diag/
# so a partial run resumes without recompressing finished datasets.
set -euo pipefail

DATA_DIR="${FQXV_DATA_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/data}"
RESULTS_DIR="${FQXV_RESULTS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/results}"
DATASETS="${FQXV_QCHUNK_DATASETS:-DRR205413 ecoli_hifi novaseq_4m}"
REPS="${FQXV_REPS:-3}"
REPREP="${FQXV_REPREP:-0}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FQXV_BIN="${FQXV_BIN:-${CARGO_TARGET_DIR:-$ROOT/target}/release/fqxv}"
NPROC="$(nproc)"

WORK="$RESULTS_DIR/qchunk_diag"
mkdir -p "$WORK"
QC_TSV="$RESULTS_DIR/qchunk_diag.tsv"
TILE_TSV="$RESULTS_DIR/qchunk_tilechunk.tsv"
SHARE_TSV="$RESULTS_DIR/qchunk_decode_share.tsv"

[[ -x "$FQXV_BIN" ]] || {
  echo "error: fqxv binary not found at $FQXV_BIN (cargo build --release -p fqxv-cli)" >&2
  exit 1
}

ver() { "$@" 2>&1 | head -n1; }

meta_header() {
  echo "# $1 — $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "# node: $(hostname)  nproc: $NPROC  slurm_job: ${SLURM_JOB_ID:-none}"
  echo "# cpu: $(awk -F': ' '/model name/{print $2; exit}' /proc/cpuinfo)"
  echo "# commit: $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
  echo "# fqxv: $(ver "$FQXV_BIN" --version)  [$FQXV_BIN]"
  echo "# datasets: $DATASETS  reps: $REPS"
}

{
  meta_header "fqxv chunked-quality encode-cost sweep (FQXV_DIAG_QCHUNK)"
  printf 'dataset\tblk\tmode\treads\tbases\tk\tbaseline\tvariant\tK\twarmup\twarm_bases\twarm_bytes\tchunk_bytes\thdr\ttotal\tdelta_pct\tarch_delta_pct\n'
} > "$QC_TSV"
{
  meta_header "fqxv ONT tile chunk-confinement sweep (FQXV_DIAG_TILECHUNK)"
  printf 'dataset\tblk\treads\tbases\tbaseline\tK\tconfined\tdelta_pct\n'
} > "$TILE_TSV"
{
  meta_header "fqxv decode share: full vs --fasta (quality share of decode)"
  printf 'dataset\tmode\tthreads\tseconds\tmb_per_s\tbytes_out\tcomp_bytes\tquality_share\treps\n'
} > "$SHARE_TSV"

fresh() {
  local out="$1"; shift
  [[ "$REPREP" == 1 ]] && return 1
  [[ -s "$out" ]] || return 1
  local d
  for d in "$@"; do [[ "$out" -nt "$d" ]] || return 1; done
  return 0
}

# median <newline-separated floats> -> stdout
median() {
  sort -g | awk '{v[NR]=$1} END{m=int((NR+1)/2); if (NR%2) print v[m]; else printf "%.3f\n", (v[m]+v[m+1])/2}'
}

# time_pipe <cmd> -> CELL_SECS (median wall seconds of REPS runs, 1 warmup).
time_pipe() {
  local cmd="$1" r t0 t1 times=""
  bash -c "set -o pipefail; $cmd" > /dev/null # warmup (also correctness pass)
  for ((r = 0; r < REPS; r++)); do
    t0="$EPOCHREALTIME"
    bash -c "set -o pipefail; $cmd" > /dev/null
    t1="$EPOCHREALTIME"
    times="$times$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b-a}')"$'\n'
  done
  CELL_SECS="$(printf '%s' "$times" | median)"
}

# Parse "[diag qchunk] k1=v1 k2=v2 ..." stderr lines into TSV rows + TOTAL rows.
# Aggregation keys on (mode, variant, K, warmup): per-dataset payload-relative
# delta = sum(total-baseline)/sum(baseline); archive-relative = /archive bytes.
emit_qchunk() { # <dataset> <log> <archive_bytes>
  awk -v ds="$1" -v arch="$3" '
    /^\[diag qchunk\]/ {
      delete f
      for (i = 3; i <= NF; i++) { split($i, kv, "="); f[kv[1]] = kv[2] }
      key = f["mode"] SUBSEP f["variant"] SUBSEP f["K"] SUBSEP f["warmup"]
      base[key] += f["baseline"]; tot[key] += f["total"]
      if (!(key in seen)) { seen[key] = 1; order[++n] = key }
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.4f\n",
        ds, f["blk"], f["mode"], f["reads"], f["bases"], f["k"], f["baseline"],
        f["variant"], f["K"], f["warmup"], f["warm_bases"], f["warm_bytes"],
        f["chunk_bytes"], f["hdr"], f["total"], f["delta_pct"],
        (f["total"] - f["baseline"]) * 100.0 / arch
    }
    END {
      for (j = 1; j <= n; j++) {
        key = order[j]; split(key, p, SUBSEP)
        printf "%s\tTOTAL\t%s\t-\t-\t-\t%d\t%s\t%s\t%s\t-\t-\t-\t-\t%d\t%.4f\t%.4f\n",
          ds, p[1], base[key], p[2], p[3], p[4], tot[key],
          (tot[key] - base[key]) * 100.0 / base[key],
          (tot[key] - base[key]) * 100.0 / arch
      }
    }' "$2"
}

emit_tilechunk() { # <dataset> <log>
  awk -v ds="$1" '
    /^\[diag tilechunk\]/ {
      delete f
      for (i = 3; i <= NF; i++) { split($i, kv, "="); f[kv[1]] = kv[2] }
      base[f["K"]] += f["baseline"]; tot[f["K"]] += f["confined"]
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n",
        ds, f["blk"], f["reads"], f["bases"], f["baseline"], f["K"],
        f["confined"], f["delta_pct"]
    }
    END {
      for (k in base)
        printf "%s\tTOTAL\t-\t-\t%d\t%s\t%d\t%.4f\n",
          ds, base[k], k, tot[k], (tot[k] - base[k]) * 100.0 / base[k]
    }' "$2"
}

for ds in $DATASETS; do
  case "$ds" in
    DRR205413) in="$DATA_DIR/DRR205413.fastq" tile=1 ;;
    ecoli_hifi) in="$DATA_DIR/ecoli_hifi.fastq" tile=0 ;;
    novaseq_4m)
      # MODE_POS control: first 4M reads of the NovaSeq R1 (the benchdata-style
      # subsample); head is record-safe on 4-line FASTQ.
      in="$WORK/novaseq_4m.fastq" tile=0
      src="$DATA_DIR/DRR174812_1.fastq"
      [[ -s "$src" ]] || { echo "== $ds: $src not staged — skipping" >&2; continue; }
      fresh "$in" "$src" || { echo "== $ds: head -4M reads"; head -n 16000000 "$src" > "$in"; }
      ;;
    *) echo "error: unknown dataset '$ds'" >&2; exit 1 ;;
  esac
  [[ -s "$in" ]] || { echo "== $ds: $in not staged — skipping" >&2; continue; }

  arc="$WORK/$ds.fqxv"
  fresh "$arc" "$in" "$FQXV_BIN" || {
    echo "== $ds: baseline compress"
    "$FQXV_BIN" compress "$in" -o "$arc" --force --quiet --threads "$NPROC"
  }
  arch_bytes="$(stat -Lc %s "$arc")"

  # --- qchunk sweep (and the probe no-op proof: archive must be byte-identical)
  qlog="$WORK/$ds.qchunk.log"
  fresh "$qlog" "$arc" || {
    echo "== $ds: FQXV_DIAG_QCHUNK compress (the ~21x quality-encode sweep)"
    FQXV_DIAG_QCHUNK=1 "$FQXV_BIN" compress "$in" -o "$arc.diag" --force --quiet \
      --threads "$NPROC" 2> "$qlog.tmp"
    cmp "$arc" "$arc.diag" || {
      echo "!! $ds: FQXV_DIAG_QCHUNK changed the archive — probe is not a no-op" >&2
      exit 1
    }
    rm -f "$arc.diag"
    mv "$qlog.tmp" "$qlog"
  }
  emit_qchunk "$ds" "$qlog" "$arch_bytes" >> "$QC_TSV"

  # --- tilechunk sweep (ONT only; separate run so its cost stays identifiable)
  if [[ "$tile" == 1 ]]; then
    tlog="$WORK/$ds.tilechunk.log"
    fresh "$tlog" "$arc" || {
      echo "== $ds: FQXV_DIAG_TILECHUNK compress (2 extra tile encodes/block)"
      FQXV_DIAG_TILECHUNK=1 "$FQXV_BIN" compress "$in" -o "$arc.tdiag" --force --quiet \
        --threads "$NPROC" 2> "$tlog.tmp"
      cmp "$arc" "$arc.tdiag" || {
        echo "!! $ds: FQXV_DIAG_TILECHUNK changed the archive — probe is not a no-op" >&2
        exit 1
      }
      rm -f "$arc.tdiag"
      mv "$tlog.tmp" "$tlog"
    }
    emit_tilechunk "$ds" "$tlog" >> "$TILE_TSV"
  fi

  # --- decode share: full vs --fasta at 1 thread (+ 8/16-thread anchors)
  echo "== $ds: decode share (full vs --fasta)"
  declare -A secs
  for cell in full:1 fasta:1 full:8 full:16; do
    mode="${cell%%:*}" n="${cell##*:}"
    [[ "$n" -le "$NPROC" ]] || continue
    flag=""
    [[ "$mode" == fasta ]] && flag="--fasta"
    cmd="'$FQXV_BIN' decompress '$arc' -Z $flag --quiet --threads $n | wc -c"
    bytes="$(bash -c "set -o pipefail; $cmd")"
    time_pipe "$cmd"
    secs["$cell"]="$CELL_SECS"
    mbps="$(awk -v b="$bytes" -v s="$CELL_SECS" 'BEGIN{printf "%.1f", (s > 0) ? b/s/1e6 : 0}')"
    share="-"
    if [[ "$cell" == "fasta:1" && -n "${secs[full:1]:-}" ]]; then
      share="$(awk -v f="${secs[full:1]}" -v a="$CELL_SECS" 'BEGIN{printf "%.4f", 1 - a/f}')"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
      "$ds" "$mode" "$n" "$CELL_SECS" "$mbps" "$bytes" "$arch_bytes" "$share" "$REPS" >> "$SHARE_TSV"
    echo "   $ds $mode t=$n ${CELL_SECS}s ${mbps} MB/s share=$share"
  done
  unset secs
done

echo "==> $QC_TSV"
echo "==> $TILE_TSV"
echo "==> $SHARE_TSV"
