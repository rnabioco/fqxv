#!/usr/bin/env bash
# fqxv vs native .sra — is an .fqxv worth keeping instead of the archive the data
# actually shipped in? `.sra` is itself a compressed columnar store (2-bit-packed
# clean bases, a quality model, spot/name columns), so for anyone pulling from
# SRA/ENA the honest question is "fqxv vs the .sra I'd otherwise keep", not just
# vs gzip/zstd. This records, per run:
#   * .sra size  — the native NCBI archive, from `sracha info --format tsv`
#                  (metadata ONLY; the multi-GB .sra is never downloaded).
#   * fqxv size  — `fqxv compress` of the same run's FASTQ, at several operating
#                  points: --max (lossless, best ratio) and the lossy quality
#                  bins bin8/bin4/bin2.
# and reports fqxv/.sra per point (< 1.0 = fqxv is smaller). Emits TSV to
# $RESULTS_DIR/sra_compare.tsv, rendered by report.py.
#
# Two phases, like run_bench.sh: the .sra metadata lookup needs the network (run
# on the login node — it is tiny and IO-bound); the fqxv compress needs a compute
# node (never the login node). The `sizes` phase caches .sra sizes so the compute
# phase needs no network:
#   pixi run -e bench bash bench/scripts/sra_compare.sh sizes       # login node: cache .sra sizes
#   srun ... pixi run bash sra_compare.sh run # compute node: compress + join
# The default `all` does both in one go (fine when the node can reach the network).
#
# Env knobs (shared with run_bench.sh where they overlap):
#   FQXV_DATA_DIR      input FASTQ dir           (default $SCRATCH/fqxv/data)
#   FQXV_RESULTS_DIR   output dir                (default $SCRATCH/fqxv/results)
#   FQXV_THREADS       threads                   (default nproc)
#   FQXV_BIN           fqxv release binary
#   SRA_PANEL          panel TSV                 (default bench/sra_panel.tsv)
#   SRA_POINTS         fqxv points               (default "max bin8 bin4 bin2")
#   SRA_INCLUDE_LARGE  include size_class=large rows (default 0)
set -euo pipefail

MODE="${1:-all}"   # sizes | run | all
case "$MODE" in sizes|run|all) ;; *) echo "usage: sra_compare.sh [sizes|run|all]" >&2; exit 2 ;; esac

DATA_DIR="${FQXV_DATA_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/data}"
RESULTS_DIR="${FQXV_RESULTS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/results}"
THREADS="${FQXV_THREADS:-$(nproc)}"
POINTS="${SRA_POINTS:-max bin8 bin4 bin2}"
INCLUDE_LARGE="${SRA_INCLUDE_LARGE:-0}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PANEL="${SRA_PANEL:-$HERE/../panels/sra_panel.tsv}"
# Same resolution as run_bench.sh: cargo honors CARGO_TARGET_DIR (set to $SCRATCH
# on this HPC), so resolve where cargo actually wrote the binary, not ROOT/target.
FQXV_BIN="${FQXV_BIN:-${CARGO_TARGET_DIR:-$(cd "$HERE/.." && pwd)/target}/release/fqxv}"

mkdir -p "$RESULTS_DIR"
SIZES="$RESULTS_DIR/sra_sizes.tsv"       # cached .sra metadata (sizes phase)
OUT="$RESULTS_DIR/sra_compare.tsv"       # the comparison table (run phase)
WORK="$RESULTS_DIR/work"; mkdir -p "$WORK"

# GNU time for wall + peak RSS; fall back to bash timing (RSS unknown = -1).
GNU_TIME=""
for c in /usr/bin/time "$(command -v time || true)"; do
  if [[ -x "$c" ]] && "$c" -f '%e %M' true >/dev/null 2>&1; then GNU_TIME="$c"; break; fi
done
measure() {  # measure CMD... -> MEAS_SECS, MEAS_RSS_KB, MEAS_RC
  local tf; tf="$(mktemp)"; MEAS_RC=0
  if [[ -n "$GNU_TIME" ]]; then
    "$GNU_TIME" -o "$tf" -f '%e %M' "$@" || MEAS_RC=$?
    read -r MEAS_SECS MEAS_RSS_KB < <(tail -n1 "$tf") || { MEAS_SECS=-1; MEAS_RSS_KB=-1; }
  else
    local t0 t1; t0="$EPOCHREALTIME"; { "$@" || MEAS_RC=$?; }; t1="$EPOCHREALTIME"
    MEAS_SECS="$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}')"; MEAS_RSS_KB=-1
  fi
  rm -f "$tf"
}

# Rows of the panel, minus comments/blanks and (unless opted in) large runs.
panel_rows() {
  grep -v '^#' "$PANEL" | awk 'NF' | while IFS=$'\t' read -r acc plat regime cls _; do
    [[ "$cls" == large && "$INCLUDE_LARGE" != 1 ]] && continue
    printf '%s\t%s\t%s\t%s\n' "$acc" "$plat" "$regime" "${cls:-small}"
  done
}

# ---- sizes phase: cache native .sra sizes from sracha (network, login node) ----
# `sracha info --format tsv` columns: accession archive_type layout nreads spots
# size_bytes platform md5. We keep size_bytes (the native .sra archive size) plus
# layout/spots. The awk `$1==acc` guard extracts exactly the data row, ignoring
# any activation/log noise on stdout.
do_sizes() {
  command -v sracha >/dev/null 2>&1 || { echo "sracha not found (run under: pixi run -e bench bash bench/scripts/sra_compare.sh sizes)" >&2; exit 1; }
  echo -e "accession\tplatform\tregime\tlayout\tspots\tsra_bytes" > "$SIZES"
  panel_rows | while IFS=$'\t' read -r acc plat regime cls; do
    local line
    line="$(sracha info "$acc" --format tsv 2>/dev/null | awk -F'\t' -v a="$acc" '$1==a{print; exit}')" || true
    if [[ -z "$line" ]]; then
      echo "  [warn] no .sra metadata for $acc (skipped)"; continue
    fi
    local layout spots sra_bytes
    layout="$(cut -f3 <<<"$line")"; spots="$(cut -f5 <<<"$line")"; sra_bytes="$(cut -f6 <<<"$line")"
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$acc" "$plat" "$regime" "$layout" "$spots" "$sra_bytes" >> "$SIZES"
    echo "  $acc  $plat  .sra=$(numfmt --to=iec "$sra_bytes")  spots=$spots  layout=$layout"
  done
  echo "==> wrote $SIZES"
}

# ---- run phase: fqxv-compress each run, join against cached .sra sizes ---------
# fqxv points. max = --max (lossless best ratio); binN = lossy quality binning.
# .sra Illumina quality is typically already 8-level, so bin8 is the like-for-like
# fidelity match; bin4/bin2 trade more fidelity for size.
fqxv_args() {  # point -> extra compress args
  case "$1" in
    max)  echo "--max" ;;
    bin8) echo "--quality-bin bin8" ;;
    bin4) echo "--quality-bin bin4" ;;
    bin2) echo "--quality-bin bin2" ;;
    *) echo "unknown point $1" >&2; return 1 ;;
  esac
}

do_run() {
  [[ -x "$FQXV_BIN" ]] || { echo "fqxv binary missing: $FQXV_BIN (cargo build --release)" >&2; exit 1; }
  # .sra sizes: prefer the cached table; fall back to a live lookup if the node
  # has network and no cache exists yet.
  declare -A SRA_B SPOTS LAYOUT
  if [[ -f "$SIZES" ]]; then
    while IFS=$'\t' read -r acc plat regime layout spots sra_bytes; do
      [[ "$acc" == accession ]] && continue
      SRA_B["$acc"]="$sra_bytes"; SPOTS["$acc"]="$spots"; LAYOUT["$acc"]="$layout"
    done < "$SIZES"
  elif command -v sracha >/dev/null 2>&1; then
    echo "[info] no $SIZES cache; looking up .sra sizes live"; do_sizes
    while IFS=$'\t' read -r acc plat regime layout spots sra_bytes; do
      [[ "$acc" == accession ]] && continue
      SRA_B["$acc"]="$sra_bytes"; SPOTS["$acc"]="$spots"; LAYOUT["$acc"]="$layout"
    done < "$SIZES"
  else
    echo "[warn] no $SIZES and no sracha: .sra column will be blank (run the sizes phase first)"
  fi

  echo -e "accession\tplatform\tregime\tlayout\tspots\tsra_bytes\tfastq_bytes\tpoint\tfqxv_bytes\tnames_bytes\tseq_bytes\tqual_bytes\tfqxv_over_sra\tfqxv_over_fastq\tc_secs\tc_rss_kb" > "$OUT"

  panel_rows | while IFS=$'\t' read -r acc plat regime cls; do
    # Resolve the FASTQ fqxv ingests. .sra bundles BOTH mates + names + spot
    # metadata, so the fair fqxv input is the whole run: pass both mates (paired,
    # per-spot interleaved — the real archive) or the single file (long read).
    local r1="$DATA_DIR/${acc}_1.fastq" r2="$DATA_DIR/${acc}_2.fastq"
    [[ -f "$r1" ]] || r1="$DATA_DIR/${acc}.fastq"
    if [[ ! -f "$r1" ]]; then echo "  [skip] $acc: $DATA_DIR/${acc}[_1].fastq missing (run fetch.sh)"; continue; fi
    local -a inputs; local fastq_bytes
    if [[ -f "$r2" ]]; then
      inputs=("$r1" "$r2"); fastq_bytes=$(( $(stat -c %s "$r1") + $(stat -c %s "$r2") ))
    else
      inputs=("$r1"); fastq_bytes="$(stat -c %s "$r1")"
    fi
    local sra_b="${SRA_B[$acc]:-}" spots="${SPOTS[$acc]:-}" layout="${LAYOUT[$acc]:--}"
    echo "==> $acc ($plat, $regime)  fastq=$(numfmt --to=iec "$fastq_bytes")  .sra=${sra_b:+$(numfmt --to=iec "$sra_b")}"

    for pt in $POINTS; do
      local comp="$WORK/${acc}.${pt}.fqxv"; rm -f "$comp"
      # shellcheck disable=SC2046  # deliberate word-split of fqxv_args
      measure "$FQXV_BIN" compress "${inputs[@]}" -o "$comp" --force --threads "$THREADS" $(fqxv_args "$pt")
      if [[ "$MEAS_RC" -ne 0 || ! -f "$comp" ]]; then
        echo "  [fail] $acc/$pt: compress exited $MEAS_RC"; rm -f "$comp"; continue
      fi
      local fqxv_bytes; fqxv_bytes="$(stat -c %s "$comp")"
      # Per-stream sizes from `fqxv info --tsv`: header then one data row of
      #   file_size reads blocks group_size seq_order quality_binning reordered \
      #   names_bytes seq_bytes qual_bytes  (indices 7/8/9).
      local names_b=-1 seq_b=-1 qual_b=-1; mapfile -t _info < <("$FQXV_BIN" info "$comp" --tsv 2>/dev/null || true)
      if [[ "${#_info[@]}" -ge 2 ]]; then
        IFS=$'\t' read -r -a _d <<<"${_info[1]}"; names_b="${_d[7]:--1}"; seq_b="${_d[8]:--1}"; qual_b="${_d[9]:--1}"
      fi
      local over_sra over_fq
      over_sra="$(awk -v f="$fqxv_bytes" -v s="${sra_b:-0}" 'BEGIN{printf (s>0)?"%.4f":"NA", (s>0)?f/s:0}')"
      over_fq="$(awk -v f="$fqxv_bytes" -v q="$fastq_bytes" 'BEGIN{printf (q>0)?"%.4f":"NA", (q>0)?f/q:0}')"
      printf '  %-6s fqxv=%-9s /.sra=%-7s /.fastq=%-7s c=%ss\n' \
        "$pt" "$(numfmt --to=iec "$fqxv_bytes")" "$over_sra" "$over_fq" "$MEAS_SECS"
      echo -e "${acc}\t${plat}\t${regime}\t${layout}\t${spots}\t${sra_b:--1}\t${fastq_bytes}\t${pt}\t${fqxv_bytes}\t${names_b}\t${seq_b}\t${qual_b}\t${over_sra}\t${over_fq}\t${MEAS_SECS}\t${MEAS_RSS_KB}" >> "$OUT"
      rm -f "$comp"
    done
  done
  echo "==> wrote $OUT"
}

case "$MODE" in
  sizes) do_sizes ;;
  run)   do_run ;;
  all)   [[ -f "$SIZES" ]] || do_sizes; do_run ;;
esac
