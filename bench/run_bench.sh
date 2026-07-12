#!/usr/bin/env bash
# Core benchmark runner: for each dataset x tool, record compressed size,
# compress/decompress wall-time, peak RSS, and a round-trip sanity check.
#
# Emits TSV to $RESULTS_DIR/results.tsv (+ per-dataset meta.tsv). Meant to run
# INSIDE an srun/sbatch allocation on one amilan node (never the login node) so
# the throughput numbers are clean. Invoke via `pixi run bash run_bench.sh`.
#
# Env knobs:
#   FQXV_DATA_DIR     input FASTQ dir            (default $SCRATCH/fqxv/data)
#   FQXV_RESULTS_DIR  output dir                 (default $SCRATCH/fqxv/results)
#   FQXV_THREADS      threads for MT tools       (default: nproc)
#   FQXV_INPUT        r1 | cat  (R1 only, or R1+R2 concatenated)  (default r1)
#   FQXV_TOOLS        space-separated subset     (default: all)
set -euo pipefail

DATA_DIR="${FQXV_DATA_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/data}"
RESULTS_DIR="${FQXV_RESULTS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/results}"
THREADS="${FQXV_THREADS:-$(nproc)}"
INPUT_MODE="${FQXV_INPUT:-r1}"
ALL_TOOLS="fqxv gzip zstd19 xz9 fqz_comp fqzcomp5 spring"
TOOLS="${FQXV_TOOLS:-$ALL_TOOLS}"
# The fqxv binary (built with `cargo build --release`).
FQXV_BIN="${FQXV_BIN:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/target/release/fqxv}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="$RESULTS_DIR/work"

# From-source baselines (built by build_tools.sh) live here; fqzcomp5 needs its
# htscodecs shared lib on LD_LIBRARY_PATH. PgRC is built but intentionally NOT
# in ALL_TOOLS: it is a sequence-only read compressor (drops names, simplifies
# quality), so it is not comparable to full-FASTQ archivers — see README.
TOOLS_DIR="${FQXV_TOOLS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/tools}"
export PATH="$TOOLS_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$TOOLS_DIR/lib:${LD_LIBRARY_PATH:-}"

mkdir -p "$RESULTS_DIR" "$WORK"
RESULTS="$RESULTS_DIR/results.tsv"
META="$RESULTS_DIR/meta.tsv"

# GNU time for wall + peak RSS; fall back to bash timing (RSS unknown = -1).
GNU_TIME=""
for c in /usr/bin/time "$(command -v time || true)"; do
  if [[ -x "$c" ]] && "$c" -f '%e %M' true >/dev/null 2>&1; then GNU_TIME="$c"; break; fi
done

# measure CMD... -> sets MEAS_SECS, MEAS_RSS_KB
measure() {
  local tf; tf="$(mktemp)"
  if [[ -n "$GNU_TIME" ]]; then
    "$GNU_TIME" -o "$tf" -f '%e %M' "$@"
    read -r MEAS_SECS MEAS_RSS_KB < "$tf"
  else
    local t0 t1; t0="$EPOCHREALTIME"; "$@"; t1="$EPOCHREALTIME"
    MEAS_SECS="$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}')"
    MEAS_RSS_KB="-1"
  fi
  rm -f "$tf"
}

fastq_records() { echo $(( $(wc -l < "$1") / 4 )); }

# --- per-tool compress/decompress. Each sets COMP (compressed path) then RT. ---
compress() {  # tool input out_prefix
  local tool="$1" in="$2" pfx="$3"
  case "$tool" in
    fqxv)     COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --threads "$THREADS" ;;
    gzip)     COMP="$pfx.gz";  measure bash -c "pigz -p $THREADS -6 -c '$in' > '$COMP'" ;;
    zstd19)   COMP="$pfx.zst"; measure bash -c "zstd -19 --long=27 -T$THREADS -q -f -o '$COMP' '$in'" ;;
    xz9)      COMP="$pfx.xz";  measure bash -c "xz -9 -T$THREADS -c '$in' > '$COMP'" ;;
    fqz_comp) COMP="$pfx.fqz"; measure bash -c "fqz_comp < '$in' > '$COMP'" ;;
    fqzcomp5) COMP="$pfx.fqz5"; measure bash -c "fqzcomp5 < '$in' > '$COMP'" ;;
    spring)   COMP="$pfx.spring"; mkdir -p "$WORK/spring_c_$$"; measure spring -c -t "$THREADS" -i "$in" -o "$COMP" -w "$WORK/spring_c_$$/" ;;
    *) echo "unknown tool $tool" >&2; return 1 ;;
  esac
}
decompress() {  # tool comp out_rt
  local tool="$1" comp="$2" rt="$3"
  case "$tool" in
    gzip)     measure bash -c "pigz -d -p $THREADS -c '$comp' > '$rt'" ;;
    zstd19)   measure bash -c "zstd -d -q -f -o '$rt' '$comp'" ;;
    xz9)      measure bash -c "xz -d -T$THREADS -c '$comp' > '$rt'" ;;
    fqz_comp) measure bash -c "fqz_comp -d < '$comp' > '$rt'" ;;
    fqzcomp5) measure bash -c "fqzcomp5 -d < '$comp' > '$rt'" ;;
    spring)   mkdir -p "$WORK/spring_d_$$"; measure spring -d -t "$THREADS" -i "$comp" -o "$rt" -w "$WORK/spring_d_$$/" ;;
  esac
}

echo -e "dataset\ttool\torig_bytes\tcomp_bytes\tratio\tc_secs\td_secs\tc_rss_kb\td_rss_kb\trt_ok" > "$RESULTS"
echo -e "dataset\torig_bytes\tn_records\tn_bases" > "$META"

mapfile -t rows < <(grep -v '^#' "$HERE/datasets.tsv" | awk 'NF')
for row in "${rows[@]}"; do
  acc="$(awk '{print $1}' <<<"$row")"
  label="$(awk '{print $2}' <<<"$row")"

  # Resolve input (R1, or R1+R2 concatenated).
  r1="$DATA_DIR/${acc}_1.fastq"
  r2="$DATA_DIR/${acc}_2.fastq"
  [[ -f "$r1" ]] || { echo "[skip] $label: $r1 missing (run fetch.sh)"; continue; }
  if [[ "$INPUT_MODE" == "cat" && -f "$r2" ]]; then
    in="$WORK/${label}.fastq"
    [[ -f "$in" ]] || cat "$r1" "$r2" > "$in"
  else
    in="$r1"
  fi

  orig_bytes="$(stat -c %s "$in")"
  nrec="$(fastq_records "$in")"
  nbases="$(awk 'NR%4==2{b+=length($0)} END{print b+0}' "$in")"
  echo -e "${label}\t${orig_bytes}\t${nrec}\t${nbases}" >> "$META"
  echo "==> $label  ($(numfmt --to=iec "$orig_bytes"), $nrec reads, $(numfmt --to=iec "$nbases") bases)"

  for tool in $TOOLS; do
    # Map tool label -> binary name (labels carry version digits / differ in case).
    case "$tool" in
      gzip) bin=pigz ;; zstd19) bin=zstd ;; xz9) bin=xz ;; pgrc) bin=PgRC ;;
      *) bin="$tool" ;;
    esac
    command -v "$bin" >/dev/null 2>&1 || { echo "  [miss] $tool ($bin)"; continue; }
    pfx="$WORK/${label}.${tool}"; rt="$WORK/${label}.${tool}.rt.fastq"
    rm -f "$pfx".* "$rt"

    compress "$tool" "$in" "$pfx"; c_secs="$MEAS_SECS"; c_rss="$MEAS_RSS_KB"
    comp_bytes="$(stat -c %s "$COMP" 2>/dev/null || echo 0)"
    decompress "$tool" "$COMP" "$rt"; d_secs="$MEAS_SECS"; d_rss="$MEAS_RSS_KB"

    rt_ok="no"
    [[ -f "$rt" && "$(fastq_records "$rt")" == "$nrec" ]] && rt_ok="yes"
    ratio="$(awk -v o="$orig_bytes" -v c="$comp_bytes" 'BEGIN{printf "%.3f", (c>0)?o/c:0}')"
    printf '  %-9s ratio=%-6s c=%ss d=%ss rss=%sK rt=%s\n' "$tool" "$ratio" "$c_secs" "$d_secs" "$c_rss" "$rt_ok"
    echo -e "${label}\t${tool}\t${orig_bytes}\t${comp_bytes}\t${ratio}\t${c_secs}\t${d_secs}\t${c_rss}\t${d_rss}\t${rt_ok}" >> "$RESULTS"
    rm -f "$pfx".* "$rt"
  done
done

echo "==> wrote $RESULTS"
