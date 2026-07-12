#!/usr/bin/env bash
# Core benchmark runner: for each dataset x tool, record compressed size,
# compress/decompress wall-time, peak RSS, a *content* round-trip check, and
# (for fqxv) per-stream byte sizes plus a thread-determinism check.
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
# fqxv, fqxv9 (level 9), fqxv-reorder (--reorder --keep-order), fqxv-bin4 (lossy
# 4-bin quality) all share one binary; the rest are external baselines.
ALL_TOOLS="fqxv fqxv9 fqxv-reorder fqxv-bin4 gzip zstd19 xz9 fqz_comp fqzcomp5 spring"
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

# Order-independent record-multiset digest: emit name<TAB>seq[<TAB>qual] per
# record, sort, and hash. This verifies *content* losslessness (names, bases,
# and — unless mode=noqual — qualities), and because it sorts, it is invariant
# to read reordering (SPRING, `fqxv --reorder`). The `+` line (record line 3) is
# deliberately excluded: fqxv normalizes it, which is the one documented, lossy-
# by-design deviation. mode=noqual drops quality for lossy-quality tools.
record_digest() {  # file mode(full|noqual)
  awk -v mode="$2" '
    NR%4==1{n=$0} NR%4==2{s=$0}
    NR%4==0{ if(mode=="noqual") print n"\t"s; else print n"\t"s"\t"$0 }
  ' "$1" | LC_ALL=C sort | md5sum | cut -d' ' -f1
}

is_fqxv() { [[ "$1" == fqxv || "$1" == fqxv-* || "$1" == fqxv[0-9] ]]; }
# Lossy tools: content check ignores quality (mode=noqual).
is_lossy() { [[ "$1" == fqxv-bin* ]]; }

# --- per-tool compress/decompress. Each sets COMP (compressed path) then RT. ---
compress() {  # tool input out_prefix
  local tool="$1" in="$2" pfx="$3"
  case "$tool" in
    fqxv)          COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --threads "$THREADS" ;;
    fqxv9)         COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" -l 9 --threads "$THREADS" ;;
    fqxv-reorder)  COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --reorder --keep-order --threads "$THREADS" ;;
    fqxv-bin4)     COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --quality-bin bin4 --threads "$THREADS" ;;
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
    fqxv|fqxv9|fqxv-reorder|fqxv-bin4) measure "$FQXV_BIN" decompress "$comp" -o "$rt" --threads "$THREADS" ;;
    gzip)     measure bash -c "pigz -d -p $THREADS -c '$comp' > '$rt'" ;;
    zstd19)   measure bash -c "zstd -d -q -f -o '$rt' '$comp'" ;;
    xz9)      measure bash -c "xz -d -T$THREADS -c '$comp' > '$rt'" ;;
    fqz_comp) measure bash -c "fqz_comp -d < '$comp' > '$rt'" ;;
    fqzcomp5) measure bash -c "fqzcomp5 -d < '$comp' > '$rt'" ;;
    spring)   mkdir -p "$WORK/spring_d_$$"; measure spring -d -t "$THREADS" -i "$comp" -o "$rt" -w "$WORK/spring_d_$$/" ;;
  esac
}

# results.tsv columns: per-stream sizes are fqxv-only (-1 for other tools);
# rt_ok is now a *content* multiset check; deterministic is a 1-thread vs
# N-thread byte-identity check (fqxv only, else n/a).
echo -e "dataset\ttool\torig_bytes\tcomp_bytes\tratio\tc_secs\td_secs\tc_rss_kb\td_rss_kb\tnames_bytes\tseq_bytes\tqual_bytes\trt_ok\tdeterministic" > "$RESULTS"
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

  # Reference content digests for this input (computed once). full = with
  # quality (lossless tools); noqual = names+bases only (lossy-quality tools).
  in_full="$(record_digest "$in" full)"
  in_noqual="$(record_digest "$in" noqual)"

  for tool in $TOOLS; do
    # Map tool label -> binary name to probe availability.
    case "$tool" in
      gzip) bin=pigz ;; zstd19) bin=zstd ;; xz9) bin=xz ;; pgrc) bin=PgRC ;;
      fqxv|fqxv9|fqxv-reorder|fqxv-bin4) bin="$FQXV_BIN" ;;
      *) bin="$tool" ;;
    esac
    if is_fqxv "$tool"; then
      [[ -x "$bin" ]] || { echo "  [miss] $tool ($bin — run: cargo build --release)"; continue; }
    else
      command -v "$bin" >/dev/null 2>&1 || { echo "  [miss] $tool ($bin)"; continue; }
    fi
    pfx="$WORK/${label}.${tool}"; rt="$WORK/${label}.${tool}.rt.fastq"
    rm -f "$pfx".* "$rt"

    compress "$tool" "$in" "$pfx"; c_secs="$MEAS_SECS"; c_rss="$MEAS_RSS_KB"
    comp_bytes="$(stat -c %s "$COMP" 2>/dev/null || echo 0)"
    decompress "$tool" "$COMP" "$rt"; d_secs="$MEAS_SECS"; d_rss="$MEAS_RSS_KB"

    # Per-stream sizes (fqxv only; others get -1).
    names_b=-1; seq_b=-1; qual_b=-1
    if is_fqxv "$tool"; then
      # info --tsv: header line then one data line of stable columns:
      #   file_size reads blocks group_size seq_order quality_binning \
      #   reordered names_bytes seq_bytes qual_bytes
      mapfile -t _info < <("$FQXV_BIN" info "$COMP" --tsv 2>/dev/null || true)
      if [[ "${#_info[@]}" -ge 2 ]]; then
        IFS=$'\t' read -r -a _d <<<"${_info[1]}"
        names_b="${_d[7]:--1}"; seq_b="${_d[8]:--1}"; qual_b="${_d[9]:--1}"
      fi
    fi

    # Content round-trip: order-independent multiset digest. Lossy-quality tools
    # are checked without quality (names + bases must still be exact).
    rt_ok="no"
    if [[ -f "$rt" ]]; then
      if is_lossy "$tool"; then
        [[ "$(record_digest "$rt" noqual)" == "$in_noqual" ]] && rt_ok="yes"
      else
        [[ "$(record_digest "$rt" full)" == "$in_full" ]] && rt_ok="yes"
      fi
    fi

    # Determinism: fqxv must be byte-identical regardless of thread count
    # (core invariant). Compare a 1-thread build of the same archive.
    deterministic="n/a"
    if is_fqxv "$tool"; then
      det1="$WORK/${label}.${tool}.det1.fqxv"
      case "$tool" in
        fqxv)         "$FQXV_BIN" compress "$in" -o "$det1" --threads 1 >/dev/null 2>&1 || true ;;
        fqxv9)        "$FQXV_BIN" compress "$in" -o "$det1" -l 9 --threads 1 >/dev/null 2>&1 || true ;;
        fqxv-reorder) "$FQXV_BIN" compress "$in" -o "$det1" --reorder --keep-order --threads 1 >/dev/null 2>&1 || true ;;
        fqxv-bin4)    "$FQXV_BIN" compress "$in" -o "$det1" --quality-bin bin4 --threads 1 >/dev/null 2>&1 || true ;;
      esac
      if [[ -f "$det1" ]] && cmp -s "$det1" "$COMP"; then deterministic="yes"; else deterministic="no"; fi
      rm -f "$det1"
    fi

    ratio="$(awk -v o="$orig_bytes" -v c="$comp_bytes" 'BEGIN{printf "%.3f", (c>0)?o/c:0}')"
    printf '  %-13s ratio=%-6s c=%ss d=%ss rss=%sK rt=%s det=%s\n' \
      "$tool" "$ratio" "$c_secs" "$d_secs" "$c_rss" "$rt_ok" "$deterministic"
    echo -e "${label}\t${tool}\t${orig_bytes}\t${comp_bytes}\t${ratio}\t${c_secs}\t${d_secs}\t${c_rss}\t${d_rss}\t${names_b}\t${seq_b}\t${qual_b}\t${rt_ok}\t${deterministic}" >> "$RESULTS"
    rm -f "$pfx".* "$rt"
  done

  # --- fqxv paired self-check: exercise per-spot interleaving of R1+R2 (the
  # container feature `cat` mode bypasses). Not part of the comparison table;
  # recorded as a `fqxv-paired` row so its size/streams/losslessness are tracked.
  if [[ -f "$r2" ]] && is_fqxv fqxv && [[ -x "$FQXV_BIN" ]] && [[ " $TOOLS " == *" fqxv "* ]]; then
    pcat="$WORK/${label}.paircat.fastq"; cat "$r1" "$r2" > "$pcat"
    p_orig="$(stat -c %s "$pcat")"
    p_full="$(record_digest "$pcat" full)"
    pfx="$WORK/${label}.fqxv-paired"; COMP="$pfx.fqxv"; rt="$WORK/${label}.fqxv-paired.rt"
    rm -f "$pfx".* "$rt"_*
    measure "$FQXV_BIN" compress "$r1" "$r2" -o "$COMP" --threads "$THREADS"; c_secs="$MEAS_SECS"; c_rss="$MEAS_RSS_KB"
    comp_bytes="$(stat -c %s "$COMP" 2>/dev/null || echo 0)"
    # Restore both mates and concatenate to compare the multiset against R1+R2.
    measure "$FQXV_BIN" decompress "$COMP" --split "$rt" --threads "$THREADS"; d_secs="$MEAS_SECS"; d_rss="$MEAS_RSS_KB"
    rt_ok="no"
    if [[ -f "${rt}_1.fastq" && -f "${rt}_2.fastq" ]]; then
      cat "${rt}_1.fastq" "${rt}_2.fastq" > "${rt}.all"
      [[ "$(record_digest "${rt}.all" full)" == "$p_full" ]] && rt_ok="yes"
    fi
    mapfile -t _info < <("$FQXV_BIN" info "$COMP" --tsv 2>/dev/null || true)
    names_b=-1; seq_b=-1; qual_b=-1
    if [[ "${#_info[@]}" -ge 2 ]]; then
      IFS=$'\t' read -r -a _d <<<"${_info[1]}"; names_b="${_d[7]:--1}"; seq_b="${_d[8]:--1}"; qual_b="${_d[9]:--1}"
    fi
    ratio="$(awk -v o="$p_orig" -v c="$comp_bytes" 'BEGIN{printf "%.3f", (c>0)?o/c:0}')"
    printf '  %-13s ratio=%-6s c=%ss d=%ss rt=%s (R1+R2 interleaved)\n' "fqxv-paired" "$ratio" "$c_secs" "$d_secs" "$rt_ok"
    echo -e "${label}-paired\tfqxv-paired\t${p_orig}\t${comp_bytes}\t${ratio}\t${c_secs}\t${d_secs}\t${c_rss}\t${d_rss}\t${names_b}\t${seq_b}\t${qual_b}\t${rt_ok}\tn/a" >> "$RESULTS"
    rm -f "$pcat" "$pfx".* "${rt}"_* "${rt}.all"
  fi
done

echo "==> wrote $RESULTS"
