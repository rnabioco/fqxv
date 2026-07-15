#!/usr/bin/env bash
#
# fqxv robustness corpus — fetch a spread of real SRA runs and hammer fqxv on
# each, hunting for crashes, round-trip corruption, and thread-nondeterminism.
#
# This is NOT the performance benchmark (run_bench.sh, a curated 4-dataset
# comparison against the field). This is a correctness net: throw a random pile
# of real-world FASTQ at fqxv and confirm every archive round-trips byte-for-
# content and builds identically regardless of thread count.
#
# For each accession we fetch FASTQ via sracha, then per compression MODE run:
#   compress -> decompress -> order-independent content round-trip
#   compress --threads 1 -> byte-compare with the many-threaded archive (det)
# and classify: PASS / FAIL_RT / FAIL_DET / FAIL_COMPRESS / FAIL_DECOMPRESS /
# ERROR_FETCH.  Modes: default (-l5 lossless), max (-l9 --order any), and
# optionally bin8 (lossy Illumina 8-level).
#
# Subcommands:
#   sample [-n N] [-s SEED] [-p PLATFORM]   write a fresh accession list (ENA)
#   fetch                                   sracha get every accession (login-node OK: IO bound)
#   run                                     compress/round-trip all fetched (COMPUTE NODE — srun/sbatch)
#   sbatch                                  submit `run` as a slurm array (one accession per task)
#   summary                                 print pass/fail tally + failing accessions
#
# Env: FQXV_BIN, FQXV_CORPUS_DIR (default $SCRATCH/fqxv/corpus), FQXV_MODES
#      (default "default max"), FQXV_THREADS (default 16), CORPUS_N, CORPUS_SEED.
#
# Typical flow:
#   pixi run bash corpus.sh sample                 # -> $CORPUS_DIR/accessions.txt
#   pixi run bash corpus.sh fetch                  # login node
#   sbatch --export=ALL corpus.sh_sbatch ...       # or: pixi run bash corpus.sh sbatch
#   pixi run bash corpus.sh summary
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CORPUS_DIR="${FQXV_CORPUS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/corpus}"
DATA_DIR="$CORPUS_DIR/data"
WORK_DIR="$CORPUS_DIR/work"
LOG_DIR="$CORPUS_DIR/logs"
ACC_LIST="$CORPUS_DIR/accessions.txt"
META_FILE="$CORPUS_DIR/metadata.tsv"
RESULTS_TSV="$CORPUS_DIR/results.tsv"
RESULTS_LOCK="$CORPUS_DIR/results.lock"

FQXV_BIN="${FQXV_BIN:-${CARGO_TARGET_DIR:-$(cd "$HERE/.." && pwd)/target}/release/fqxv}"
THREADS="${FQXV_THREADS:-16}"
MODES="${FQXV_MODES:-default max}"
FETCH_THREADS="${FQXV_FETCH_THREADS:-8}"
CORPUS_N="${CORPUS_N:-20}"
CORPUS_SEED="${CORPUS_SEED:-}"

mkdir -p "$CORPUS_DIR" "$DATA_DIR" "$WORK_DIR" "$LOG_DIR"

# ---------- content round-trip helpers ----------
# Order-independent record-multiset digest of (name, seq, qual), excluding the
# `+` line (fqxv normalizes it — the one documented lossy-by-design deviation).
# Uses the fqdigest Rust tool (single O(n) pass, bounded memory, no sort), built
# on demand from bench/fqdigest.rs.
FQDIGEST="${FQDIGEST:-${SCRATCH:-$HOME/scratch}/fqxv/tools/bin/fqdigest}"
FQDIGEST_SRC="$HERE/fqdigest.rs"
ensure_fqdigest() {
  if [[ ! -x "$FQDIGEST" || "$FQDIGEST_SRC" -nt "$FQDIGEST" ]]; then
    mkdir -p "$(dirname "$FQDIGEST")"
    rustc -O --edition 2021 "$FQDIGEST_SRC" -o "$FQDIGEST"
  fi
}
ensure_fqdigest

record_digest() {  # file...
  "$FQDIGEST" "$@"
}
# Same, but pass each quality byte through fqxv's bin table first — the expected
# content of a correct lossy round-trip. fqdigest's --bin mirrors QualityBinning::apply.
record_digest_binned() {  # scheme file...
  local scheme="$1"; shift
  "$FQDIGEST" --bin "$scheme" "$@"
}

# Discover the FASTQ inputs sracha wrote for an accession: numbered members
# (_1/_2/_3/_4, interleaved per spot) if present, else the bare single-end file.
inputs_for() {  # accession -> echoes space-separated paths (empty if none)
  local acc="$1" f found=()
  for f in "$DATA_DIR/${acc}"_[0-9].fastq; do [[ -f "$f" ]] && found+=("$f"); done
  if [[ ${#found[@]} -eq 0 && -f "$DATA_DIR/${acc}.fastq" ]]; then found=("$DATA_DIR/${acc}.fastq"); fi
  echo "${found[@]:-}"
}

compress_args() {  # mode -> echoes fqxv compress flags
  case "$1" in
    default) echo "" ;;
    max)     echo "--max" ;;
    bin8)    echo "--quality-bin bin8" ;;
    bin4)    echo "--quality-bin bin4" ;;
    bin2)    echo "--quality-bin bin2" ;;
    l9)      echo "-l 9" ;;
    reorder) echo "--order any" ;;
    *) echo "" ;;
  esac
}

record() {  # accession mode status note comp_bytes secs
  ( flock -x 200
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$1" "$2" "$3" "$4" "$5" "$6" >> "$RESULTS_TSV"
  ) 200>"$RESULTS_LOCK"
}

# ---------- one accession: fetch (if needed) + all modes ----------
process_accession() {
  local acc="$1"
  local log="$LOG_DIR/${acc}.log"
  : > "$log"

  # Fetch if not already present (idempotent; fetch subcommand pre-populates).
  local inputs
  inputs=$(inputs_for "$acc")
  if [[ -z "$inputs" ]]; then
    echo "=== sracha get $acc ===" >> "$log"
    if ! sracha get "$acc" --output-dir "$DATA_DIR" --threads "$FETCH_THREADS" \
            --split split-3 --no-gzip --no-progress >> "$log" 2>&1; then
      echo "  $acc: ERROR_FETCH"
      record "$acc" "-" "ERROR_FETCH" "sracha get failed" "" ""
      return
    fi
    inputs=$(inputs_for "$acc")
  fi
  if [[ -z "$inputs" ]]; then
    echo "  $acc: ERROR_FETCH (no fastq produced)"
    record "$acc" "-" "ERROR_FETCH" "no fastq after get" "" ""
    return
  fi
  # shellcheck disable=SC2206
  local INPUTS=($inputs)
  local reads; reads=$(( $(wc -l < "${INPUTS[0]}") / 4 ))
  echo "  $acc: inputs=${#INPUTS[@]} reads/member~$reads" | tee -a "$log"

  # Reference digests over the concatenated inputs (order-independent).
  local ref_lossless ref_bin8=""
  ref_lossless=$(record_digest "${INPUTS[@]}")

  local wd="$WORK_DIR/$acc"; mkdir -p "$wd"
  local mode
  for mode in $MODES; do
    local comp="$wd/$mode.fqxv" rt="$wd/$mode.rt.fastq" args t0 rc secs bytes
    args=$(compress_args "$mode")
    echo "=== compress ($mode) fqxv $args ===" >> "$log"
    t0=$SECONDS
    # shellcheck disable=SC2086
    if ! "$FQXV_BIN" compress "${INPUTS[@]}" -o "$comp" --threads "$THREADS" $args >> "$log" 2>&1; then
      rc=$?; secs=$((SECONDS-t0))
      echo "  $acc/$mode: FAIL_COMPRESS (rc=$rc)"
      record "$acc" "$mode" "FAIL_COMPRESS" "rc=$rc" "" "$secs"
      continue
    fi
    secs=$((SECONDS-t0))
    bytes=$(stat -c%s "$comp" 2>/dev/null || echo "")

    # Decompress (interleaved to stdout captured to file).
    echo "=== decompress ($mode) ===" >> "$log"
    if ! "$FQXV_BIN" decompress "$comp" -o "$rt" --threads "$THREADS" >> "$log" 2>&1; then
      rc=$?
      echo "  $acc/$mode: FAIL_DECOMPRESS (rc=$rc)"
      record "$acc" "$mode" "FAIL_DECOMPRESS" "rc=$rc" "$bytes" "$secs"
      rm -f "$comp" "$rt"; continue
    fi

    # Content round-trip.
    local got want
    got=$(record_digest "$rt")
    if [[ "$mode" == bin8 || "$mode" == bin4 || "$mode" == bin2 ]]; then
      want=$(record_digest_binned "$mode" "${INPUTS[@]}")
    else
      want="$ref_lossless"
    fi
    if [[ "$got" != "$want" ]]; then
      echo "  $acc/$mode: FAIL_RT (want=$want got=$got)"
      echo "content digest mismatch: want=$want got=$got" >> "$log"
      record "$acc" "$mode" "FAIL_RT" "want=$want got=$got" "$bytes" "$secs"
      rm -f "$rt"; continue   # keep $comp for post-mortem
    fi
    rm -f "$rt"

    # Determinism: single-thread archive must be byte-identical.
    echo "=== determinism ($mode, --threads 1) ===" >> "$log"
    local comp1="$wd/$mode.t1.fqxv"
    # shellcheck disable=SC2086
    if ! "$FQXV_BIN" compress "${INPUTS[@]}" -o "$comp1" --threads 1 $args >> "$log" 2>&1; then
      rc=$?
      echo "  $acc/$mode: FAIL_COMPRESS (threads=1, rc=$rc)"
      record "$acc" "$mode" "FAIL_COMPRESS" "threads=1 rc=$rc" "$bytes" "$secs"
      rm -f "$comp"; continue
    fi
    if ! cmp -s "$comp" "$comp1"; then
      echo "  $acc/$mode: FAIL_DET (threads 1 vs $THREADS differ)"
      echo "determinism: --threads 1 archive differs from --threads $THREADS" >> "$log"
      record "$acc" "$mode" "FAIL_DET" "t1!=t$THREADS bytes $(stat -c%s "$comp1")!=$bytes" "$bytes" "$secs"
      continue   # keep both archives for post-mortem
    fi
    rm -f "$comp" "$comp1"

    echo "  $acc/$mode: PASS (${bytes}B, ${secs}s)"
    record "$acc" "$mode" "PASS" "-" "$bytes" "$secs"
  done
  rmdir "$wd" 2>/dev/null || true
}

# ---------- summary ----------
print_summary() {
  [[ -f "$RESULTS_TSV" ]] || { echo "no results yet ($RESULTS_TSV)"; return; }
  echo "=== corpus summary ($RESULTS_TSV) ==="
  awk -F'\t' 'NR>1{c[$3]++; tot++} END{ printf "  total rows: %d\n", tot; for(k in c) printf "  %-16s %d\n", k, c[k] }' "$RESULTS_TSV" | sort
  echo
  echo "non-PASS:"
  awk -F'\t' 'NR>1 && $3!="PASS"{printf "  %-14s %-8s %-16s %s\n", $1, $2, $3, $4}' "$RESULTS_TSV"
}

# ---------- dispatch ----------
cmd="${1:-}"; shift || true
case "$cmd" in
  sample)
    SAMPLE_EXTRA=()
    while [[ $# -gt 0 ]]; do case "$1" in
      -n) CORPUS_N="$2"; shift 2 ;; -s|--seed) CORPUS_SEED="$2"; shift 2 ;;
      *) SAMPLE_EXTRA+=("$1"); shift ;; esac; done
    ARGS=(-n "$CORPUS_N"); [[ -n "$CORPUS_SEED" ]] && ARGS+=(-s "$CORPUS_SEED")
    [[ ${#SAMPLE_EXTRA[@]} -gt 0 ]] && ARGS+=("${SAMPLE_EXTRA[@]}")
    bash "$HERE/sample_accessions.sh" "${ARGS[@]}" > "$ACC_LIST" 2> "$META_FILE"
    echo "wrote $(grep -cv '^$' "$ACC_LIST") accessions to $ACC_LIST"
    cat "$META_FILE" >&2
    ;;
  fetch)
    [[ -f "$ACC_LIST" ]] || { echo "no $ACC_LIST — run 'corpus.sh sample' first" >&2; exit 1; }
    mapfile -t accs < <(awk 'NF && !/^#/' "$ACC_LIST")
    echo "==> fetching ${#accs[@]} accession(s) to $DATA_DIR"
    for acc in "${accs[@]}"; do
      if [[ -n "$(inputs_for "$acc")" ]]; then echo "  [skip] $acc"; continue; fi
      echo "  [get ] $acc"
      sracha get "$acc" --output-dir "$DATA_DIR" --threads "$FETCH_THREADS" \
        --split split-3 --no-gzip --no-progress || echo "  [FAIL] $acc" >&2
    done
    ;;
  run)
    [[ -x "$FQXV_BIN" ]] || { echo "no fqxv binary at $FQXV_BIN" >&2; exit 1; }
    [[ -f "$RESULTS_TSV" ]] || printf 'accession\tmode\tstatus\tnote\tbytes\tsecs\n' > "$RESULTS_TSV"
    if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
      acc=$(awk 'NF && !/^#/' "$ACC_LIST" | sed -n "${SLURM_ARRAY_TASK_ID}p")
      [[ -z "$acc" ]] && { echo "no accession at index $SLURM_ARRAY_TASK_ID" >&2; exit 1; }
      echo "# array task $SLURM_ARRAY_TASK_ID on $(hostname) -> $acc"
      process_accession "$acc"
    else
      mapfile -t accs < <(awk 'NF && !/^#/' "$ACC_LIST")
      # Skip accessions already fully recorded (all modes present & none pending).
      i=0
      for acc in "${accs[@]}"; do
        i=$((i+1)); echo "[$i/${#accs[@]}] $acc"
        process_accession "$acc"
      done
      print_summary
    fi
    ;;
  sbatch)
    total=$(awk 'NF && !/^#/' "$ACC_LIST" | wc -l)
    [[ "$total" -ge 1 ]] || { echo "no accessions in $ACC_LIST" >&2; exit 1; }
    printf 'accession\tmode\tstatus\tnote\tbytes\tsecs\n' > "$RESULTS_TSV"
    conc="${CORPUS_CONCURRENCY:-6}"
    echo "submitting array 1-$total%$conc"
    sbatch --job-name=fqxv-corpus --comment=fqxv-corpus --partition="${SLURM_PARTITION:-amilan}" \
      --qos="${SLURM_QOS:-normal}" --array="1-${total}%${conc}" \
      --cpus-per-task="$THREADS" --mem="${CORPUS_MEM:-32G}" --time="${CORPUS_TIME:-02:00:00}" \
      --output="$LOG_DIR/slurm-%A_%a.out" \
      --wrap="cd '$HERE' && FQXV_CORPUS_DIR='$CORPUS_DIR' FQXV_BIN='$FQXV_BIN' FQDIGEST='$FQDIGEST' FQXV_MODES='$MODES' FQXV_THREADS='$THREADS' pixi run bash corpus.sh run"
    ;;
  build-digest)
    command -v rustc >/dev/null || { echo "rustc not on PATH (need the rust toolchain)" >&2; exit 1; }
    mkdir -p "$(dirname "$FQDIGEST")"
    echo "rustc -O --edition 2021 $FQDIGEST_SRC -> $FQDIGEST"
    rustc -O --edition 2021 "$FQDIGEST_SRC" -o "$FQDIGEST"
    "$FQDIGEST" --help >/dev/null 2>&1 && echo "ok: $FQDIGEST"
    ;;
  summary) print_summary ;;
  *)
    sed -n '3,40p' "$0" | sed 's/^# \{0,1\}//'
    exit 0
    ;;
esac
