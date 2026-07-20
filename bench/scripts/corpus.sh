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
# ERROR_FETCH.  Modes: default (-l5 lossless), max (-l9 --order any), the Illumina
# lossy bins (bin8/bin4/bin2), and the long-read lossy bins (ont, hifi).
#
# By default the accession list spans every platform fqxv compresses — Illumina,
# Oxford Nanopore, PacBio/HiFi, and MGI/BGI — so the long-read overlap seq codec
# and the ONT/HiFi paths are exercised, not just Illumina short reads.
#
# Subcommands:
#   sample [-n N] [-s SEED] [-p PLATFORM]   write a fresh accession list (ENA);
#                                           no -p = fan out across CORPUS_PLATFORMS
#   fetch                                   sracha get every accession (login-node OK: IO bound)
#   run                                     compress/round-trip all fetched (COMPUTE NODE — srun/sbatch)
#   sbatch                                  submit `run` as a slurm array (one accession per task)
#   summary                                 print pass/fail tally + failing accessions
#
# Env: FQXV_BIN, FQXV_CORPUS_DIR (default $SCRATCH/fqxv/corpus), FQXV_MODES
#      (default "default max"), FQXV_THREADS (default 16), CORPUS_N, CORPUS_SEED,
#      CORPUS_PLATFORMS (default "ILLUMINA OXFORD_NANOPORE PACBIO_SMRT BGISEQ"),
#      CORPUS_HIFI_MODELS (PacBio instrument filter; default "Revio,Sequel II,
#      Sequel IIe" for HiFi — set empty to include CLR).
#
# Typical flow:
#   pixi run bash corpus.sh sample                 # -> $CORPUS_DIR/accessions.txt
#   pixi run bash corpus.sh fetch                  # login node
#   sbatch --export=ALL corpus.sh_sbatch ...       # or: pixi run bash corpus.sh sbatch
#   pixi run bash corpus.sh summary
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Repo root manifest: pixi is consolidated to one manifest at the repo root, with
# the harness deps under the `bench` environment (see pixi.toml).
ROOT_MANIFEST="$(cd "$HERE/../.." && pwd)/pixi.toml"
CORPUS_DIR="${FQXV_CORPUS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/corpus}"
DATA_DIR="$CORPUS_DIR/data"
WORK_DIR="$CORPUS_DIR/work"
LOG_DIR="$CORPUS_DIR/logs"
ACC_LIST="$CORPUS_DIR/accessions.txt"
META_FILE="$CORPUS_DIR/metadata.tsv"
RESULTS_TSV="$CORPUS_DIR/results.tsv"
RESULTS_LOCK="$CORPUS_DIR/results.lock"

FQXV_BIN="${FQXV_BIN:-${CARGO_TARGET_DIR:-$(cd "$HERE/../.." && pwd)/target}/release/fqxv}"
THREADS="${FQXV_THREADS:-16}"
# Second thread count for the once-per-accession determinism check. It only has to
# DIFFER from THREADS (same count wouldn't test thread-count invariance); we use
# THREADS/2 so the extra compress stays near the multi-thread speed instead of the
# old per-mode `--threads 1` re-compress that dominated wall-clock (30-60 min/mode
# on 1 Gbase ONT). Exhaustive 1-vs-N determinism lives in the crate proptests.
DET_THREADS="${FQXV_DET_THREADS:-$(( THREADS/2 > 1 ? THREADS/2 : 1 ))}"
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
FQDIGEST_SRC="$HERE/../tools/fqdigest.rs"
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
    ont)     echo "--quality-bin ont" ;;    # Nanopore 4-level lossy
    hifi)    echo "--quality-bin hifi" ;;    # PacBio HiFi 5-level lossy
    l9)      echo "-l 9" ;;
    reorder) echo "--order any" ;;
    *) echo "" ;;
  esac
}

# True for the lossy quality-binning modes (their round-trip is checked against a
# same-binned reference digest, not the lossless one). Mirrors compress_args.
is_binned_mode() {  # mode -> exit 0 if lossy-binned
  case "$1" in bin8|bin4|bin2|ont|hifi) return 0 ;; *) return 1 ;; esac
}

# Per-platform base-count window for sampling. Long-read platforms (Nanopore,
# PacBio) get a wider ceiling so a run actually carries genuine multi-kb reads
# and enough of them to exercise the long-read overlap seq codec; short-read
# platforms keep the modest default so files stay scratch-friendly.
platform_bases() {  # ENA platform -> echoes "--min-bases N --max-bases N"
  case "$1" in
    OXFORD_NANOPORE|PACBIO_SMRT) echo "--min-bases 50000000 --max-bases 2000000000" ;;
    *)                           echo "--min-bases 20000000 --max-bases 600000000" ;;
  esac
}

# Instrument-model filter per platform. ENA has no direct CCS/HiFi flag, so we
# narrow PacBio to HiFi-capable instruments (Revio is HiFi-only; Sequel II/IIe
# emit CCS) — the closest proxy for HiFi rather than legacy CLR subread runs.
# Echoes a comma-separated model list (may contain spaces) or nothing. Set
# CORPUS_HIFI_MODELS="" to sample any PacBio (CLR included).
platform_models() {  # ENA platform -> echoes comma-separated instrument models (or empty)
  case "$1" in
    PACBIO_SMRT) echo "${CORPUS_HIFI_MODELS-Revio,Sequel II,Sequel IIe}" ;;
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
  # Determinism is a property of the parallel pipeline and is essentially mode-
  # independent, so we check it ONCE per accession (on the first mode) rather than
  # re-compressing every mode single-threaded. We keep that mode's THREADS archive
  # and diff it against a DET_THREADS re-compress after the loop.
  local det_mode; det_mode="${MODES%% *}"   # first mode in the list
  local keep=""                             # kept THREADS archive of det_mode
  local mode
  for mode in $MODES; do
    local comp="$wd/$mode.fqxv" rt="$wd/$mode.rt.fastq" args t0 rc secs bytes
    args=$(compress_args "$mode")
    echo "=== compress ($mode) fqxv $args ===" >> "$log"
    t0=$SECONDS
    # shellcheck disable=SC2086
    if ! "$FQXV_BIN" compress "${INPUTS[@]}" -o "$comp" --force --threads "$THREADS" $args >> "$log" 2>&1; then
      rc=$?; secs=$((SECONDS-t0))
      echo "  $acc/$mode: FAIL_COMPRESS (rc=$rc)"
      record "$acc" "$mode" "FAIL_COMPRESS" "rc=$rc" "" "$secs"
      continue
    fi
    secs=$((SECONDS-t0))
    bytes=$(stat -c%s "$comp" 2>/dev/null || echo "")

    # Decompress (interleaved to stdout captured to file).
    echo "=== decompress ($mode) ===" >> "$log"
    if ! "$FQXV_BIN" decompress "$comp" -o "$rt" --force --threads "$THREADS" >> "$log" 2>&1; then
      rc=$?
      echo "  $acc/$mode: FAIL_DECOMPRESS (rc=$rc)"
      record "$acc" "$mode" "FAIL_DECOMPRESS" "rc=$rc" "$bytes" "$secs"
      rm -f "$comp" "$rt"; continue
    fi

    # Content round-trip.
    local got want
    got=$(record_digest "$rt")
    if is_binned_mode "$mode"; then
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

    echo "  $acc/$mode: PASS (${bytes}B, ${secs}s)"
    record "$acc" "$mode" "PASS" "-" "$bytes" "$secs"

    # Keep the det_mode archive for the single determinism check below.
    if [[ "$mode" == "$det_mode" ]]; then mv "$comp" "$wd/keep.fqxv"; keep="$wd/keep.fqxv"; else rm -f "$comp"; fi
  done

  # One determinism check per accession: re-compress det_mode at a DIFFERENT
  # (still parallel) thread count and require byte-identical output.
  if [[ -n "$keep" ]]; then
    echo "=== determinism ($det_mode, --threads $DET_THREADS vs $THREADS) ===" >> "$log"
    local compd="$wd/det.fqxv" dargs; dargs=$(compress_args "$det_mode")
    # shellcheck disable=SC2086
    if ! "$FQXV_BIN" compress "${INPUTS[@]}" -o "$compd" --force --threads "$DET_THREADS" $dargs >> "$log" 2>&1; then
      rc=$?
      echo "  $acc/det: FAIL_COMPRESS (threads=$DET_THREADS, rc=$rc)"
      record "$acc" "det:$det_mode" "FAIL_COMPRESS" "threads=$DET_THREADS rc=$rc" "" ""
    elif ! cmp -s "$keep" "$compd"; then
      echo "  $acc/det: FAIL_DET ($DET_THREADS vs $THREADS differ)" | tee -a "$log"
      record "$acc" "det:$det_mode" "FAIL_DET" "t$DET_THREADS != t$THREADS" "" ""
    else
      echo "  $acc/det: PASS ($DET_THREADS==$THREADS)"
      record "$acc" "det:$det_mode" "PASS" "-" "" ""
    fi
    rm -f "$keep" "$compd"
  fi
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
    # Default corpus spans every platform fqxv can compress so ONT and PacBio/HiFi
    # long reads are exercised alongside Illumina/MGI short reads, not just the
    # historical Illumina-only set. `-p PLATFORM` narrows it to one platform.
    SAMPLE_EXTRA=(); SAMPLE_P=""
    while [[ $# -gt 0 ]]; do case "$1" in
      -n) CORPUS_N="$2"; shift 2 ;; -s|--seed) CORPUS_SEED="$2"; shift 2 ;;
      -p|--platform) SAMPLE_P="$2"; shift 2 ;;
      *) SAMPLE_EXTRA+=("$1"); shift ;; esac; done

    # If the caller pinned a platform, honor it as a single-platform corpus;
    # otherwise fan out across CORPUS_PLATFORMS (Illumina, Nanopore, PacBio, MGI).
    if [[ -n "$SAMPLE_P" ]]; then
      PLATS=("$SAMPLE_P")
    else
      # shellcheck disable=SC2206
      PLATS=(${CORPUS_PLATFORMS:-ILLUMINA OXFORD_NANOPORE PACBIO_SMRT BGISEQ})
    fi

    : > "$ACC_LIST"; : > "$META_FILE"
    nplat=${#PLATS[@]}; idx=0
    for plat in "${PLATS[@]}"; do
      # Split N as evenly as possible across platforms (remainder to the first few).
      per=$(( CORPUS_N / nplat )); (( idx < CORPUS_N % nplat )) && per=$((per+1))
      idx=$((idx+1))
      [[ "$per" -lt 1 ]] && continue
      # Deterministic, platform-distinct seed when a base seed is given; else let
      # sample_accessions pick its own $RANDOM per platform.
      ARGS=(-n "$per" -p "$plat")
      [[ -n "$CORPUS_SEED" ]] && ARGS+=(-s "${CORPUS_SEED}-${plat}")
      # shellcheck disable=SC2206
      ARGS+=($(platform_bases "$plat"))
      pmodels=$(platform_models "$plat")
      [[ -n "$pmodels" ]] && ARGS+=(--instrument-model "$pmodels")
      [[ ${#SAMPLE_EXTRA[@]} -gt 0 ]] && ARGS+=("${SAMPLE_EXTRA[@]}")
      echo "# === $plat (n=$per) ===" | tee -a "$ACC_LIST" >> "$META_FILE"
      bash "$HERE/sample_accessions.sh" "${ARGS[@]}" >> "$ACC_LIST" 2>> "$META_FILE"
    done
    echo "wrote $(grep -cv '^\(#\|$\)' "$ACC_LIST") accessions across ${nplat} platform(s) to $ACC_LIST"
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
    # No %concurrency cap: let Slurm schedule every array task as resources free
    # up (each task is one accession, independent, and writes results under flock).
    # Set CORPUS_CONCURRENCY to reimpose a throttle (e.g. to be a good neighbor on
    # a busy partition); empty/unset = unthrottled.
    conc="${CORPUS_CONCURRENCY:-}"
    array_spec="1-${total}${conc:+%${conc}}"
    echo "submitting array $array_spec"
    sbatch --job-name=fqxv-corpus --comment=fqxv-corpus --partition="${SLURM_PARTITION:-amilan}" \
      --qos="${SLURM_QOS:-normal}" --array="$array_spec" \
      --cpus-per-task="$THREADS" --mem="${CORPUS_MEM:-32G}" --time="${CORPUS_TIME:-02:00:00}" \
      --output="$LOG_DIR/slurm-%A_%a.out" \
      --wrap="cd '$HERE' && FQXV_CORPUS_DIR='$CORPUS_DIR' FQXV_BIN='$FQXV_BIN' FQDIGEST='$FQDIGEST' FQXV_MODES='$MODES' FQXV_THREADS='$THREADS' pixi run -e bench --manifest-path '$ROOT_MANIFEST' bash corpus.sh run"
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
