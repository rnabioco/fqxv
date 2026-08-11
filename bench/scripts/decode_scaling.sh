#!/usr/bin/env bash
# Decode thread-scaling benchmark: fqxv decompress vs gzip / pigz / zstd.
#
# Motivation (BINSEQ, doi:10.1371/journal.pcbi.1014181): a sequencing format is
# an *analysis substrate*, and gzip-FASTQ consumers saturate at 2-4 threads
# because DEFLATE inflation is inherently serial, while block-based formats
# scale near-linearly with cores. fqxv's container is independent,
# parallel-decodable blocks, so its decode should scale the same way — this
# script measures that curve explicitly against gzip/pigz/zstd decoding the
# same FASTQ, and records each tool's compressed size next to its decode speed
# so the ratio-vs-throughput tradeoff reads from one table.
#
# For each dataset x tool x thread count it times a FULL decompression with the
# output forced through a pipe to `wc -c` — never to /dev/null (tools can
# fast-path a null sink; this repo has been burned by that) and never to disk
# (that would time the filesystem). The wc -c byte count doubles as a
# correctness check: every tool at every thread count must reproduce exactly
# the same number of FASTQ bytes.
#
# Compressed inputs are synthesized once (idempotent, mtime-checked) from the
# staged FASTQ in FQXV_DATA_DIR, and only for the tools actually selected:
#   .fqxv      built by THIS build of fqxv (paired datasets as one
#              spot-interleaved archive; decode streams interleaved FASTQ
#              with -Z) — the default operating point
#   .max.fqxv  the same, compressed with --max (the best-ratio preset; the
#              same block-parallel decode path, so the scaling claim is
#              measured on the archive people actually keep)
#   .gz        pigz -6 (the harness's "gzip" baseline of record; pigz emits
#              standard DEFLATE, so `gzip -dc` decodes it)
#   .zst       zstd -19 --long=27 (toolsets.sh's zstd19 point)
#
# Meant to run INSIDE a dedicated sbatch allocation
# (slurm/decode_scaling.sbatch) under the bench pixi env:
#   pixi run -e bench bash scripts/decode_scaling.sh
#
# Env knobs:
#   FQXV_DATA_DIR         staged FASTQ dir       (default $SCRATCH/fqxv/data)
#   FQXV_RESULTS_DIR      output dir             (default $SCRATCH/fqxv/results)
#   FQXV_BIN              fqxv binary            (default <repo>/target/release/fqxv)
#   FQXV_DECODE_DATASETS  accessions             (default "DRR174812 DRR205413 SRR2627175")
#   FQXV_THREAD_STEPS     thread counts          (default "1 2 4 8 16 32", capped at nproc)
#   FQXV_DECODE_TOOLS     subset of "fqxv fqxv-max gzip pigz zstd"
#   FQXV_REPS             timed repetitions      (default 3; the median is reported)
#   FQXV_REPREP=1         force re-synthesis of the compressed inputs
#
# Output: $FQXV_RESULTS_DIR/decode_scaling.tsv — `# `-prefixed metadata header
# (node, CPU, commit, tool versions), then one row per cell:
#   dataset  tool  threads  seconds  mb_per_s  bytes_out  bytes_ok  comp_bytes  reps
set -euo pipefail

DATA_DIR="${FQXV_DATA_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/data}"
RESULTS_DIR="${FQXV_RESULTS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/results}"
DATASETS="${FQXV_DECODE_DATASETS:-DRR174812 DRR205413 SRR2627175}"
STEPS="${FQXV_THREAD_STEPS:-1 2 4 8 16 32}"
TOOLS="${FQXV_DECODE_TOOLS:-fqxv fqxv-max gzip pigz zstd}"
REPS="${FQXV_REPS:-3}"
REPREP="${FQXV_REPREP:-0}"
# Same binary resolution as run_bench.sh: cargo honors CARGO_TARGET_DIR, so
# resolve the location cargo actually wrote to, not a hardcoded ROOT/target.
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FQXV_BIN="${FQXV_BIN:-${CARGO_TARGET_DIR:-$ROOT/target}/release/fqxv}"

WORK="$RESULTS_DIR/decode_scaling"
OUT_TSV="$RESULTS_DIR/decode_scaling.tsv"
mkdir -p "$WORK"

[[ -x "$FQXV_BIN" ]] || {
  echo "error: fqxv binary not found at $FQXV_BIN (cargo build --release -p fqxv-cli)" >&2
  exit 1
}

# hyperfine is a user-level cargo install on this cluster, not a pixi dep.
command -v hyperfine >/dev/null 2>&1 || export PATH="$HOME/.cargo/bin:$PATH"
HYPERFINE="$(command -v hyperfine || true)"
# The hyperfine path parses its JSON export with python3; without python3 fall
# back to the manual timing loop rather than half-working.
if [[ -n "$HYPERFINE" ]] && ! command -v python3 >/dev/null 2>&1; then HYPERFINE=""; fi

NPROC="$(nproc)"
RUN_STEPS=""
for n in $STEPS; do [[ "$n" -le "$NPROC" ]] && RUN_STEPS="$RUN_STEPS $n"; done
RUN_STEPS="${RUN_STEPS# }"
[[ -n "$RUN_STEPS" ]] || { echo "error: no thread step <= nproc ($NPROC)" >&2; exit 1; }

ver() { "$@" 2>&1 | head -n1; }  # pigz prints --version to stderr

if [[ -n "$HYPERFINE" ]]; then
  TIMER_DESC="$(ver "$HYPERFINE" --version) (median of $REPS runs, 1 warmup)"
else
  TIMER_DESC="manual wall-clock loop (median of $REPS runs, 1 warmup)"
fi

{
  echo "# fqxv decode thread-scaling — $(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "# node: $(hostname)  nproc: $NPROC  slurm_job: ${SLURM_JOB_ID:-none}"
  echo "# cpu: $(awk -F': ' '/model name/{print $2; exit}' /proc/cpuinfo)"
  echo "# commit: $(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)"
  echo "# fqxv: $(ver "$FQXV_BIN" --version)  [$FQXV_BIN]"
  echo "# gzip: $(ver gzip --version)"
  echo "# pigz: $(ver pigz --version)"
  echo "# zstd: $(ver zstd --version)"
  echo "# timer: $TIMER_DESC"
  echo "# thread_steps: $RUN_STEPS  tools: $TOOLS"
  printf 'dataset\ttool\tthreads\tseconds\tmb_per_s\tbytes_out\tbytes_ok\tcomp_bytes\treps\n'
} > "$OUT_TSV"

# fresh <out> <dep>... -> true when <out> exists, is non-empty, and is newer
# than every dependency (and re-prep was not forced).
fresh() {
  local out="$1"; shift
  [[ "$REPREP" == 1 ]] && return 1
  [[ -s "$out" ]] || return 1
  local d
  for d in "$@"; do [[ "$out" -nt "$d" ]] || return 1; done
  return 0
}

# time_cell <cmd-with-pipe-to-wc> -> CELL_SECS (median wall seconds).
time_cell() {
  local cmd="$1"
  if [[ -n "$HYPERFINE" ]]; then
    local json; json="$(mktemp)"
    "$HYPERFINE" --shell=bash --warmup 1 --runs "$REPS" --style none \
      --export-json "$json" "$cmd" >/dev/null
    CELL_SECS="$(python3 - "$json" <<'PY'
import json, sys
print(f"{json.load(open(sys.argv[1]))['results'][0]['median']:.3f}")
PY
    )"
    rm -f "$json"
  else
    local r t0 t1 times=""
    bash -c "$cmd"   # warmup
    for ((r = 0; r < REPS; r++)); do
      t0="$EPOCHREALTIME"; bash -c "$cmd"; t1="$EPOCHREALTIME"
      times="$times$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.3f", b-a}')"$'\n'
    done
    CELL_SECS="$(printf '%s' "$times" | sort -g | awk -v n="$REPS" '
      {v[NR]=$1} END{m=int((n+1)/2); if (n%2) print v[m]; else printf "%.3f", (v[m]+v[m+1])/2}')"
  fi
}

FAILED_BYTES=0

for acc in $DATASETS; do
  # Resolve staged input: paired _1/_2, else single (sracha writes single-end
  # runs as either <acc>.fastq or <acc>_0.fastq).
  r1="$DATA_DIR/${acc}_1.fastq" r2="$DATA_DIR/${acc}_2.fastq" single=""
  if [[ -f "$r1" && -f "$r2" ]]; then
    paired=1
  elif [[ -f "$DATA_DIR/$acc.fastq" ]]; then
    paired=0 single="$DATA_DIR/$acc.fastq"
  elif [[ -f "$DATA_DIR/${acc}_0.fastq" ]]; then
    paired=0 single="$DATA_DIR/${acc}_0.fastq"
  else
    echo "== $acc: not staged in $DATA_DIR — skipping" >&2
    continue
  fi

  # The plain-FASTQ reference every tool must reproduce byte-for-byte. Two
  # normalizations make that byte count tool-independent:
  #   - the `+` separator line (record line 3, addressed by LINE NUMBER — a
  #     quality line may legitimately start with '+') is stripped of any
  #     repeated header: fqxv drops it by design (as SPRING/fqz_comp do), so
  #     the byte-stream tools must decode the same normalized content or the
  #     counts could never agree (SRA-era runs carry `+NAME length=N` lines);
  #   - paired datasets are concatenated R1+R2: fqxv stores them as one
  #     spot-interleaved archive and -Z streams the same records interleaved,
  #     so the byte COUNT matches cat(R1,R2) though the record order differs.
  plain="$WORK/$acc.fastq"
  if [[ "$paired" == 1 ]]; then
    fresh "$plain" "$r1" "$r2" ||
      { echo "== $acc: cat R1+R2, normalize '+'"; cat "$r1" "$r2" | sed '3~4s/^+.*/+/' > "$plain"; }
  else
    fresh "$plain" "$single" ||
      { echo "== $acc: normalize '+'"; sed '3~4s/^+.*/+/' "$single" > "$plain"; }
  fi
  expected="$(stat -Lc %s "$plain")"

  arc="$WORK/$acc.fqxv" maxarc="$WORK/$acc.max.fqxv"
  gz="$WORK/$acc.fastq.gz" zst="$WORK/$acc.fastq.zst"
  # fqxv_prep <archive> [--max] — the .fqxv must come from THIS binary (stale
  # archives on scratch predate the current container format), so the binary
  # itself is a dependency. Only the artifacts a selected tool needs are built.
  fqxv_prep() {
    local out="$1"; shift
    if [[ "$paired" == 1 ]]; then
      fresh "$out" "$r1" "$r2" "$FQXV_BIN" ||
        { echo "== $acc: fqxv compress $* (paired)"; "$FQXV_BIN" compress "$r1" "$r2" -o "$out" --force --quiet --threads "$NPROC" "$@"; }
    else
      fresh "$out" "$plain" "$FQXV_BIN" ||
        { echo "== $acc: fqxv compress $*"; "$FQXV_BIN" compress "$plain" -o "$out" --force --quiet --threads "$NPROC" "$@"; }
    fi
  }
  has_tool() { [[ " $TOOLS " == *" $1 "* ]]; }
  ! has_tool fqxv || fqxv_prep "$arc"
  ! has_tool fqxv-max || fqxv_prep "$maxarc" --max
  if has_tool gzip || has_tool pigz; then
    fresh "$gz" "$plain" || { echo "== $acc: pigz -6"; pigz -6 -p "$NPROC" -c "$plain" > "$gz"; }
  fi
  if has_tool zstd; then
    fresh "$zst" "$plain" || { echo "== $acc: zstd -19 --long=27"; zstd -19 --long=27 -T"$NPROC" -q -f -o "$zst" "$plain"; }
  fi

  for tool in $TOOLS; do
    case "$tool" in
      gzip) steps=1 comp="$gz" ;;   # no thread knob: the serial baseline, one row
      pigz) steps="$RUN_STEPS" comp="$gz" ;;
      zstd) steps="$RUN_STEPS" comp="$zst" ;;
      fqxv) steps="$RUN_STEPS" comp="$arc" ;;
      fqxv-max) steps="$RUN_STEPS" comp="$maxarc" ;;
      *) echo "error: unknown tool '$tool'" >&2; exit 1 ;;
    esac
    comp_bytes="$(stat -Lc %s "$comp")"
    for n in $steps; do
      case "$tool" in
        fqxv | fqxv-max) cmd="'$FQXV_BIN' decompress '$comp' -Z --quiet --threads $n" ;;
        gzip) cmd="gzip -dc '$comp'" ;;
        pigz) cmd="pigz -dc -p $n '$comp'" ;;
        zstd) cmd="zstd -dc --long=27 -T$n '$comp'" ;;
      esac
      # Correctness pass (also the first cache warm): the piped byte count
      # must equal the plain FASTQ size for every tool and thread count.
      bytes="$(set -o pipefail; bash -c "set -o pipefail; $cmd | wc -c")"
      if [[ "$bytes" == "$expected" ]]; then bytes_ok=yes; else
        bytes_ok=no; FAILED_BYTES=1
        echo "!! $acc/$tool -T$n: byte count $bytes != expected $expected" >&2
      fi
      time_cell "set -o pipefail; $cmd | wc -c > /dev/null"
      mbps="$(awk -v b="$bytes" -v s="$CELL_SECS" 'BEGIN{printf "%.1f", (s > 0) ? b/s/1e6 : 0}')"
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$acc" "$tool" "$n" "$CELL_SECS" "$mbps" "$bytes" "$bytes_ok" "$comp_bytes" "$REPS" >> "$OUT_TSV"
      echo "   $acc $tool t=$n ${CELL_SECS}s ${mbps} MB/s ($bytes_ok)"
    done
  done
done

echo "==> results: $OUT_TSV"
if [[ "$FAILED_BYTES" == 1 ]]; then
  echo "!! BYTE-COUNT MISMATCHES above — those cells are not comparable" >&2
  exit 1
fi
