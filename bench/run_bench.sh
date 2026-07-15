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
# fqxv, fqxv9 (level 9), fqxv-reorder (--order any), fqxv-max (--max, i.e.
# `-l 9 --order any` — the advertised best-ratio preset), and the lossy quality
# points fqxv-bin8/bin4/bin2 all share one binary; the rest are external
# baselines. spring-illbin (`-q ill_bin`, Illumina 8-level) and spring-binary
# (`-q binary`, 2-level) are SPRING's lossy quality modes — the only field tools
# with Illumina-comparable binning, so they are the like-for-like lossy rivals to
# fqxv-bin8 and fqxv-bin2 (fqz_comp/fqzcomp5 have no Illumina binning mode).
ALL_TOOLS="fqxv fqxv9 fqxv-reorder fqxv-max fqxv-shuffle fqxv-bin8 fqxv-bin4 fqxv-bin2 fqxv-reorder-bin8 fqxv-reorder-bin4 fqxv-reorder-bin2 gzip zstd19 xz9 fqz_comp fqzcomp5 spring spring-illbin spring-binary colord"
TOOLS="${FQXV_TOOLS:-$ALL_TOOLS}"
# The fqxv binary (built with `cargo build --release`). Cargo honors
# CARGO_TARGET_DIR (set to $SCRATCH on this HPC), so the build lands there, NOT
# in ROOT/target — resolve the same location cargo actually wrote to, else the
# harness silently measures a stale ROOT/target/release leftover.
FQXV_BIN="${FQXV_BIN:-${CARGO_TARGET_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)/target}/release/fqxv}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="$RESULTS_DIR/work"

# From-source baselines (built by build_tools.sh) live here; fqzcomp5 needs its
# htscodecs shared lib on LD_LIBRARY_PATH. PgRC is built but intentionally NOT
# in ALL_TOOLS: it is a sequence-only read compressor (drops names, simplifies
# quality), so it is not comparable to full-FASTQ archivers — see README.
TOOLS_DIR="${FQXV_TOOLS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/tools}"
export PATH="$TOOLS_DIR/bin:$PATH"
export LD_LIBRARY_PATH="$TOOLS_DIR/lib:${LD_LIBRARY_PATH:-}"

# Rust content-digest tool (single O(n) streaming pass, bounded memory, no sort) —
# the record-multiset round-trip check. Built on demand from bench/fqdigest.rs.
FQDIGEST="${FQDIGEST:-$TOOLS_DIR/bin/fqdigest}"
FQDIGEST_SRC="$HERE/fqdigest.rs"
ensure_fqdigest() {
  if [[ ! -x "$FQDIGEST" || "$FQDIGEST_SRC" -nt "$FQDIGEST" ]]; then
    mkdir -p "$(dirname "$FQDIGEST")"
    rustc -O --edition 2021 "$FQDIGEST_SRC" -o "$FQDIGEST"
  fi
}
ensure_fqdigest

mkdir -p "$RESULTS_DIR" "$WORK"
# Execution mode (for the parallel Slurm-array harness — see submit_parallel.sh):
#   default        one node, all datasets x tools, sequential (clean throughput)
#   FQXV_PREP_ONLY=1        compute per-dataset meta + input digests only, into
#                          $RESULTS_DIR/prep/<label>.env, then exit (the shared,
#                          expensive record-sort done once instead of per cell)
#   FQXV_PART_TAG=<tag>     one array cell: run FQXV_ONLY_DATASET x FQXV_TOOLS,
#                          reuse the prep digests, append rows to a private part
#                          file $RESULTS_DIR/parts/results.<tag>.tsv (no header)
PREP_ONLY="${FQXV_PREP_ONLY:-0}"
PART_TAG="${FQXV_PART_TAG:-}"
if [[ -n "$PART_TAG" ]]; then
  mkdir -p "$RESULTS_DIR/parts"
  RESULTS="$RESULTS_DIR/parts/results.$PART_TAG.tsv"
  META="$RESULTS_DIR/parts/meta.$PART_TAG.tsv"
  : > "$RESULTS"
  : > "$META"
else
  RESULTS="$RESULTS_DIR/results.tsv"
  META="$RESULTS_DIR/meta.tsv"
fi

# GNU time for wall + peak RSS; fall back to bash timing (RSS unknown = -1).
GNU_TIME=""
for c in /usr/bin/time "$(command -v time || true)"; do
  if [[ -x "$c" ]] && "$c" -f '%e %M' true >/dev/null 2>&1; then GNU_TIME="$c"; break; fi
done

# measure CMD... -> sets MEAS_SECS, MEAS_RSS_KB, MEAS_RC (the command's exit
# code). Never returns non-zero itself, so a tool that fails — e.g. fqz_comp,
# which cannot parse long-read (Nanopore) FASTQ — is captured in MEAS_RC and
# recorded rather than aborting the whole run under `set -e`.
measure() {
  local tf; tf="$(mktemp)"; MEAS_RC=0
  if [[ -n "$GNU_TIME" ]]; then
    "$GNU_TIME" -o "$tf" -f '%e %M' "$@" || MEAS_RC=$?
    # On failure GNU time prepends a "Command exited with non-zero status N"
    # line before the "%e %M" line, so read the *last* line for the metrics.
    read -r MEAS_SECS MEAS_RSS_KB < <(tail -n1 "$tf") || { MEAS_SECS=-1; MEAS_RSS_KB=-1; }
  else
    local t0 t1; t0="$EPOCHREALTIME"; { "$@" || MEAS_RC=$?; }; t1="$EPOCHREALTIME"
    MEAS_SECS="$(awk -v a="$t0" -v b="$t1" 'BEGIN{printf "%.2f", b-a}')"
    MEAS_RSS_KB="-1"
  fi
  rm -f "$tf"
}

fastq_records() { echo $(( $(wc -l < "$1") / 4 )); }

# Order-independent record-multiset digest via fqdigest: hashes each
# (name, sequence[, quality]) record and sums the hashes, so the result is
# invariant to read reordering (SPRING, `fqxv --order any`) in one streaming pass.
# Verifies *content* losslessness. The `+` line (record line 3) is excluded:
# fqxv normalizes it, the one documented lossy-by-design deviation. mode=noqual
# drops quality (for lossy-quality tools), via `--no-qual`.
record_digest() {  # file mode(full|noqual|nonames)
  case "$2" in
    noqual)  "$FQDIGEST" --no-qual  "$1" ;;
    nonames) "$FQDIGEST" --no-names "$1" ;;
    *)       "$FQDIGEST"            "$1" ;;
  esac
}

# Like record_digest with quality, but first pass every quality byte through
# fqxv's bin table `scheme` (bin8|bin4|bin2) — i.e. the *expected* content of a
# correct lossy round-trip. fqdigest's `--bin` mirrors QualityBinning::apply.
record_digest_binned() {  # file scheme(bin8|bin4|bin2)
  "$FQDIGEST" --bin "$2" "$1"
}

# Per-base quality distortion of a lossy round-trip vs the original: mean absolute
# error, RMSE, and % of bases whose quality changed. Records are matched by name
# (order-independent, so it holds for read-reordering tools like SPRING). Prints
# "mae rmse pct"; "-1 -1 -1" if nothing matched. O(quality bytes) memory, like
# the digest sort, so only run for lossy tools.
qual_distortion() {  # orig rt  ->  "mae rmse pct"
  awk '
    BEGIN{ for(i=0;i<256;i++) ord[sprintf("%c",i)]=i }
    NR==FNR{ if(FNR%4==1) name=$0; else if(FNR%4==0) oq[name]=$0; next }
    { if(FNR%4==1) name=$0
      else if(FNR%4==0){ o=oq[name]; r=$0; m=length(r); if(length(o)<m)m=length(o)
        for(i=1;i<=m;i++){ d=ord[substr(r,i,1)]-ord[substr(o,i,1)]; if(d<0)d=-d; sa+=d; sq+=d*d; if(d>0)ch++; n++ } } }
    END{ if(n==0) print "-1 -1 -1"; else printf "%.4f %.4f %.4f\n", sa/n, sqrt(sq/n), 100.0*ch/n }
  ' "$1" "$2"
}

is_fqxv() { [[ "$1" == fqxv || "$1" == fqxv-* || "$1" == fqxv[0-9] ]]; }
# Lossy-quality tools (quality changed on purpose). Covers both the plain
# fqxv-bin* points and the fqxv-reorder-bin* combos (reorder + binning).
is_lossy() { [[ "$1" == fqxv*bin* || "$1" == spring-illbin || "$1" == spring-binary ]]; }
# The exact fqxv bin table a tool applies, for the binned-expected round-trip;
# `none` for tools whose internal table we don't assert (the SPRING rivals).
# Match on the bin suffix so fqxv-bin8 and fqxv-reorder-bin8 map alike.
bin_scheme() {
  case "$1" in
    *bin8) echo bin8 ;;
    *bin4) echo bin4 ;;
    *bin2) echo bin2 ;;
    *) echo none ;;
  esac
}
# Name-lossy tools renumber reads: the retained content is the seq+qual multiset,
# so verify that (`--no-names`), not the names. (`fqxv --order shuffle`.)
is_name_lossy() { [[ "$1" == *shuffle* ]]; }

# --- per-tool compress/decompress. Each sets COMP (compressed path) then RT. ---
compress() {  # tool input out_prefix
  local tool="$1" in="$2" pfx="$3"
  case "$tool" in
    fqxv)          COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --force --threads "$THREADS" ;;
    fqxv9)         COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --force -l 9 --threads "$THREADS" ;;
    fqxv-reorder)  COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --force --order any --threads "$THREADS" ;;
    # fqxv-max: the advertised best-ratio preset (`--max` == `-l 9 --order any`):
    # deepest sequence context AND read reordering together.
    fqxv-max)      COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --force --max --threads "$THREADS" ;;
    # fqxv-shuffle: best-ratio RENUMBER preset (`-l 9 --order shuffle`) — the
    # apples-to-apples point vs SPRING, which also renumbers/reorders. Reads come
    # back as a seq+qual multiset with fresh names (name-lossy); verified below
    # with a `--no-names` digest, as SPRING's own reordering is verified by a
    # (order-independent) multiset digest.
    fqxv-shuffle)  COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --force -l 9 --order shuffle --threads "$THREADS" ;;
    fqxv-bin8)     COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --force --quality-bin bin8 --threads "$THREADS" ;;
    fqxv-bin4)     COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --force --quality-bin bin4 --threads "$THREADS" ;;
    fqxv-bin2)     COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --force --quality-bin bin2 --threads "$THREADS" ;;
    # reorder + binning combined — the like-for-like rivals to SPRING's lossy
    # modes (spring-illbin vs fqxv-reorder-bin8, spring-binary vs -bin2), since
    # SPRING always reorders. The plain fqxv-bin* rows keep original order.
    fqxv-reorder-bin8) COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --force --order any --quality-bin bin8 --threads "$THREADS" ;;
    fqxv-reorder-bin4) COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --force --order any --quality-bin bin4 --threads "$THREADS" ;;
    fqxv-reorder-bin2) COMP="$pfx.fqxv"; measure "$FQXV_BIN" compress "$in" -o "$COMP" --force --order any --quality-bin bin2 --threads "$THREADS" ;;
    gzip)     COMP="$pfx.gz";  measure bash -c "pigz -p $THREADS -6 -c '$in' > '$COMP'" ;;
    zstd19)   COMP="$pfx.zst"; measure bash -c "zstd -19 --long=27 -T$THREADS -q -f -o '$COMP' '$in'" ;;
    xz9)      COMP="$pfx.xz";  measure bash -c "xz -9 -T$THREADS -c '$in' > '$COMP'" ;;
    fqz_comp) COMP="$pfx.fqz"; measure bash -c "fqz_comp < '$in' > '$COMP'" ;;
    fqzcomp5) COMP="$pfx.fqz5"; measure bash -c "fqzcomp5 < '$in' > '$COMP'" ;;
    spring)   COMP="$pfx.spring"; mkdir -p "$WORK/spring_c_$$"; measure spring -c -t "$THREADS" -i "$in" -o "$COMP" -w "$WORK/spring_c_$$/" ;;
    # CoLoRd long-read SOTA, lossless quality (`-q org`). Compresses sequence and
    # quality; the meaningful bar for our ONT streams.
    colord)   COMP="$pfx.colord"; measure bash -c "rm -f '$COMP'; colord compress-ont -t $THREADS -q org '$in' '$COMP'" ;;
    # spring-illbin: Illumina 8-level binning (like-for-like vs fqxv-bin8).
    spring-illbin) COMP="$pfx.spring"; mkdir -p "$WORK/spring_c_$$"; measure spring -c -t "$THREADS" -q ill_bin -i "$in" -o "$COMP" -w "$WORK/spring_c_$$/" ;;
    # spring-binary thr=25 high=37 low=15 mirrors fqxv-bin2 (q<25 -> 15, else 37).
    spring-binary) COMP="$pfx.spring"; mkdir -p "$WORK/spring_c_$$"; measure spring -c -t "$THREADS" -q binary 25 37 15 -i "$in" -o "$COMP" -w "$WORK/spring_c_$$/" ;;
    *) echo "unknown tool $tool" >&2; return 1 ;;
  esac
}
decompress() {  # tool comp out_rt
  local tool="$1" comp="$2" rt="$3"
  case "$tool" in
    fqxv|fqxv9|fqxv-reorder|fqxv-max|fqxv-shuffle|fqxv-bin8|fqxv-bin4|fqxv-bin2|fqxv-reorder-bin8|fqxv-reorder-bin4|fqxv-reorder-bin2) measure "$FQXV_BIN" decompress "$comp" -o "$rt" --force --threads "$THREADS" ;;
    gzip)     measure bash -c "pigz -d -p $THREADS -c '$comp' > '$rt'" ;;
    zstd19)   measure bash -c "zstd -d -q -f -o '$rt' '$comp'" ;;
    xz9)      measure bash -c "xz -d -T$THREADS -c '$comp' > '$rt'" ;;
    fqz_comp) measure bash -c "fqz_comp -d < '$comp' > '$rt'" ;;
    fqzcomp5) measure bash -c "fqzcomp5 -d < '$comp' > '$rt'" ;;
    spring|spring-illbin|spring-binary)   mkdir -p "$WORK/spring_d_$$"; measure spring -d -t "$THREADS" -i "$comp" -o "$rt" -w "$WORK/spring_d_$$/" ;;
    colord)   measure bash -c "rm -f '$rt'; colord decompress '$comp' '$rt'" ;;
  esac
}

# results.tsv columns: per-stream sizes are fqxv-only (-1 for other tools);
# rt_ok is now a *content* multiset check; deterministic is a 1-thread vs
# N-thread byte-identity check (fqxv only, else n/a). Part files carry no header
# (the merge step adds one); prep writes meta.tsv but not results.tsv.
if [[ -z "$PART_TAG" ]]; then
  echo -e "dataset\torig_bytes\tn_records\tn_bases" > "$META"
  [[ "$PREP_ONLY" == 1 ]] || echo -e "dataset\ttool\torig_bytes\tcomp_bytes\tratio\tc_secs\td_secs\tc_rss_kb\td_rss_kb\tnames_bytes\tseq_bytes\tqual_bytes\trt_ok\tdeterministic\tqual_mae\tqual_rmse\tqual_pct_changed" > "$RESULTS"
fi

mapfile -t rows < <(grep -v '^#' "$HERE/datasets.tsv" | awk 'NF')
for row in "${rows[@]}"; do
  acc="$(awk '{print $1}' <<<"$row")"
  label="$(awk '{print $2}' <<<"$row")"

  # Array cells process a single dataset.
  [[ -n "${FQXV_ONLY_DATASET:-}" && "$label" != "$FQXV_ONLY_DATASET" ]] && continue

  # Resolve input (R1, or R1+R2 concatenated). Single-end runs (e.g. Nanopore)
  # are written by `sracha --split split-3` as `${acc}.fastq` with no `_1`
  # suffix, so fall back to that name when the paired R1 is absent.
  r1="$DATA_DIR/${acc}_1.fastq"
  r2="$DATA_DIR/${acc}_2.fastq"
  [[ -f "$r1" ]] || r1="$DATA_DIR/${acc}.fastq"
  [[ -f "$r1" ]] || { echo "[skip] $label: $DATA_DIR/${acc}[_1].fastq missing (run fetch.sh)"; continue; }
  if [[ "$INPUT_MODE" == "cat" && -f "$r2" ]]; then
    in="$WORK/${label}.fastq"
    [[ -f "$in" ]] || cat "$r1" "$r2" > "$in"
  else
    in="$r1"
  fi

  # Meta + reference content digests. The digest is an expensive full-file
  # record sort, so in the parallel harness it is computed once by the prep
  # phase and every array cell just sources it. full = with quality (lossless
  # tools); noqual = names+bases only (lossy-quality tools).
  prep="$RESULTS_DIR/prep/$label.env"
  if [[ -z "$PREP_ONLY" || "$PREP_ONLY" != 1 ]] && [[ -f "$prep" ]]; then
    # shellcheck disable=SC1090
    source "$prep"
  else
    orig_bytes="$(stat -c %s "$in")"
    nrec="$(fastq_records "$in")"
    nbases="$(awk 'NR%4==2{b+=length($0)} END{print b+0}' "$in")"
    in_full="$(record_digest "$in" full)"
    in_noqual="$(record_digest "$in" noqual)"
    in_nonames="$(record_digest "$in" nonames)"
    if [[ "$PREP_ONLY" == 1 ]]; then
      mkdir -p "$RESULTS_DIR/prep"
      {
        printf "in=%q\n" "$in"
        printf "orig_bytes=%q\n" "$orig_bytes"
        printf "nrec=%q\n" "$nrec"
        printf "nbases=%q\n" "$nbases"
        printf "in_full=%q\n" "$in_full"
        printf "in_noqual=%q\n" "$in_noqual"
        printf "in_nonames=%q\n" "$in_nonames"
      } > "$prep"
    fi
  fi
  # meta.tsv is authored by the sequential run or the prep phase, never by cells.
  [[ -z "$PART_TAG" ]] && echo -e "${label}\t${orig_bytes}\t${nrec}\t${nbases}" >> "$META"
  echo "==> $label  ($(numfmt --to=iec "$orig_bytes"), $nrec reads, $(numfmt --to=iec "$nbases") bases)"
  # Prep phase stops here — no compression, just the shared digests + meta.
  [[ "$PREP_ONLY" == 1 ]] && continue

  for tool in $TOOLS; do
    # Map tool label -> binary name to probe availability.
    case "$tool" in
      gzip) bin=pigz ;; zstd19) bin=zstd ;; xz9) bin=xz ;; pgrc) bin=PgRC ;;
      fqxv|fqxv9|fqxv-reorder|fqxv-max|fqxv-shuffle|fqxv-bin8|fqxv-bin4|fqxv-bin2|fqxv-reorder-bin8|fqxv-reorder-bin4|fqxv-reorder-bin2) bin="$FQXV_BIN" ;;
      spring-illbin|spring-binary) bin=spring ;;
      *) bin="$tool" ;;
    esac
    if is_fqxv "$tool"; then
      [[ -x "$bin" ]] || { echo "  [miss] $tool ($bin — run: cargo build --release)"; continue; }
    else
      command -v "$bin" >/dev/null 2>&1 || { echo "  [miss] $tool ($bin)"; continue; }
    fi
    pfx="$WORK/${label}.${tool}"; rt="$WORK/${label}.${tool}.rt.fastq"
    rm -f "$pfx".* "$rt"

    compress "$tool" "$in" "$pfx"; c_secs="$MEAS_SECS"; c_rss="$MEAS_RSS_KB"; c_rc="$MEAS_RC"
    # A failed compressor may leave a partial/garbage file; don't report its size
    # as a real ratio. Record 0 bytes + rt=no so report.py ranks it last.
    if [[ "$c_rc" -ne 0 ]]; then
      echo "  [fail] $tool: compress exited $c_rc (recorded rt=no, continuing)"
      rm -f "$COMP"; comp_bytes=0
    else
      comp_bytes="$(stat -c %s "$COMP" 2>/dev/null || echo 0)"
    fi
    # Only attempt decompress when compress produced a real archive; a tool that
    # cannot handle this data (e.g. fqz_comp on long reads) is left as rt=no.
    d_secs=0; d_rss=-1
    if [[ "$c_rc" -eq 0 && "$comp_bytes" -gt 0 ]]; then
      decompress "$tool" "$COMP" "$rt"; d_secs="$MEAS_SECS"; d_rss="$MEAS_RSS_KB"
      [[ "$MEAS_RC" -ne 0 ]] && echo "  [fail] $tool: decompress exited $MEAS_RC"
    fi

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

    # Content round-trip + quality distortion. For fqxv lossy tools we know the
    # exact bin table, so the round-trip verifies the *full* binned content
    # (names + bases + input-quality-through-that-table). Competitor lossy tools
    # (spring-*) are checked on names + bases only — we do not assert their
    # internal table — but still get distortion metrics vs the original quality.
    rt_ok="no"; qmae=-1; qrmse=-1; qpct=-1
    if [[ -f "$rt" ]]; then
      if is_lossy "$tool"; then
        scheme="$(bin_scheme "$tool")"
        if [[ "$scheme" != none ]]; then
          [[ "$(record_digest "$rt" full)" == "$(record_digest_binned "$in" "$scheme")" ]] && rt_ok="yes"
        else
          [[ "$(record_digest "$rt" noqual)" == "$in_noqual" ]] && rt_ok="yes"
        fi
        read -r qmae qrmse qpct < <(qual_distortion "$in" "$rt")
      elif is_name_lossy "$tool"; then
        # Renumbered reads: seq+qual multiset must be preserved exactly.
        [[ "$(record_digest "$rt" nonames)" == "$in_nonames" ]] && rt_ok="yes"
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
        fqxv)         "$FQXV_BIN" compress "$in" -o "$det1" --force --threads 1 >/dev/null 2>&1 || true ;;
        fqxv9)        "$FQXV_BIN" compress "$in" -o "$det1" --force -l 9 --threads 1 >/dev/null 2>&1 || true ;;
        fqxv-reorder) "$FQXV_BIN" compress "$in" -o "$det1" --force --order any --threads 1 >/dev/null 2>&1 || true ;;
        fqxv-max)     "$FQXV_BIN" compress "$in" -o "$det1" --force --max --threads 1 >/dev/null 2>&1 || true ;;
        fqxv-shuffle) "$FQXV_BIN" compress "$in" -o "$det1" --force -l 9 --order shuffle --threads 1 >/dev/null 2>&1 || true ;;
        fqxv-bin8)    "$FQXV_BIN" compress "$in" -o "$det1" --force --quality-bin bin8 --threads 1 >/dev/null 2>&1 || true ;;
        fqxv-bin4)    "$FQXV_BIN" compress "$in" -o "$det1" --force --quality-bin bin4 --threads 1 >/dev/null 2>&1 || true ;;
        fqxv-bin2)    "$FQXV_BIN" compress "$in" -o "$det1" --force --quality-bin bin2 --threads 1 >/dev/null 2>&1 || true ;;
        fqxv-reorder-bin8) "$FQXV_BIN" compress "$in" -o "$det1" --force --order any --quality-bin bin8 --threads 1 >/dev/null 2>&1 || true ;;
        fqxv-reorder-bin4) "$FQXV_BIN" compress "$in" -o "$det1" --force --order any --quality-bin bin4 --threads 1 >/dev/null 2>&1 || true ;;
        fqxv-reorder-bin2) "$FQXV_BIN" compress "$in" -o "$det1" --force --order any --quality-bin bin2 --threads 1 >/dev/null 2>&1 || true ;;
      esac
      if [[ -f "$det1" ]] && cmp -s "$det1" "$COMP"; then deterministic="yes"; else deterministic="no"; fi
      rm -f "$det1"
    fi

    ratio="$(awk -v o="$orig_bytes" -v c="$comp_bytes" 'BEGIN{printf "%.3f", (c>0)?o/c:0}')"
    if is_lossy "$tool"; then
      printf '  %-13s ratio=%-6s c=%ss d=%ss rss=%sK rt=%s det=%s  Δq mae=%s rmse=%s chg=%s%%\n' \
        "$tool" "$ratio" "$c_secs" "$d_secs" "$c_rss" "$rt_ok" "$deterministic" "$qmae" "$qrmse" "$qpct"
    else
      printf '  %-13s ratio=%-6s c=%ss d=%ss rss=%sK rt=%s det=%s\n' \
        "$tool" "$ratio" "$c_secs" "$d_secs" "$c_rss" "$rt_ok" "$deterministic"
    fi
    echo -e "${label}\t${tool}\t${orig_bytes}\t${comp_bytes}\t${ratio}\t${c_secs}\t${d_secs}\t${c_rss}\t${d_rss}\t${names_b}\t${seq_b}\t${qual_b}\t${rt_ok}\t${deterministic}\t${qmae}\t${qrmse}\t${qpct}" >> "$RESULTS"
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
    measure "$FQXV_BIN" compress "$r1" "$r2" -o "$COMP" --force --threads "$THREADS"; c_secs="$MEAS_SECS"; c_rss="$MEAS_RSS_KB"
    comp_bytes="$(stat -c %s "$COMP" 2>/dev/null || echo 0)"
    # Restore both mates and concatenate to compare the multiset against R1+R2.
    measure "$FQXV_BIN" decompress "$COMP" --split "$rt" --force --threads "$THREADS"; d_secs="$MEAS_SECS"; d_rss="$MEAS_RSS_KB"
    rt_ok="no"
    # `--split` (PR #44) writes BGZF mates named <prefix>_R1.fastq.gz /
    # _R2.fastq.gz — decompress before digesting. record_digest sorts, so the
    # concat order across mates is irrelevant.
    shopt -s nullglob
    _parts=( "${rt}"_*.fastq.gz "${rt}"_*.fastq )
    shopt -u nullglob
    if [[ "${#_parts[@]}" -ge 2 ]]; then
      : > "${rt}.all"
      for _p in "${_parts[@]}"; do
        case "$_p" in *.gz) zcat "$_p" ;; *) cat "$_p" ;; esac >> "${rt}.all"
      done
      [[ "$(record_digest "${rt}.all" full)" == "$p_full" ]] && rt_ok="yes"
    fi
    mapfile -t _info < <("$FQXV_BIN" info "$COMP" --tsv 2>/dev/null || true)
    names_b=-1; seq_b=-1; qual_b=-1
    if [[ "${#_info[@]}" -ge 2 ]]; then
      IFS=$'\t' read -r -a _d <<<"${_info[1]}"; names_b="${_d[7]:--1}"; seq_b="${_d[8]:--1}"; qual_b="${_d[9]:--1}"
    fi
    ratio="$(awk -v o="$p_orig" -v c="$comp_bytes" 'BEGIN{printf "%.3f", (c>0)?o/c:0}')"
    printf '  %-13s ratio=%-6s c=%ss d=%ss rt=%s (R1+R2 interleaved)\n' "fqxv-paired" "$ratio" "$c_secs" "$d_secs" "$rt_ok"
    echo -e "${label}-paired\tfqxv-paired\t${p_orig}\t${comp_bytes}\t${ratio}\t${c_secs}\t${d_secs}\t${c_rss}\t${d_rss}\t${names_b}\t${seq_b}\t${qual_b}\t${rt_ok}\tn/a\t-1\t-1\t-1" >> "$RESULTS"
    rm -f "$pcat" "$pfx".* "${rt}"_* "${rt}.all"
  fi
done

echo "==> wrote $RESULTS"
