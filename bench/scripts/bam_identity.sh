#!/usr/bin/env bash
# Alignment-level round-trip proof: do the BAMs match before vs after fqxv?
#
# The strongest end-to-end fidelity check we can run: push reads all the way
# through a real aligner (bwa mem) and compare the resulting alignments, rather
# than trusting fqxv's internal round-trip digest. Comparison uses `bamcmp`
# (bench/bamcmp.rs, built on demand) for order-independent multiset digests in a
# single streaming pass — no `samtools sort` and no `sort | md5sum`. Each record
# is hashed and the hashes summed, so the digest survives read reordering AND
# renaming while preserving multiplicity:
#
#   content : whole record        (QNAME + FLAG..QUAL + tags) — order/name/all
#   body    : record minus QNAME  (FLAG..QUAL + tags)         — invariant to rename
#   place   : FLAG..SEQ only       (no QNAME, no QUAL)         — invariant to QUAL too
#   coord   : md5 of the *coordinate-sorted* record stream — the actual file a
#             user gets from `samtools sort` (order-SENSITIVE, so it exposes
#             tie-order changes among equal-position reads).
#
# Modes:
#   lossless        default codec. Reads are byte-identical, so ALL of content/
#                   body/place/coord must equal the original -> byte-identical BAM.
#   reorder-any     --order any: permits reordering. For SRA-style names the
#                   library restores original order automatically (order_changed
#                   should read no) -> BAM identical.
#   reorder-shuffle --order shuffle: renumber+drop order where names are a saving;
#                   falls back to order-preserving otherwise.
#   reorder-forced  a genuine whole-record permutation of the decompressed reads
#                   (shuf). fqxv itself does not reorder real SRA output (the
#                   modes above read order_changed=no), so this isolates the
#                   ALIGNER: bwa mem is order-sensitive — it realigns a small
#                   fraction of reads differently when their file position
#                   changes (deterministic; not threading, not fixed by -K). The
#                   note column reports that fraction (~1% on this set), and
#                   fqdigest confirms the read multiset is identical, so the
#                   difference is bwa's, not fqxv's. Takeaway: preserving read
#                   order (fqxv's default, and its reorder modes on real data) is
#                   what guarantees a reproducible BAM.
#   bin8/bin4/bin2  lossy --quality-bin. Only quality changes: `place` must match
#                   (reads don't move) while `body`/`content` differ by QUAL. Per-
#                   base Phred distortion is reported (bamcmp qualdelta).
#
# Dataset needs a reference genome in datasets.tsv (col 7). Defaults to the small
# E. coli MiSeq set (runs in a few minutes); any labelled dataset with a
# reference works.
#
#   pixi run bash bam_identity.sh                 # default dataset (ecoli_miseq)
#   pixi run bash bam_identity.sh ecoli_miseq     # explicit
#
# Env: FQXV_DATA_DIR, FQXV_RESULTS_DIR, FQXV_THREADS, FQXV_BIN, FQXV_BAMCMP;
# FQXV_QUALITY_BINS overrides the bin sweep; FQXV_READS subsamples to N reads
# (0 = all); FQXV_BAM_DATASET sets the default label.
set -euo pipefail

DATA_DIR="${FQXV_DATA_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/data}"
RESULTS_DIR="${FQXV_RESULTS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/results}"
THREADS="${FQXV_THREADS:-$(nproc)}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
FQXV_BIN="${FQXV_BIN:-$ROOT/target/release/fqxv}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BAMCMP="${FQXV_BAMCMP:-${SCRATCH:-$HOME/scratch}/fqxv/tools/bin/bamcmp}"
BINS="${FQXV_QUALITY_BINS:-bin8 bin4 bin2}"
READS="${FQXV_READS:-4000000}"          # subsample cap; 0 = all reads
REFDIR="$RESULTS_DIR/refs"
WORK="$RESULTS_DIR/bam_identity"
OUT="$RESULTS_DIR/bam_identity.tsv"
DATASET="${1:-${FQXV_BAM_DATASET:-ecoli_miseq}}"
TMP="$WORK/tmp"
mkdir -p "$REFDIR" "$WORK" "$TMP"

FQDIGEST="${FQDIGEST:-${SCRATCH:-$HOME/scratch}/fqxv/tools/bin/fqdigest}"

# Build the digest tools on demand (single-file rustc, like the corpus check).
build_tool() {  # src bin
  [[ -x "$2" ]] && return
  command -v rustc >/dev/null || { echo "rustc not on PATH (need the rust toolchain)" >&2; exit 1; }
  mkdir -p "$(dirname "$2")"
  echo "==> building $(basename "$2") -> $2"
  rustc -O --edition 2021 "$1" -o "$2"
}
build_tool "$HERE/../tools/bamcmp.rs" "$BAMCMP"
build_tool "$HERE/../tools/fqdigest.rs" "$FQDIGEST"

prepare_ref() {  # src -> echoes local .fa path
  local src="$1" base fa
  base="$(basename "$src")"; base="${base%.gz}"
  fa="$REFDIR/$base"
  if [[ ! -f "$fa" ]]; then
    if [[ -f "$src" ]]; then
      if [[ "$src" == *.gz ]]; then gzip -dc "$src" > "$fa"; else cp "$src" "$fa"; fi
    elif [[ "$src" == *.gz ]]; then curl -fsSL "$src" | gzip -dc > "$fa"
    else curl -fsSL "$src" > "$fa"; fi
  fi
  [[ -f "$fa.fai" ]] || samtools faidx "$fa"
  [[ -f "$fa.bwt" ]] || bwa index "$fa" >/dev/null 2>&1
  echo "$fa"
}

align() {  # fastq ref out.bam
  bwa mem -t "$THREADS" "$2" "$1" 2>/dev/null | samtools view -b -@ "$THREADS" -o "$3" -
}

# Order-independent content/body/place digests; sets DC/DB/DP.
read_digests() {  # in.bam
  local o; o="$(samtools view "$1" 2>/dev/null | "$BAMCMP" digest)"
  DC="$(awk '$1=="content"{print $2}' <<<"$o")"
  DB="$(awk '$1=="body"{print $2}' <<<"$o")"
  DP="$(awk '$1=="place"{print $2}' <<<"$o")"
}

# md5 of the coordinate-sorted record stream (order-sensitive; the real file).
coord_md5() {  # in.bam tag
  samtools sort -@ "$THREADS" -T "$TMP/c.$2.$$" -O bam "$1" 2>/dev/null \
    | samtools view - 2>/dev/null | md5sum | awk '{print $1}'
}

mapped_frac() {  # in.bam
  samtools flagstat "$1" 2>/dev/null \
    | awk '/ in total/{t=$1} / mapped \(/{m=$1} END{printf "%.6f", (t>0)?m/t:0}'
}

order_changed() {  # orig.fastq other.fastq -> yes/no
  if cmp -s <(awk 'NR%4==1' "$1") <(awk 'NR%4==1' "$2"); then echo no; else echo yes; fi
}

yn() { [[ "$1" == "$2" ]] && echo YES || echo NO; }

# --- resolve dataset ----------------------------------------------------------
row="$(grep -v '^#' "$HERE/../panels/datasets.tsv" | awk -v L="$DATASET" '$2==L{print; exit}')"
[[ -n "$row" ]] || { echo "no dataset labelled '$DATASET' in datasets.tsv" >&2; exit 1; }
acc="$(awk '{print $1}' <<<"$row")"; src="$(awk '{print $7}' <<<"$row")"
[[ "$src" != "-" && -n "$src" ]] || { echo "dataset '$DATASET' has no reference (col 7)"; exit 1; }
[[ -x "$FQXV_BIN" ]] || { echo "fqxv binary missing: $FQXV_BIN" >&2; exit 1; }
r1="$DATA_DIR/${acc}_1.fastq"; [[ -f "$r1" ]] || r1="$DATA_DIR/${acc}.fastq"
[[ -f "$r1" ]] || { echo "reads missing for $acc (run fetch.sh $acc)" >&2; exit 1; }

w="$WORK/$DATASET"; mkdir -p "$w"
echo "==> dataset=$DATASET acc=$acc  bwa threads=$THREADS"
ref="$(prepare_ref "$src")"; echo "==> reference: $ref"

orig="$w/orig.fastq"
if [[ "$READS" -gt 0 ]]; then head -n "$((READS*4))" "$r1" > "$orig"; else cp "$r1" "$orig"; fi
echo "==> using $(( $(wc -l < "$orig") / 4 )) reads"

echo "==> aligning ORIGINAL reads"
align "$orig" "$ref" "$w/orig.bam"
read_digests "$w/orig.bam"; base_c="$DC"; base_b="$DB"; base_p="$DP"
base_coord="$(coord_md5 "$w/orig.bam" orig)"
base_map="$(mapped_frac "$w/orig.bam")"
echo "    content=$base_c body=$base_b place=$base_p mapped=$base_map"

: > "$OUT"
echo -e "dataset\tmode\torder_changed\tcontent_id\tbody_id\tplace_id\tcoord_id\tqual_mean_abs\tqual_rmse\tqual_max\tqual_pct_changed\tmapped_frac\tnote" >> "$OUT"

NOTE="-"  # per-mode annotation (set before evaluate; reset each call)

# Digest an already-aligned variant ($w/$mode.bam), compare to base, append row.
evaluate() {  # mode fastq
  local mode="$1" fq="$2"
  read_digests "$w/$mode.bam"
  local coord map ordc c_id b_id p_id co_id
  coord="$(coord_md5 "$w/$mode.bam" "$mode")"
  map="$(mapped_frac "$w/$mode.bam")"; ordc="$(order_changed "$orig" "$fq")"
  c_id="$(yn "$DC" "$base_c")"; b_id="$(yn "$DB" "$base_b")"
  p_id="$(yn "$DP" "$base_p")"; co_id="$(yn "$coord" "$base_coord")"
  local qma=- qrm=- qmx=- qpct=-
  if [[ "$mode" == bin* ]]; then
    local qn qch
    read -r qn qch qma qrm qmx < <("$BAMCMP" qualdelta "$orig" "$fq")
    qpct="$(awk -v c="$qch" -v n="$qn" 'BEGIN{printf "%.4f",(n>0)?100*c/n:0}')"
  fi
  printf '    order_changed=%s  content=%s body=%s place=%s coord=%s  mapped=%s\n' \
    "$ordc" "$c_id" "$b_id" "$p_id" "$co_id" "$map"
  [[ "$mode" == bin* ]] && printf '    qualΔ  mean|abs|=%s  rmse=%s  max=%s  changed=%s%%\n' \
    "$qma" "$qrm" "$qmx" "$qpct"
  [[ "$NOTE" != "-" ]] && printf '    note: %s\n' "$NOTE"
  echo -e "${DATASET}\t${mode}\t${ordc}\t${c_id}\t${b_id}\t${p_id}\t${co_id}\t${qma}\t${qrm}\t${qmx}\t${qpct}\t${map}\t${NOTE}" >> "$OUT"
  NOTE="-"
}

# For a forced whole-record permutation, bwa itself is order-sensitive (it
# realigns some reads differently at batch boundaries — deterministic, but
# order-dependent). Since reorder-forced keeps read names, pair primary
# alignments by QNAME and report the fraction that differ, so the row is
# quantified rather than a bare NO. Also confirm the read multiset is identical
# (fqdigest) — proving the difference is the aligner's, not fqxv's.
aligner_order_effect() {  # returns a NOTE string
  local dorig dperm diff n st="$TMP/oe.$$"; mkdir -p "$st"
  dorig="$("$FQDIGEST" "$orig")"; dperm="$("$FQDIGEST" "$w/reorder-forced.fastq")"
  local same="reads-differ"; [[ "$dorig" == "$dperm" ]] && same="reads-identical"
  samtools view -F 0x900 "$w/orig.bam" 2>/dev/null | cut -f1-6 \
    | LC_ALL=C sort -k1,1 -S 1G -T "$st" > "$st/a"
  samtools view -F 0x900 "$w/reorder-forced.bam" 2>/dev/null | cut -f1-6 \
    | LC_ALL=C sort -k1,1 -S 1G -T "$st" > "$st/b"
  read -r diff n < <(paste "$st/a" "$st/b" | awk -F'\t' \
    '{n++; if($2!=$8||$3!=$9||$4!=$10||$5!=$11||$6!=$12) d++} END{printf "%d %d\n", d+0, n+0}')
  rm -rf "$st"
  awk -v s="$same" -v d="$diff" -v n="$n" \
    'BEGIN{printf "%s; bwa-order-effect=%d/%d reads (%.4f%%)", s, d, n, (n>0)?100*d/n:0}'
}

# fqxv round-trip -> $w/$mode.fastq, then evaluate.
roundtrip() {  # mode compress-args...
  local mode="$1"; shift
  echo "==> $mode"
  "$FQXV_BIN" compress "$orig" -o "$w/$mode.fqxv" --force "$@" --threads "$THREADS" >/dev/null 2>&1
  "$FQXV_BIN" decompress "$w/$mode.fqxv" -o "$w/$mode.fastq" --force --threads "$THREADS" >/dev/null 2>&1
  align "$w/$mode.fastq" "$ref" "$w/$mode.bam"
  evaluate "$mode" "$w/$mode.fastq"
}

roundtrip lossless
roundtrip reorder-any     --order any
roundtrip reorder-shuffle --order shuffle

# reorder-forced: a real whole-record permutation of the lossless reads (records
# are 4-line groups; paste/tr keeps them intact). Deterministic random source so
# the run is reproducible.
echo "==> reorder-forced (shuf of decompressed reads)"
paste - - - - < "$w/lossless.fastq" \
  | shuf --random-source=<(yes fqxv-seed-1337) \
  | tr '\t' '\n' > "$w/reorder-forced.fastq"
align "$w/reorder-forced.fastq" "$ref" "$w/reorder-forced.bam"
NOTE="$(aligner_order_effect)"   # quantifies bwa's own order-sensitivity
evaluate reorder-forced "$w/reorder-forced.fastq"

for b in $BINS; do roundtrip "$b" --quality-bin "$b"; done

echo; echo "==> wrote $OUT"; echo
column -t -s $'\t' "$OUT"
