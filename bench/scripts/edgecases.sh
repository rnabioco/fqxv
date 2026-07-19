#!/usr/bin/env bash
#
# fqxv edge-case fuzzer — synthetic adversarial FASTQ that targets the documented
# invariants and risky code paths directly, where real SRA data (corpus.sh) can't
# reach. Two kinds of case:
#
#   rt      valid FASTQ that MUST round-trip losslessly (+ build identically at
#           --threads 1). Same content check as corpus.sh.
#   reject  malformed FASTQ that fqxv must handle *safely*: either a clean error
#           (nonzero exit, no panic/signal) or, if it chooses to accept the input,
#           a lossless round-trip. A crash (panic/SIGSEGV) or silent corruption
#           (accepted but round-trip differs) is a bug.
#
# Cases live in $EDGE_DIR/in, results in $EDGE_DIR/results.tsv, per-case logs in
# $EDGE_DIR/logs. Meant to run inside an srun/sbatch allocation (compute node).
#
#   pixi run bash corpus.sh build-digest      # once, for the fast round-trip hash
#   srun ... bash edgecases.sh                # gen + run (default)
#   bash edgecases.sh gen                     # just (re)generate inputs
#   bash edgecases.sh run                     # just run checks over existing inputs
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EDGE_DIR="${EDGE_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/edge}"
IN="$EDGE_DIR/in"; WORK="$EDGE_DIR/work"; LOGS="$EDGE_DIR/logs"
MANIFEST="$IN/manifest.tsv"
RESULTS="$EDGE_DIR/results.tsv"

FQXV_BIN="${FQXV_BIN:-${CARGO_TARGET_DIR:-$(cd "$HERE/.." && pwd)/target}/release/fqxv}"
FQDIGEST="${FQDIGEST:-${SCRATCH:-$HOME/scratch}/fqxv/tools/bin/fqdigest}"
FQDIGEST_SRC="$HERE/../tools/fqdigest.rs"
THREADS="${FQXV_THREADS:-8}"
MODES="${FQXV_MODES:-default max maxkeep bin8 bin4 bin2}"

mkdir -p "$IN" "$WORK" "$LOGS"

# Order-independent content digest via the fqdigest Rust tool, built on demand.
if [[ ! -x "$FQDIGEST" || "$FQDIGEST_SRC" -nt "$FQDIGEST" ]]; then
  mkdir -p "$(dirname "$FQDIGEST")"
  rustc -O --edition 2024 "$FQDIGEST_SRC" -o "$FQDIGEST"
fi
digest() {  # file...
  "$FQDIGEST" "$@"
}
digest_bin() {  # scheme file...
  local s="$1"; shift
  "$FQDIGEST" --bin "$s" "$@"
}

compress_args() { case "$1" in
  default) echo "" ;; max) echo "--max" ;;
  maxkeep) echo "--order any --keep-order" ;;   # exercise the permutation / order-restore path
  bin8) echo "--quality-bin bin8" ;; bin4) echo "--quality-bin bin4" ;; bin2) echo "--quality-bin bin2" ;;
  *) echo "" ;; esac; }
# The bin table a mode applies (for the binned-expected round-trip digest).
bin_of() { case "$1" in bin8|bin4|bin2) echo "$1" ;; *) echo "" ;; esac; }

# ------------------------------------------------------------------ generate
gen() {
  command -v python3 >/dev/null || { echo "need python3 to generate inputs" >&2; exit 1; }
  rm -f "$IN"/*.fastq "$MANIFEST" 2>/dev/null
  python3 - "$IN" "$MANIFEST" <<'PY'
import os, sys
IN, MAN = sys.argv[1], sys.argv[2]
man = []  # (id, kind, files_csv, desc)

def write(name, records, trailing_nl=True, eol="\n"):
    """records: list of (name, seq, plus, qual). Write raw so we control EOL/newline."""
    path = os.path.join(IN, name)
    parts = []
    for (n, s, p, q) in records:
        parts.append(eol.join([n, s, p, q]))
    data = eol.join(parts)
    if trailing_nl and records:
        data += eol
    with open(path, "w", newline="") as f:
        f.write(data)
    return name

def rec(name, seq, qual, plus="+"):
    return (name, seq, plus, qual)

# deterministic pseudo-sequences/qualities (no RNG → reproducible)
BASES = "ACGT"
def seq_of(n, alpha=BASES, phase=0):
    return "".join(alpha[(i+phase) % len(alpha)] for i in range(n))
def qual_of(n, lo=33, hi=73, phase=0):
    span = hi - lo + 1
    return "".join(chr(lo + ((i+phase) % span)) for i in range(n))

def simple(nrec, rlen=36, namebase="@r"):
    return [rec(f"{namebase}{i}", seq_of(rlen, phase=i), qual_of(rlen, phase=i)) for i in range(nrec)]

# --- rt cases: must round-trip losslessly ---
# 1 empty file (zero records)
with open(os.path.join(IN,"empty.fastq"),"w") as f: pass
man.append(("empty","rt","empty.fastq","zero records"))

# 2 single record
write("single.fastq", simple(1)); man.append(("single","rt","single.fastq","one record"))

# 3 zero-length reads (empty seq + empty qual) mixed with normal
write("zerolen.fastq", [rec("@z0","",""), rec("@z1", seq_of(20), qual_of(20)), rec("@z2","","")])
man.append(("zerolen","rt","zerolen.fastq","empty seq/qual records"))

# 4 non-ACGT: N, IUPAC ambiguity codes, '.' and '-'
write("nonacgt.fastq", [
    rec("@n0","ACGTNRYSWKMBDHVN","!"*16),
    rec("@n1","....----NNNN.ACGT","I"*17),
    rec("@n2", seq_of(30, alpha="ACGTN"), qual_of(30)),
]); man.append(("nonacgt","rt","nonacgt.fastq","N + IUPAC + . - bases (exception list)"))

# 5 all-N reads
write("alln.fastq", [rec(f"@x{i}","N"*40, "#"*40) for i in range(5)])
man.append(("alln","rt","alln.fastq","all-N sequences"))

# 6 full quality range Phred+33 33..126 (QMAX stress)
q = "".join(chr(c) for c in range(33,127))
write("qualrange.fastq", [rec("@q0", seq_of(len(q)), q)])
man.append(("qualrange","rt","qualrange.fastq","quality 33..126 (QMAX)"))

# 7 constant quality
write("qualconst.fastq", [rec(f"@c{i}", seq_of(50), "I"*50) for i in range(10)])
man.append(("qualconst","rt","qualconst.fastq","constant quality"))

# 8 wildly variable read lengths in one file
write("varlen.fastq", [rec(f"@v{i}", seq_of(L), qual_of(L)) for i,L in enumerate([1,2,5,50,500,5000,1,37,123])])
man.append(("varlen","rt","varlen.fastq","variable read lengths"))

# 9 one very long read (long-read style)
L=200000; write("longread.fastq", [rec("@long1", seq_of(L), qual_of(L))])
man.append(("longread","rt","longread.fastq","single 200kb read"))

# 10 very long / weird read names
write("longnames.fastq", [
    rec("@"+("TOKEN:"*80)+"end 1:N:0:ACGTACGT", seq_of(30), qual_of(30)),
    rec("@"+"x"*500, seq_of(30), qual_of(30)),
]); man.append(("longnames","rt","longnames.fastq","500-char / many-token names"))

# 11 tokenizer-hostile but PRESERVABLE names: no delimiters, empty comment, huge
# counter, multiple internal spaces, trailing space *within* a non-empty desc.
write("weirdnames.fastq", [
    rec("@nodelim", seq_of(20), qual_of(20)),
    rec("@read 999999999999999999", seq_of(20), qual_of(20)),   # counter > i64 -> non-numeric token
    rec("@", seq_of(20), qual_of(20)),                          # name is just '@'
    rec("@a  two  spaces", seq_of(20), qual_of(20)),            # multiple internal spaces (preserved)
    rec("@id desc ", seq_of(20), qual_of(20)),                  # trailing space but desc non-empty (preserved)
]); man.append(("weirdnames","rt","weirdnames.fastq","preservable delimiter/counter edge names"))

# 11b/c UNPRESERVABLE headers — filed as rnabioco/fqxv#49 (name/description split
# drops a trailing separator and rewrites tab->space). Tracked as xfail so the
# known bug doesn't mask new regressions; flip to `rt` when #49 is fixed.
write("trailspace.fastq", [rec("@a:b:c:d:e:f:g ", seq_of(20), qual_of(20))])
man.append(("trailspace","xfail","trailspace.fastq","#49 trailing space dropped"))
write("tabsep.fastq", [rec("@id\tdesc:1", seq_of(20), qual_of(20))])
man.append(("tabsep","xfail","tabsep.fastq","#49 tab separator -> space"))

# 12 homopolymer runs (context model)
write("homopolymer.fastq", [
    rec("@h0","A"*100,"I"*100), rec("@h1","G"*100,"I"*100),
    rec("@h2", ("AC"*50),"I"*100), rec("@h3", ("A"*50+"T"*50),"I"*100),
]); man.append(("homopolymer","rt","homopolymer.fastq","homopolymer/low-complexity"))

# 13 quality lines starting with '@' and '+' (FASTQ parser ambiguity)
write("qualat.fastq", [
    rec("@p0", seq_of(20), "@"+ "I"*19),
    rec("@p1", seq_of(20), "+"+ "I"*19),
    rec("@p2", seq_of(20), "@+@+@+@+@+@+@+@+@+@+"),
]); man.append(("qualat","rt","qualat.fastq","quality starts with @ or +"))

# 14 many duplicate identical records (reorder/dedup path)
write("dup.fastq", [rec("@d", "ACGTACGTACGTACGT", "IIIIIIIIIIIIIIII") for _ in range(200)])
man.append(("dup","rt","dup.fastq","200 identical reads"))

# 15 lowercase / soft-masked bases
write("lowercase.fastq", [rec("@lc0","acgtACGTacgtNnNn"*2, "I"*32)])
man.append(("lowercase","rt","lowercase.fastq","lowercase soft-masked bases"))

# 16 no trailing newline
write("notrailnl.fastq", simple(3), trailing_nl=False)
man.append(("notrailnl","rt","notrailnl.fastq","no trailing newline"))

# 17 CRLF line endings (digest strips CR on both sides — tests content modulo EOL)
write("crlf.fastq", simple(3), eol="\r\n")
man.append(("crlf","rt","crlf.fastq","CRLF line endings"))

# 18 '+' line repeats the name (fqxv drops it by design; digest excludes line 3)
write("plusname.fastq", [rec(f"@r{i}", seq_of(30), qual_of(30), plus=f"+r{i}") for i in range(5)])
man.append(("plusname","rt","plusname.fastq","+ line repeats name (normalized)"))

# 19 length-1 reads
write("onebase.fastq", [rec(f"@o{i}", BASES[i%4], "I") for i in range(8)])
man.append(("onebase","rt","onebase.fastq","single-base reads"))

# 20 paired: mates with different read lengths per spot
write("pairdifflen_1.fastq", [rec(f"@s{i}/1", seq_of(50), qual_of(50)) for i in range(20)])
write("pairdifflen_2.fastq", [rec(f"@s{i}/2", seq_of(75), qual_of(75)) for i in range(20)])
man.append(("pairdifflen","rt","pairdifflen_1.fastq,pairdifflen_2.fastq","paired, R1 50bp / R2 75bp"))

# --- reject cases: must error cleanly OR round-trip; never crash/corrupt ---
# R1 quality length != sequence length
write("bad_quallen.fastq", [rec("@b0", seq_of(30), qual_of(20))])
man.append(("bad_quallen","reject","bad_quallen.fastq","qual shorter than seq"))

# R2 truncated final record (missing quality line)
with open(os.path.join(IN,"truncated.fastq"),"w") as f:
    f.write("\n".join(["@t0", seq_of(20), "+", qual_of(20), "@t1", seq_of(20), "+"])+"\n")
man.append(("truncated","reject","truncated.fastq","truncated final record"))

# R3 not FASTQ at all (FASTA)
with open(os.path.join(IN,"fasta.fastq"),"w") as f:
    f.write(">seq1\nACGTACGTACGT\n>seq2\nTTTTGGGG\n")
man.append(("fasta","reject","fasta.fastq","FASTA, not FASTQ"))

# R4 record not starting with '@'
with open(os.path.join(IN,"noat.fastq"),"w") as f:
    f.write("\n".join(["r0", seq_of(20), "+", qual_of(20)])+"\n")
man.append(("noat","reject","noat.fastq","header missing @"))

# R5 paired mate-count mismatch (R1 has 4, R2 has 3)
write("paircount_1.fastq", simple(4, namebase="@m"))
write("paircount_2.fastq", simple(3, namebase="@m"))
man.append(("paircount","reject","paircount_1.fastq,paircount_2.fastq","paired count mismatch"))

# --- more structural corner cases (exercise grouping/reorder/tokenizer paths) ---

# 21/22 single-cell: 3- and 4-file per-spot interleave (10x-style: short barcode +
# cDNA + index reads). Round-trips as a multiset over all members.
def sc(n, lens, tag):
    return [write(f"{tag}_{k+1}.fastq",
                  [rec(f"@sc{i}", seq_of(L), qual_of(L)) for i in range(n)])
            for k, L in enumerate(lens)]
f3 = sc(30, [16, 90, 8], "sc3")
man.append(("singlecell3","rt",",".join(f3),"3-file per-spot interleave (bc/cDNA/index)"))
f4 = sc(30, [16, 90, 8, 8], "sc4")
man.append(("singlecell4","rt",",".join(f4),"4-file per-spot interleave (R1/R2/I1/I2)"))

# 23 single stream that is ALREADY interleaved paired (auto-detect should see g=2)
inter = []
for i in range(25):
    inter.append(rec(f"@spot{i}/1", seq_of(40), qual_of(40)))
    inter.append(rec(f"@spot{i}/2", seq_of(40, phase=1), qual_of(40, phase=1)))
write("interleaved.fastq", inter)
man.append(("interleaved","rt","interleaved.fastq","auto-detected interleaved paired"))

# 24 SRA-style counter names (delta/counter tokenizer + reorder name regeneration)
write("sracounter.fastq", [rec(f"@DRR000001.{i} {i} length=40", seq_of(40), qual_of(40)) for i in range(1,201)])
man.append(("sracounter","rt","sracounter.fastq","SRA counter names (delta/regen path)"))

# 25 many reads -> cross multiple compression blocks (block-boundary handling)
write("manyreads.fastq", simple(20000, rlen=60))
man.append(("manyreads","rt","manyreads.fastq","20k reads spanning blocks"))

# 26 reverse-complement-heavy: each read followed by its RC (reorder RC clustering)
COMP = {"A":"T","C":"G","G":"C","T":"A","N":"N"}
def rc(s): return "".join(COMP[c] for c in reversed(s))
rcrecs = []
for i in range(100):
    s = seq_of(60, phase=i)
    rcrecs.append(rec(f"@f{i}", s, qual_of(60)))
    rcrecs.append(rec(f"@r{i}", rc(s), qual_of(60)))
write("rcheavy.fastq", rcrecs)
man.append(("rcheavy","rt","rcheavy.fastq","reverse-complement pairs (reorder)"))

# 27 RNA 'U' + extended IUPAC + lowercase soft-masking in the sequence stream
write("rna_iupac.fastq", [
    rec("@u0","ACGUacguNnRYSWKMBDHVryswkm","I"*26),
    rec("@u1", seq_of(40, alpha="ACGUN"), qual_of(40)),
]); man.append(("rna_iupac","rt","rna_iupac.fastq","U/IUPAC/lowercase bases"))

# 28 N-dense reads (exception-list heavy)
write("ndense.fastq", [rec(f"@nd{i}", ("N"*8+"ACGT")*4, "I"*48) for i in range(20)])
man.append(("ndense","rt","ndense.fastq","N-dense (exception list)"))

# 29 blank line between records -> malformed, must reject
with open(os.path.join(IN,"blankline.fastq"),"w") as f:
    f.write("\n".join(["@b0", seq_of(20), "+", qual_of(20), "", "@b1", seq_of(20), "+", qual_of(20)])+"\n")
man.append(("blankline","reject","blankline.fastq","blank line between records"))

# 30 empty quality but non-empty sequence -> qual/seq length mismatch, reject
with open(os.path.join(IN,"zeroqual.fastq"),"w") as f:
    f.write("\n".join(["@z0", seq_of(20), "+", ""])+"\n")
man.append(("zeroqual","reject","zeroqual.fastq","empty qual, non-empty seq"))

with open(MAN,"w") as f:
    for row in man: f.write("\t".join(row)+"\n")
print(f"generated {len(man)} cases -> {IN}")
PY
}

# ------------------------------------------------------------------ run
crashed() {  # log -> 0 if a panic/signal signature is present
  grep -qiE "panicked|RUST_BACKTRACE|SIGSEGV|SIGABRT|core dumped|stack overflow" "$1"
}

run_rt() {  # id files_csv
  local id="$1" csv="$2"; IFS=',' read -ra rel <<< "$csv"
  local ins=(); local r; for r in "${rel[@]}"; do ins+=("$IN/$r"); done
  local ref; ref=$(digest "${ins[@]}")
  local mode
  for mode in $MODES; do
    local comp="$WORK/$id.$mode.fqxv" rt="$WORK/$id.$mode.rt" log="$LOGS/$id.$mode.log" args
    args=$(compress_args "$mode")
    # shellcheck disable=SC2086
    "$FQXV_BIN" compress "${ins[@]}" -o "$comp" --force --threads "$THREADS" $args > "$log" 2>&1
    local rc=$?
    if [[ $rc -ne 0 ]]; then
      crashed "$log" && { echo -e "$id\t$mode\tFAIL_CRASH\tcompress rc=$rc"; continue; }
      echo -e "$id\t$mode\tFAIL_COMPRESS\trc=$rc"; continue
    fi
    if ! "$FQXV_BIN" decompress "$comp" -o "$rt" --force --threads "$THREADS" >> "$log" 2>&1; then
      rc=$?; crashed "$log" && { echo -e "$id\t$mode\tFAIL_CRASH\tdecompress rc=$rc"; continue; }
      echo -e "$id\t$mode\tFAIL_DECOMPRESS\trc=$rc"; rm -f "$comp"; continue
    fi
    local got want
    got=$(digest "$rt")
    if [[ "$mode" == bin* ]]; then want=$(digest_bin "$mode" "${ins[@]}"); else want="$ref"; fi
    if [[ "$got" != "$want" ]]; then echo -e "$id\t$mode\tFAIL_RT\twant=$want got=$got"; continue; fi
    # determinism
    local comp1="$WORK/$id.$mode.t1.fqxv"
    # shellcheck disable=SC2086
    "$FQXV_BIN" compress "${ins[@]}" -o "$comp1" --force --threads 1 $args >> "$log" 2>&1
    if ! cmp -s "$comp" "$comp1"; then echo -e "$id\t$mode\tFAIL_DET\tt1 != t$THREADS"; rm -f "$rt"; continue; fi
    rm -f "$comp" "$comp1" "$rt"
    echo -e "$id\t$mode\tPASS\t-"
  done
}

run_reject() {  # id files_csv
  local id="$1" csv="$2"; IFS=',' read -ra rel <<< "$csv"
  local ins=(); local r; for r in "${rel[@]}"; do ins+=("$IN/$r"); done
  local comp="$WORK/$id.fqxv" rt="$WORK/$id.rt" log="$LOGS/$id.log"
  # shellcheck disable=SC2086
  "$FQXV_BIN" compress "${ins[@]}" -o "$comp" --force --threads "$THREADS" > "$log" 2>&1
  local rc=$?
  if [[ $rc -ge 128 ]] || crashed "$log"; then echo -e "$id\treject\tFAIL_CRASH\trc=$rc"; return; fi
  if [[ $rc -ne 0 ]]; then echo -e "$id\treject\tPASS_REJECT\tclean error rc=$rc"; rm -f "$comp"; return; fi
  # Accepted the malformed input: that's allowed only if it round-trips losslessly.
  if ! "$FQXV_BIN" decompress "$comp" -o "$rt" --force --threads "$THREADS" >> "$log" 2>&1; then
    crashed "$log" && { echo -e "$id\treject\tFAIL_CRASH\tdecompress"; return; }
    echo -e "$id\treject\tFAIL_DECOMPRESS\taccepted then failed to decode"; return
  fi
  local got want; got=$(digest "$rt"); want=$(digest "${ins[@]}")
  if [[ "$got" == "$want" ]]; then echo -e "$id\treject\tPASS_ACCEPT\taccepted + lossless"
  else echo -e "$id\treject\tFAIL_CORRUPT\taccepted but round-trip differs"; fi
  rm -f "$comp" "$rt"
}

# Known-failing case (a filed bug). Runs the default-mode round-trip and reports
# XFAIL when it fails (expected) or XPASS when it unexpectedly passes (bug fixed —
# promote it back to an `rt` case). Keeps filed bugs from masking new regressions.
run_xfail() {  # id files_csv note
  local id="$1" csv="$2" note="$3"; IFS=',' read -ra rel <<< "$csv"
  local ins=(); local r; for r in "${rel[@]}"; do ins+=("$IN/$r"); done
  local comp="$WORK/$id.fqxv" rt="$WORK/$id.rt" log="$LOGS/$id.log"
  if ! "$FQXV_BIN" compress "${ins[@]}" -o "$comp" --force --threads "$THREADS" > "$log" 2>&1; then
    crashed "$log" && { echo -e "$id\tdefault\tFAIL_CRASH\t$note"; return; }
    echo -e "$id\tdefault\tXFAIL\t$note"; return
  fi
  "$FQXV_BIN" decompress "$comp" -o "$rt" --force --threads "$THREADS" >> "$log" 2>&1 || { echo -e "$id\tdefault\tXFAIL\t$note"; return; }
  if [[ "$(digest "$rt")" == "$(digest "${ins[@]}")" ]]; then echo -e "$id\tdefault\tXPASS\t$note (now passes — promote to rt)"
  else echo -e "$id\tdefault\tXFAIL\t$note"; fi
  rm -f "$comp" "$rt"
}

run() {
  [[ -x "$FQXV_BIN" ]] || { echo "no fqxv binary at $FQXV_BIN" >&2; exit 1; }
  [[ -f "$MANIFEST" ]] || { echo "no manifest — run 'edgecases.sh gen' first" >&2; exit 1; }
  printf 'case\tmode\tstatus\tnote\n' > "$RESULTS"
  while IFS=$'\t' read -r id kind files desc; do
    [[ -z "$id" ]] && continue
    case "$kind" in
      rt)     run_rt "$id" "$files" ;;
      reject) run_reject "$id" "$files" ;;
      xfail)  run_xfail "$id" "$files" "$desc" ;;
    esac
  done < "$MANIFEST" | tee -a "$RESULTS"
  echo
  echo "=== tally ==="
  awk -F'\t' 'NR>1{c[$3]++} END{for(k in c)printf "  %-14s %d\n",k,c[k]}' "$RESULTS" | sort
  # XFAIL is a known/filed bug (expected); XPASS and every FAIL_* are actionable.
  echo "actionable (excludes expected XFAIL):"
  awk -F'\t' 'NR>1 && $3!~/^PASS/ && $3!="XFAIL"{printf "  %-14s %-8s %-14s %s\n",$1,$2,$3,$4}' "$RESULTS"
}

case "${1:-all}" in
  gen) gen ;;
  run) run ;;
  all) gen; run ;;
  *) echo "usage: edgecases.sh [gen|run|all]" >&2; exit 2 ;;
esac
