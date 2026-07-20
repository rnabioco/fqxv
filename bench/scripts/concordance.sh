#!/usr/bin/env bash
# Phase-2 fidelity check: does lossy quality binning change variant calls?
#
# For each dataset with a reference genome (datasets.tsv `reference` column),
# align the ORIGINAL (lossless) reads and each fqxv `--quality-bin` output to the
# reference, call variants on each, and report how the binned calls concord with
# the lossless baseline (SNP + indel recall / precision). This is the fidelity
# half of the lossy tradeoff — read it next to the compression ratios from
# run_bench.sh / report.py.
#
# HEAVY (alignment + variant calling); NOT part of run_bench.sh. Run inside an
# srun/sbatch allocation, in the pixi env:
#   pixi run bash concordance.sh                # every dataset with a reference
#   pixi run bash concordance.sh ecoli_miseq    # just one
#
# Binned FASTQs are produced by the real codec (compress --quality-bin then
# decompress), so the concordance reflects exactly what fqxv would store. Env
# knobs mirror run_bench.sh: FQXV_DATA_DIR, FQXV_RESULTS_DIR, FQXV_THREADS,
# FQXV_BIN; FQXV_QUALITY_BINS overrides the bin sweep (default "bin8 bin4 bin2").
set -euo pipefail

DATA_DIR="${FQXV_DATA_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/data}"
RESULTS_DIR="${FQXV_RESULTS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/results}"
THREADS="${FQXV_THREADS:-$(nproc)}"
FQXV_BIN="${FQXV_BIN:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)/target/release/fqxv}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BINS="${FQXV_QUALITY_BINS:-bin8 bin4 bin2}"
REFDIR="$RESULTS_DIR/refs"
WORK="$RESULTS_DIR/concordance"
OUT="$RESULTS_DIR/concordance.tsv"
ONLY="${1:-}"

mkdir -p "$REFDIR" "$WORK"
echo -e "dataset\tbin\tbase_snps\tbin_snps\tsnp_shared\tsnp_recall\tsnp_precision\tbase_indels\tbin_indels\tindel_shared\tindel_recall\tindel_precision" > "$OUT"

# Count non-header records in a (plain) VCF; 0 if missing.
count_vcf() { [[ -f "$1" ]] && awk '!/^#/{n++} END{print n+0}' "$1" || echo 0; }

# Resolve a reference: `src` may be a URL (downloaded, gunzipped) or a local
# FASTA path. Lands it in REFDIR and builds the bwa + samtools indexes (cached).
# Echoes the local .fa path.
prepare_ref() {  # src
  local src="$1" base fa
  base="$(basename "$src")"; base="${base%.gz}"
  fa="$REFDIR/$base"
  if [[ ! -f "$fa" ]]; then
    if [[ -f "$src" ]]; then
      if [[ "$src" == *.gz ]]; then gzip -dc "$src" > "$fa"; else cp "$src" "$fa"; fi
    elif [[ "$src" == *.gz ]]; then
      curl -fsSL "$src" | gzip -dc > "$fa"
    else
      curl -fsSL "$src" > "$fa"
    fi
  fi
  [[ -f "$fa.fai" ]] || samtools faidx "$fa"
  [[ -f "$fa.bwt" ]] || bwa index "$fa" >/dev/null 2>&1
  echo "$fa"
}

# Align a FASTQ and call variants -> bgzipped, indexed VCF at "$pfx.vcf.gz".
call_variants() {  # fastq ref pfx
  local fq="$1" ref="$2" pfx="$3"
  bwa mem -t "$THREADS" "$ref" "$fq" 2>/dev/null \
    | samtools sort -@ "$THREADS" -o "$pfx.bam" - 2>/dev/null
  samtools index "$pfx.bam"
  bcftools mpileup -f "$ref" "$pfx.bam" 2>/dev/null \
    | bcftools call -mv -Oz -o "$pfx.vcf.gz" 2>/dev/null
  bcftools index "$pfx.vcf.gz"
}

# Concordance of a variant type between baseline and binned calls. Prints
# "base_total bin_total shared".
isec_counts() {  # base.vcf.gz bin.vcf.gz type(snps|indels) dir
  local base="$1" bin="$2" type="$3" dir="$4"
  rm -rf "$dir"; mkdir -p "$dir"
  bcftools view -v "$type" -Oz -o "$dir/base.vcf.gz" "$base" 2>/dev/null; bcftools index "$dir/base.vcf.gz"
  bcftools view -v "$type" -Oz -o "$dir/bin.vcf.gz" "$bin" 2>/dev/null; bcftools index "$dir/bin.vcf.gz"
  # isec writes 0000=base-only, 0001=bin-only, 0002=shared(from base).
  bcftools isec -p "$dir/isec" "$dir/base.vcf.gz" "$dir/bin.vcf.gz" >/dev/null 2>&1
  local ob ib sh
  ob="$(count_vcf "$dir/isec/0000.vcf")"
  ib="$(count_vcf "$dir/isec/0001.vcf")"
  sh="$(count_vcf "$dir/isec/0002.vcf")"
  echo "$((ob + sh)) $((ib + sh)) $sh"
}

ratio() { awk -v s="$1" -v t="$2" 'BEGIN{printf "%.4f", (t>0)?s/t:0}'; }

run_one() {  # acc label ref_src
  local acc="$1" label="$2" src="$3"
  if [[ "$src" == "-" || -z "$src" ]]; then echo "[skip] $label: no reference"; return; fi
  local r1="$DATA_DIR/${acc}_1.fastq"
  [[ -f "$r1" ]] || { echo "[skip] $label: $r1 missing (run fetch.sh)"; return; }
  [[ -x "$FQXV_BIN" ]] || { echo "[skip] $label: fqxv binary missing ($FQXV_BIN)"; return; }

  local ref w
  ref="$(prepare_ref "$src")"
  w="$WORK/$label"; mkdir -p "$w"

  echo "==> $label: lossless baseline calls"
  call_variants "$r1" "$ref" "$w/base"

  for b in $BINS; do
    echo "==> $label: $b"
    # The exact binned FASTQ the codec would store (compress then decompress).
    "$FQXV_BIN" compress "$r1" -o "$w/$b.fqxv" --force --quality-bin "$b" --threads "$THREADS" >/dev/null 2>&1
    "$FQXV_BIN" decompress "$w/$b.fqxv" -o "$w/$b.fastq" --force --threads "$THREADS" >/dev/null 2>&1
    call_variants "$w/$b.fastq" "$ref" "$w/$b"

    read -r bs bn sh < <(isec_counts "$w/base.vcf.gz" "$w/$b.vcf.gz" snps "$w/isec_snp_$b")
    read -r bi ni si < <(isec_counts "$w/base.vcf.gz" "$w/$b.vcf.gz" indels "$w/isec_indel_$b")
    local sr sp ir ip
    sr="$(ratio "$sh" "$bs")"; sp="$(ratio "$sh" "$bn")"
    ir="$(ratio "$si" "$bi")"; ip="$(ratio "$si" "$ni")"
    printf '  %-6s SNP recall=%s prec=%s (%s/%s shared)   INDEL recall=%s prec=%s (%s/%s)\n' \
      "$b" "$sr" "$sp" "$sh" "$bs" "$ir" "$ip" "$si" "$bi"
    echo -e "${label}\t${b}\t${bs}\t${bn}\t${sh}\t${sr}\t${sp}\t${bi}\t${ni}\t${si}\t${ir}\t${ip}" >> "$OUT"
  done
}

mapfile -t rows < <(grep -v '^#' "$HERE/../panels/datasets.tsv" | awk 'NF')
for row in "${rows[@]}"; do
  acc="$(awk '{print $1}' <<<"$row")"
  label="$(awk '{print $2}' <<<"$row")"
  src="$(awk '{print $7}' <<<"$row")"
  [[ -n "$ONLY" && "$label" != "$ONLY" ]] && continue
  run_one "$acc" "$label" "$src"
done

echo "==> wrote $OUT"
