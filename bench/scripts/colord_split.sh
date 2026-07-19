#!/usr/bin/env bash
# Re-derive CoLoRd's seq/qual split ADDITIVELY, and fqxv's the same way.
#   size(-q org)  = total lossless archive
#   size(-q none) = archive with quality discarded (Q1 for all bases)
#   quality       = org - none      <- marginal cost of real qualities
# Additive by construction: none + quality == org, exactly.
set -uo pipefail
D=$HOME/scratch/fqxv/data
O=$HOME/scratch/fqxv/colord_split
mkdir -p "$O"
CO="pixi exec -s colord -c bioconda -c conda-forge -- colord"
FQXV=/beevol/home/jhessel/devel/rnabioco/fqxv/target/release/fqxv
for DS in ecoli_ont ecoli_hifi; do
  IN=$D/$DS.fastq
  [ -f "$IN" ] || { echo "$DS: missing"; continue; }
  echo "##### $DS  ($(stat -c %s "$IN") bytes)"
  for Q in org none; do
    rm -f "$O/$DS.$Q.cord"
    /usr/bin/time -f "  colord -q $Q: %e s" $CO compress-ont -t 64 -q $Q "$IN" "$O/$DS.$Q.cord" >/dev/null 2>>"$O/$DS.$Q.err"
    echo "  colord -q $Q  = $(stat -c %s "$O/$DS.$Q.cord") bytes"
  done
  ORG=$(stat -c %s "$O/$DS.org.cord"); NONE=$(stat -c %s "$O/$DS.none.cord")
  echo "  => CoLoRd quality (org - none) = $((ORG-NONE))  | non-quality (none) = $NONE  | total = $ORG"
  echo "  => check: none + quality = $((NONE + ORG-NONE))  (must equal total $ORG)"
  rm -f "$O/$DS.fqxv"
  $FQXV compress "$IN" -o "$O/$DS.fqxv" -f --quiet 2>/dev/null
  echo "  fqxv total = $(stat -c %s "$O/$DS.fqxv") bytes"
  $FQXV inspect "$O/$DS.fqxv" 2>/dev/null | head -20
done
