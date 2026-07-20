#!/usr/bin/env bash
# Single source of truth for the benchmark tool matrix, sourced by BOTH drivers:
# the sequential one (run_bench.sh) and the parallel Slurm one
# (submit_parallel.sh).
#
# Why this file exists: the two drivers used to keep independent tool lists and
# they drifted. The parallel driver — the one actually used for full runs —
# defaulted to a short-read set that omitted `fqxv-max` and `fqxv-shuffle` (the
# configurations that beat SPRING), both fair SPRING lossy baselines
# (`spring-illbin`, `spring-binary`), and the whole `fqxv-reorder-bin*` stack.
# A default parallel run therefore could not reproduce the project's own
# headline results, and the omission was invisible because the sequential
# driver's list was complete. Add a tool HERE and every run picks it up.
#
# Selection is by platform because most tools are platform-specific: SPRING is
# Illumina-only, CoLoRd is long-read-only, and the reorder codec auto-disables
# above `REORDER_MAX_MEAN_LEN` so `fqxv-reorder*` on long reads would just
# duplicate the plain `fqxv` rows on a separate node.

# --- Illumina / short read ------------------------------------------------
# The full field matrix. fqxv-bin8/4/2 are the lossy ladder; spring-illbin
# (8-level) and spring-binary (2-level) are the like-for-like lossy rivals to
# fqxv-bin8 and fqxv-bin2 — the only field tools with Illumina-comparable
# binning. fqxv-binont/binhifi are excluded: their cutpoints are calibrated for
# long-read quality distributions and mean nothing here.
FQXV_TOOLSET_SHORT="fqxv fqxv9 fqxv-max fqxv-shuffle fqxv-reorder \
fqxv-bin8 fqxv-bin4 fqxv-bin2 \
fqxv-reorder-bin8 fqxv-reorder-bin4 fqxv-reorder-bin2 \
gzip zstd19 xz9 fqz_comp fqzcomp5 spring spring-illbin spring-binary"

# --- Oxford Nanopore ------------------------------------------------------
# fqxv-binont carries the ONT-calibrated cutpoints; the generic bin8/4/2 ladder
# is kept alongside it to show what the Illumina-calibrated tables cost on ONT.
# fqz_comp cannot parse long reads and is expected to record `rt=no` — the row
# is kept deliberately, because dropping it would make an untested tool
# indistinguishable from an inapplicable one.
FQXV_TOOLSET_ONT="fqxv fqxv9 fqxv-max \
fqxv-bin8 fqxv-bin4 fqxv-bin2 fqxv-binont \
gzip zstd19 xz9 fqz_comp colord colord-lossy"

# --- PacBio HiFi ----------------------------------------------------------
# Same shape as ONT (fqz_comp included for the same reason). fqxv-binont is kept
# on HiFi as the aggressive-binning point; fqxv-binhifi is the calibrated one.
FQXV_TOOLSET_HIFI="fqxv fqxv9 fqxv-max \
fqxv-bin8 fqxv-bin4 fqxv-bin2 fqxv-binhifi fqxv-binont \
gzip zstd19 xz9 fqz_comp colord colord-lossy"

# Echo the tool set for a datasets.tsv platform column value.
fqxv_toolset_for_platform() {
  case "$1" in
    MinION | GridION | PromethION | *[Nn]anopore* | ONT) echo "$FQXV_TOOLSET_ONT" ;;
    SequelII | Sequel* | Revio | *[Hh]iFi* | PacBio*) echo "$FQXV_TOOLSET_HIFI" ;;
    # MGI/BGISEQ is short-read, so it takes the full short-read field matrix.
    # Listed explicitly rather than left to the `*)` fallback: it reads as a
    # deliberate choice, and it is the one short-read platform where a field
    # tool may legitimately fail — SPRING assumes fixed-length Illumina reads,
    # so a variable-length MGI run can record `rt=no`. That row is kept for the
    # same reason ONT keeps its fqz_comp row: dropping it would make an
    # untested tool indistinguishable from an inapplicable one.
    BGISEQ* | MGISEQ* | DNBSEQ* | *[Mm][Gg][Ii]*) echo "$FQXV_TOOLSET_SHORT" ;;
    *) echo "$FQXV_TOOLSET_SHORT" ;;
  esac
}

# Union of every set, order-preserving and de-duplicated. Used where a
# platform-independent list is needed (documentation, `build_tools.sh` probes).
# PgRC is intentionally absent from all sets: it is a sequence-only read
# compressor (drops names, simplifies quality), so it is not comparable to full
# FASTQ archivers — see bench/README.md.
FQXV_TOOLSET_ALL="$(
  printf '%s\n' $FQXV_TOOLSET_SHORT $FQXV_TOOLSET_ONT $FQXV_TOOLSET_HIFI |
    awk '!seen[$0]++' | tr '\n' ' ' | sed 's/ $//'
)"
