#!/usr/bin/env bash
#
# Sample random SRA run accessions from ENA's portal API — the source of the
# fqxv robustness corpus (see corpus.sh). Deliberately random so the corpus
# spans the messy long tail of real FASTQ (odd read lengths, wide/narrow
# quality alphabets, Ns, empty reads, long names) that curated datasets miss.
#
# Usage:
#   bash sample_accessions.sh [-n N] [-s SEED] [-p PLATFORM] \
#        [--min-bases N] [--max-bases N] [--library-strategy S] [--pool N]
#
# Defaults: N=20, SEED=$RANDOM, PLATFORM=ILLUMINA, bases 20M-600M (modest files
# so a ~20-accession corpus fetches quickly and fits scratch).
#
# Prints one run_accession per line to stdout; a "# <all columns>" metadata
# line per pick to stderr so the caller can log platform/strategy/layout/bases.
#
# Adapted from rnabioco/sracha-rs validation/sample_accessions.sh.
set -uo pipefail

N=20
SEED=""
PLATFORM="ILLUMINA"
MIN_BASES=20000000
MAX_BASES=600000000
POOL=5000
LIBRARY_STRATEGY=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        -n)                   N="$2"; shift 2 ;;
        -s|--seed)            SEED="$2"; shift 2 ;;
        -p|--platform)        PLATFORM="$2"; shift 2 ;;
        --min-bases)          MIN_BASES="$2"; shift 2 ;;
        --max-bases)          MAX_BASES="$2"; shift 2 ;;
        --pool)               POOL="$2"; shift 2 ;;
        --library-strategy)   LIBRARY_STRATEGY="$2"; shift 2 ;;
        -h|--help)
            sed -n '3,16p' "$0" | sed 's/^# \{0,1\}//'
            exit 0
            ;;
        *) echo "unknown arg: $1" >&2; exit 2 ;;
    esac
done

[[ -z "$SEED" ]] && SEED="$RANDOM$RANDOM"

if [[ "$PLATFORM" == "all" ]]; then
    QUERY="base_count>=${MIN_BASES} AND base_count<=${MAX_BASES}"
else
    QUERY="instrument_platform=${PLATFORM} AND base_count>=${MIN_BASES} AND base_count<=${MAX_BASES}"
fi
[[ -n "$LIBRARY_STRATEGY" ]] && QUERY="${QUERY} AND library_strategy=${LIBRARY_STRATEGY}"

echo "# sampling N=${N} from ENA pool=${POOL} query=\"${QUERY}\" seed=${SEED}" >&2

TMP=$(mktemp)
SEED_SOURCE=$(mktemp)
trap 'rm -f "$TMP" "$SEED_SOURCE"' EXIT

# ENA portal API. POST with --data-urlencode to dodge shell-escape headaches.
HTTP_STATUS=$(curl -sS -o "$TMP" -w '%{http_code}' \
    -X POST "https://www.ebi.ac.uk/ena/portal/api/search" \
    --data-urlencode "result=read_run" \
    --data-urlencode "query=${QUERY}" \
    --data-urlencode "fields=run_accession,instrument_platform,instrument_model,library_strategy,library_layout,base_count" \
    --data-urlencode "limit=${POOL}" \
    --data-urlencode "format=tsv")

if [[ "$HTTP_STATUS" != "200" ]]; then
    echo "ERROR: ENA query failed with HTTP $HTTP_STATUS" >&2
    head -5 "$TMP" >&2
    exit 1
fi

LINES=$(wc -l < "$TMP")
if [[ "$LINES" -lt 2 ]]; then
    echo "ERROR: ENA returned empty result set" >&2
    exit 1
fi
echo "# pool=$((LINES - 1)) records returned" >&2

# Deterministic shuffle: seed a fixed-key AES-CTR keystream as shuf's random
# source so a given seed always yields the same corpus.
openssl enc -aes-256-ctr -pass "pass:${SEED}" -nosalt < /dev/zero 2>/dev/null \
    | head -c 1048576 > "$SEED_SOURCE"

tail -n +2 "$TMP" \
    | shuf -n "$N" --random-source="$SEED_SOURCE" \
    | awk -F'\t' '{print $1; print "# "$0 > "/dev/stderr"}'
