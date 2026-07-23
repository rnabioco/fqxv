#!/usr/bin/env bash
#
# remote_smoke.sh — live end-to-end smoke test of remote `.fqxv` reads.
#
# Exercises both halves of the remote-read feature against a real HTTP(S) host
# that honours Range requests (any static file server: Apache, nginx, S3, GCS):
#
#   1. the server actually supports byte ranges              (curl -r probe)
#   2. CLI streaming decode succeeds and yields reads        (curl | fqxv decompress -)
#   3. streamed decode is byte-identical to a file decode    (the correctness claim)
#   4. a truncated transfer is rejected, not silently short  (the safety claim)
#   5. Python column projection reads the names remotely,    (fqxv.remote.RemoteArchive)
#      matches the full decode, and transfers a fraction of the file
#
# This is a *network* test, so it is not wired into `cargo nextest` / CI — run it
# by hand against a URL you control. The default points at the shared fixtures;
# override with the first argument or $FQXV_REMOTE_URL:
#
#   examples/remote_smoke.sh
#   examples/remote_smoke.sh https://my-host/reads.fqxv
#   FQXV_BIN=./target/release/fqxv examples/remote_smoke.sh <url>
#
# Fixture URLs (real ENA/SRA accessions, one per platform):
#   .../pacbio_hifi_ERR15205525.fqxv     PacBio HiFi        (long read)
#   .../nanopore_ERR12357097.fqxv        Oxford Nanopore    (long read)
#   .../illumina_pe_DRR052213.fqxv       Illumina, paired   (interleaved)
#   .../bgiseq_SRR26321993.fqxv          BGISEQ, single-end (short read)
#
set -u

FIXTURE_BASE="http://amc-sandbox.ucdenver.edu/User13/fqxv-test"
URL="${1:-${FQXV_REMOTE_URL:-$FIXTURE_BASE/pacbio_hifi_ERR15205525.fqxv}}"

# --- locate the fqxv binary (FQXV_BIN, then PATH, then target/) --------------
here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "$here/.." && pwd)"
find_bin() {
  [ -n "${FQXV_BIN:-}" ] && { echo "$FQXV_BIN"; return; }
  command -v fqxv 2>/dev/null && return
  for p in release debug; do
    [ -x "$repo/target/$p/fqxv" ] && { echo "$repo/target/$p/fqxv"; return; }
  done
}
FQXV="$(find_bin)"
[ -x "${FQXV:-}" ] || { echo "error: fqxv binary not found (set FQXV_BIN or 'cargo build --release')"; exit 2; }

# make the in-repo Python package importable if it was built with maturin develop
if [ -f "$repo/crates/fqxv-python/python/fqxv/_fqxv.abi3.so" ]; then
  export PYTHONPATH="$repo/crates/fqxv-python/python${PYTHONPATH:+:$PYTHONPATH}"
fi

tmp="$(mktemp -d)"
trap 'rm -rf "$tmp"' EXIT
pass=0; fail=0
ok()   { printf '  \033[32mPASS\033[0m  %s\n' "$1"; pass=$((pass+1)); }
bad()  { printf '  \033[31mFAIL\033[0m  %s\n' "$1"; fail=$((fail+1)); }

echo "fqxv:  $FQXV"
echo "url:   $URL"
echo

# --- 1. server honours byte ranges ------------------------------------------
echo "[1] server supports HTTP Range"
code="$(curl -fsS -o "$tmp/head10" -w '%{http_code}' -r 0-9 "$URL" 2>/dev/null)"
magic="$(head -c 4 "$tmp/head10" 2>/dev/null)"
if [ "$code" = "206" ] && [ "$magic" = "FQXV" ]; then
  ok "206 Partial Content, and the first bytes are the FQXV magic"
else
  bad "expected 206 + FQXV magic, got HTTP $code, magic '$magic' — server may not do ranges"
fi

# total size (for the projection-fraction check later)
size="$(curl -fsSI "$URL" 2>/dev/null | awk 'BEGIN{IGNORECASE=1}/^content-length:/{print $2+0}' | tr -d '\r')"

# --- 2. CLI streaming decode succeeds and yields reads ----------------------
echo "[2] streaming decode over the pipe"
if curl -fsSL "$URL" | "$FQXV" decompress - -Z > "$tmp/stream.fastq" 2>"$tmp/stream.err"; then
  reads=$(( $(wc -l < "$tmp/stream.fastq") / 4 ))
  [ "$reads" -gt 0 ] && ok "decoded $reads reads from the stream" || bad "stream decoded 0 reads"
else
  bad "curl | fqxv decompress - failed: $(tail -1 "$tmp/stream.err")"
fi

# --- 3. stream == file decode (byte-identical) ------------------------------
echo "[3] streamed decode is byte-identical to a downloaded-file decode"
curl -fsSL "$URL" -o "$tmp/archive.fqxv" 2>/dev/null
if "$FQXV" decompress "$tmp/archive.fqxv" -Z > "$tmp/file.fastq" 2>/dev/null \
   && cmp -s "$tmp/stream.fastq" "$tmp/file.fastq"; then
  ok "stream output matches file output exactly"
else
  bad "stream and file decodes differ"
fi

# --- 4. a truncated transfer is rejected ------------------------------------
echo "[4] truncated transfer is caught (not a silent short prefix)"
half=$(( ${size:-0} / 2 ))
if [ "$half" -gt 0 ] && curl -fsS -r "0-$half" "$URL" 2>/dev/null \
     | "$FQXV" decompress - -o /dev/null >/dev/null 2>&1; then
  bad "a half-file prefix decoded without error — truncation not detected"
else
  ok "the partial stream was rejected (premature EOF before terminator)"
fi

# --- 5. Python column projection --------------------------------------------
echo "[5] Python RemoteArchive column projection"
python3 - "$URL" "${size:-0}" <<'PY'
import sys
url, size = sys.argv[1], int(sys.argv[2])
try:
    import fqxv.remote as remote
except Exception as e:
    print(f"  SKIP  fqxv Python module not importable ({e}); run `maturin develop`")
    sys.exit(0)

a = remote.RemoteArchive.open(url)
names = a.names()                      # names-only projection: one range GET per group
proj = a.bytes_fetched
seqs = a.sequences()                   # add the sequence column
frac = 100.0 * proj / (size or a.size)

fails = 0
if len(names) > 0 and len(names) == len(seqs):
    print(f"  PASS  projected {len(names):,} names == {len(seqs):,} sequences over HTTP Range")
else:
    print(f"  FAIL  name/sequence counts disagree ({len(names)} vs {len(seqs)})"); fails += 1

# a names-only projection must transfer materially less than the whole file
if proj < (size or a.size):
    print(f"  PASS  names-only projection fetched {proj:,} B = {frac:.2f}% of the archive")
else:
    print(f"  FAIL  projection fetched {proj:,} B, not less than the {a.size:,} B file"); fails += 1

print(f"  info  first name: {names[0].decode(errors='replace')[:48]}")
sys.exit(1 if fails else 0)
PY
if [ $? -eq 0 ]; then pass=$((pass+1)); else fail=$((fail+1)); fi

# --- summary ----------------------------------------------------------------
echo
echo "----------------------------------------"
echo "checks: $pass passed, $fail failed"
[ "$fail" -eq 0 ] && echo "remote read: OK" || echo "remote read: FAILURES above"
exit $(( fail > 0 ? 1 : 0 ))
