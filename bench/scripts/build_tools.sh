#!/usr/bin/env bash
# Build the from-source baseline compressors (PgRC2, fqzcomp5) that aren't in
# bioconda, using the pixi toolchain. Binaries land in $FQXV_TOOLS_DIR/bin,
# which run_bench.sh prepends to PATH.
#
#   pixi run build-tools          # or: pixi run bash build_tools.sh
set -euo pipefail

TOOLS_DIR="${FQXV_TOOLS_DIR:-${SCRATCH:-$HOME/scratch}/fqxv/tools}"
BIN="$TOOLS_DIR/bin"
SRC="$TOOLS_DIR/src"
JOBS="$(nproc)"
mkdir -p "$BIN" "$SRC"

echo "==> build tools into $BIN (CC=${CC:-cc}, CXX=${CXX:-c++})"

# --- PgRC2 (kowallus/PgRC) : CMake, produces the `PgRC` binary ----------------
if [[ -x "$BIN/PgRC" ]]; then
  echo "  [skip] PgRC (already built)"
else
  echo "  [build] PgRC"
  cd "$SRC"; rm -rf PgRC
  git clone --depth 1 https://github.com/kowallus/PgRC.git
  cd PgRC
  cmake -DCMAKE_BUILD_TYPE=Release . >/dev/null
  make -j"$JOBS" PgRC
  cp -f PgRC "$BIN/"
fi

# --- fqzcomp5 (jkbonfield/fqzcomp5) : Make, bundles htscodecs submodule --------
if [[ -x "$BIN/fqzcomp5" ]]; then
  echo "  [skip] fqzcomp5 (already built)"
else
  echo "  [build] fqzcomp5"
  cd "$SRC"; rm -rf fqzcomp5
  git clone --recursive --depth 1 https://github.com/jkbonfield/fqzcomp5.git
  cd fqzcomp5
  # htscodecs ships as an autotools project; bootstrap + configure it first.
  if [[ -d htscodecs ]]; then
    ( cd htscodecs && autoreconf -i >/dev/null 2>&1 && ./configure >/dev/null && make -j"$JOBS" ) || {
      echo "  [warn] htscodecs autotools build failed; trying fqzcomp5 make anyway" >&2
    }
  fi
  make -j"$JOBS" || make
  cp -f fqzcomp5 "$BIN/"
  # fqzcomp5 dynamically links libhtscodecs; stage the shared lib so
  # run_bench.sh can find it via LD_LIBRARY_PATH=$TOOLS_DIR/lib.
  mkdir -p "$TOOLS_DIR/lib"
  find . -name 'libhtscodecs.so*' -exec cp -a {} "$TOOLS_DIR/lib/" \;
fi

echo "==> done:"
ls -l "$BIN"
"$BIN/PgRC" 2>&1 | head -3 || true
"$BIN/fqzcomp5" 2>&1 | head -3 || true
