#!/usr/bin/env bash
# Build ped-sim (GPL-3.0; run as a binary, never vendored) at a pinned tag into .pedsim/.
# Skips the build when .pedsim/ped-sim already exists at the pinned commit.
set -euo pipefail
TAG=v1.4.7
COMMIT=44d0c7abf6206479d5d6df926ecca086a8bbee32
DIR="$PIXI_PROJECT_ROOT/.pedsim"
if [[ -x "$DIR/ped-sim" && "$(git -C "$DIR" rev-parse HEAD 2>/dev/null)" == "$COMMIT" ]]; then
  exit 0
fi
rm -rf "$DIR"
git clone --quiet --depth 1 --branch "$TAG" https://github.com/williamslab/ped-sim.git "$DIR"
got="$(git -C "$DIR" rev-parse HEAD)"
if [[ "$got" != "$COMMIT" ]]; then
  echo "ped-sim $TAG resolved to $got, expected $COMMIT; refusing to build" >&2
  exit 1
fi
P="$CONDA_PREFIX"
make -C "$DIR" -f Makefile-gsl --quiet \
  GPP="$P/bin/x86_64-conda-linux-gnu-g++" GCC="$P/bin/x86_64-conda-linux-gnu-gcc" \
  CFLAGS="-Wall -O2 -DUSEGSL -I$P/include" \
  LIBS="-L$P/lib -Wl,-rpath,$P/lib -lz -lgsl -lgslcblas"
