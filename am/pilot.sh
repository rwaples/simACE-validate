#!/usr/bin/env bash
# Pilot runs for issue #4: simACE assortative mating vs xftsim LinearAssortativeMatingRegime.
# Pilot seeds (GLOSSARY) come from 9000001.., a range never scored.
# Usage: bash am/pilot.sh OUTDIR [JOBS]
# Writes OUTDIR/xft.jsonl and OUTDIR/simace.jsonl, one JSON line per run.
set -euo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
out=$1
jobs=${2:-4}
mkdir -p "$out/tmp/sim" "$out/tmp/xft"
xpy=$here/.pixi/envs/xftsim/bin/python
dpy=$here/.pixi/envs/default/bin/python
export PATH=$here/.pixi/envs/default/bin:$PATH
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
N=10000
G=6
R=40
tasks=$out/tmp/tasks.txt
: > "$tasks"
sim() { # seed r traits
  echo "$dpy $here/am/run_simace.py $1 $N $G $2 $out/tmp $3 > $out/tmp/sim/s$1_r$2_t$3.json" >> "$tasks"
}
xft() { # seed m r traits
  echo "$xpy $here/am/run_xft.py $1 $N $2 $G $3 $4 > $out/tmp/xft/s$1_m$2_r$3_t$4.json" >> "$tasks"
}
for i in $(seq 1 $R); do
  s=$((9000000 + i))
  # single-trait assortment
  for r in 0 0.5 0.8; do sim $s $r 1; xft $s 1000 $r 1; done
  for m in 250 4000; do for r in 0.5 0.8; do xft $s $m $r 1; done; done
  # two-trait assortment, every mate correlation r (feasible for r < 0.5 with independent traits)
  for r in 0.2 0.4; do sim $s $r 2; xft $s 1000 $r 2; done
done
xargs -P "$jobs" -I{} bash -c "{}" < "$tasks"
cat "$out"/tmp/sim/*.json > "$out/simace.jsonl"
cat "$out"/tmp/xft/*.json > "$out/xft.jsonl"
wc -l "$out"/*.jsonl
