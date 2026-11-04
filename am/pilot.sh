#!/usr/bin/env bash
# Pilot runs for issue #4: simACE assort1 vs xftsim LinearAssortativeMatingRegime.
# Pilot seeds (GLOSSARY) come from 9000001.., a range never scored.
# Usage: bash am/pilot.sh OUTDIR [JOBS]
# Writes OUTDIR/xft.jsonl and OUTDIR/simace.jsonl, one JSON line per run.
set -euo pipefail
here=$(cd "$(dirname "$0")/.." && pwd)
out=$1
jobs=${2:-4}
mkdir -p "$out/tmp"
xpy=$here/.pixi/envs/xftsim/bin/python
dpy=$here/.pixi/envs/default/bin/python
export PATH=$here/.pixi/envs/default/bin:$PATH
export OMP_NUM_THREADS=1 NUMBA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
N=10000
G=6
R=40
tasks=$out/tmp/tasks.txt
: > "$tasks"
for i in $(seq 1 $R); do
  s=$((9000000 + i))
  for r in 0 0.5 0.8; do
    echo "$dpy $here/am/run_simace.py $s $N $G $r $out/tmp >> $out/simace.jsonl" >> "$tasks"
    echo "$xpy $here/am/run_xft.py $s $N 1000 $G $r >> $out/xft.jsonl" >> "$tasks"
  done
  for m in 250 4000; do for ra in 0.5 0.8; do
    echo "$xpy $here/am/run_xft.py $s $N $m $G $ra >> $out/xft.jsonl" >> "$tasks"; done
  done
done
# each line appends one short JSON line (< PIPE_BUF), so concurrent appends do not interleave
xargs -P "$jobs" -I{} bash -c "{}" < "$tasks"
wc -l "$out"/*.jsonl
