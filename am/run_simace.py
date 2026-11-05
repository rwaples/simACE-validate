"""One simACE run for the assortative-mating pilot (issue #4); prints one JSON line.

Runs in the default pixi env against the pinned simace. Both traits have
A = 0.5, C = 0, E = 0.5 and are independent (rA = rC = rE = 0). No MZ twins,
mating_lambda 0.5, G_sim = G_ped = g.

- TRAITS = 1: mate correlation assort1 = r on trait 1 (single-trait
  assortment); trait 2 is a passenger (assort2 = 0).
- TRAITS = 2: both-trait assortment with every mate correlation r:
  assort1 = assort2 = r and assort_matrix [[r, r], [r, r]], the law of
  xftsim mating on the sum of two standardised traits (run_xft.py).

C > 0 (one trait only) sets trait 1's household C to C and its E to 0.5 - C,
matching run_xft.py's household C.

Usage: python -m am.run_simace SEED N G R OUTDIR [TRAITS [C]]
The output's generation g is the offspring of g rounds of mating, as in
run_xft.py: simACE's recorded generation k follows k + 1 matings
(ADR 0001 §3), so it is reported as g = k + 1. Founders (g = 0) are not
recorded by simACE.
"""

import json
import subprocess
import sys
import time
from pathlib import Path

import polars as pl

from am.simace_stats import pedigree_rows

seed, n, g = (int(a) for a in sys.argv[1:4])
r = float(sys.argv[4])
traits = int(sys.argv[6]) if len(sys.argv) > 6 else 1
c = float(sys.argv[7]) if len(sys.argv) > 7 else 0.0
out = Path(sys.argv[5]) / f"simace_s{seed}_n{n}_r{r}_t{traits}_c{c}.parquet"
out.parent.mkdir(parents=True, exist_ok=True)

am = ["--assort1", str(r)]
if traits == 2:
    am = ["--assort1", str(r), "--assort2", str(r), "--assort-matrix", json.dumps([[r, r], [r, r]])]

t0 = time.perf_counter()
subprocess.run(
    [
        "simace",
        "simulate",
        "-q",
        "--seed",
        str(seed),
        "--N",
        str(n),
        "--G-ped",
        str(g),
        "--G-sim",
        str(g),
        "--mating-lambda",
        "0.5",
        "--p-mztwin",
        "0",
        "--A1",
        "0.5",
        "--C1",
        str(c),
        "--E1",
        str(0.5 - c),
        "--A2",
        "0.5",
        "--C2",
        "0",
        "--E2",
        "0.5",
        # simace's CLI defaults are rA 0.5 and rC 0.3; the matched configuration has independent traits
        "--rA",
        "0",
        "--rC",
        "0",
        "--rE",
        "0",
        *am,
        "--output-pedigree",
        str(out),
    ],
    check=True,
)
ped = pl.read_parquet(out)
version = subprocess.run(["simace", "--version"], capture_output=True, text=True, check=True).stdout.strip()
rows = pedigree_rows(ped, traits)
out.unlink()

print(
    json.dumps(
        {
            "tool": "simace",
            "version": version,
            "seed": seed,
            "n": n,
            "g": g,
            "r": r,
            "traits": traits,
            "C": c,
            "h2": 0.5,
            "wall_s": round(time.perf_counter() - t0, 2),
            "rows": rows,
        }
    )
)
