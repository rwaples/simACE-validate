"""One simACE run for the assortative-mating pilot (issue #4); prints one JSON line.

Runs in the default pixi env against the pinned simace. Trait 1: A1 = 0.5,
C1 = 0, E1 = 0.5, mate correlation assort1 = r (single-trait assortment);
trait 2 is an independent passenger (A2 = E2 = 0.5, assort2 = 0). No MZ
twins, mating_lambda 0.5, G_sim = G_ped = g.

Usage: run_simace.py SEED N G R OUTDIR
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

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stats import generation_stats

seed, n, g = (int(a) for a in sys.argv[1:4])
r = float(sys.argv[4])
out = Path(sys.argv[5]) / f"simace_s{seed}_n{n}_r{r}.parquet"
out.parent.mkdir(parents=True, exist_ok=True)

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
        "0",
        "--E1",
        "0.5",
        "--A2",
        "0.5",
        "--C2",
        "0",
        "--E2",
        "0.5",
        "--assort1",
        str(r),
        "--output-pedigree",
        str(out),
    ],
    check=True,
)
ped = pl.read_parquet(out).sort("id")
version = subprocess.run(["simace", "--version"], capture_output=True, text=True, check=True).stdout.strip()

rows = []
prev = None
for k in range(g):
    cur = ped.filter(pl.col("generation") == k)
    a = cur["A1"].to_numpy()
    y = cur["liability1"].to_numpy()
    if prev is None:
        rows.append({"generation": 1, "var_A": float(np.var(a, ddof=1)), "var_L": float(np.var(y, ddof=1))})
    else:
        pid = prev["id"].to_numpy()
        pos = np.searchsorted(pid, cur["mother"].to_numpy())
        fos = np.searchsorted(pid, cur["father"].to_numpy())
        assert (pid[pos] == cur["mother"].to_numpy()).all()
        assert (pid[fos] == cur["father"].to_numpy()).all()
        stats = generation_stats(a, y, pos, fos, prev["A1"].to_numpy(), prev["liability1"].to_numpy())
        rows.append({"generation": k + 1, **stats})
    prev = cur
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
            "h2": 0.5,
            "wall_s": round(time.perf_counter() - t0, 2),
            "rows": rows,
        }
    )
)
