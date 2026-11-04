"""One simACE run for the assortative-mating pilot (issue #4); prints one JSON line.

Runs in the default pixi env against the pinned simace. Both traits have
A = 0.5, C = 0, E = 0.5 and are independent (rA = rC = rE = 0). No MZ twins,
mating_lambda 0.5, G_sim = G_ped = g.

- TRAITS = 1: mate correlation assort1 = r on trait 1 (single-trait
  assortment); trait 2 is a passenger (assort2 = 0).
- TRAITS = 2: both-trait assortment with every mate correlation r:
  assort1 = assort2 = r and assort_matrix [[r, r], [r, r]], the law of
  xftsim mating on the sum of two standardised traits (run_xft.py).

Usage: run_simace.py SEED N G R OUTDIR [TRAITS]
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
from stats import cross_trait_stats, generation_stats

seed, n, g = (int(a) for a in sys.argv[1:4])
r = float(sys.argv[4])
traits = int(sys.argv[6]) if len(sys.argv) > 6 else 1
out = Path(sys.argv[5]) / f"simace_s{seed}_n{n}_r{r}_t{traits}.parquet"
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
        "0",
        "--E1",
        "0.5",
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
ped = pl.read_parquet(out).sort("id")
version = subprocess.run(["simace", "--version"], capture_output=True, text=True, check=True).stdout.strip()


def values(frame: pl.DataFrame) -> list[tuple[np.ndarray, np.ndarray]]:
    """(A, liability) of each scored trait."""
    return [(frame[f"A{t}"].to_numpy(), frame[f"liability{t}"].to_numpy()) for t in range(1, traits + 1)]


def suffixed(d: dict, t: int) -> dict:
    """Trait 1 keys unchanged; trait 2 keys get `_t2`, as in run_xft.py."""
    return d if t == 0 else {f"{k}_t{t + 1}": v for k, v in d.items()}


rows = []
prev = None
for k in range(g):
    cur = ped.filter(pl.col("generation") == k)
    v = values(cur)
    row = {"generation": k + 1}
    if prev is None:
        for t, (a, y) in enumerate(v):
            row.update(suffixed({"var_A": float(np.var(a, ddof=1)), "var_L": float(np.var(y, ddof=1))}, t))
    else:
        pid = prev["id"].to_numpy()
        pos = np.searchsorted(pid, cur["mother"].to_numpy())
        fos = np.searchsorted(pid, cur["father"].to_numpy())
        assert (pid[pos] == cur["mother"].to_numpy()).all()
        assert (pid[fos] == cur["father"].to_numpy()).all()
        pv = values(prev)
        for t in range(traits):
            row.update(suffixed(generation_stats(v[t][0], v[t][1], pos, fos, pv[t][0], pv[t][1]), t))
        if traits == 2:
            row.update(cross_trait_stats(v, pv, pos, fos))
    rows.append(row)
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
            "traits": traits,
            "h2": 0.5,
            "wall_s": round(time.perf_counter() - t0, 2),
            "rows": rows,
        }
    )
)
