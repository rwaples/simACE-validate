"""Follow-up of a pilot signal: simACE parent-offspring correlations with household C, random mating.

The household-C pilot (am/pilot/summary.md) gave simACE r_FO_L at C = 0.2,
r = 0, g = 5 of 0.2440 (z -3.67 against 0.25). This reruns simACE alone at the
same configuration on fresh pilot seeds and reports, by g, the full-sib,
mother-offspring and father-offspring liability correlations with their z
against predict1 (C adds to full sibs only; household = mother).

Usage: python -m am.followup_c FIRST_SEED N_SEEDS WORKDIR
"""

import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl

from am.simace_stats import pedigree_rows
from am.stats import predict1

first, n_seeds = int(sys.argv[1]), int(sys.argv[2])
work = Path(sys.argv[3])
work.mkdir(parents=True, exist_ok=True)
N, G, C = 10000, 6, 0.2
STATS = ("r_FS_L", "r_MO_L", "r_FO_L", "r_PO_L")


def one(seed: int) -> list[dict]:
    """Per-generation statistics of one simACE pedigree."""
    p = work / f"c_{seed}.parquet"
    cmd = ["simace", "simulate", "-q", "--seed", str(seed), "--N", str(N), "--G-ped", str(G), "--G-sim", str(G)]
    cmd += ["--mating-lambda", "0.5", "--p-mztwin", "0", "--A1", "0.5", "--C1", str(C), "--E1", str(0.5 - C)]
    cmd += ["--A2", "0.5", "--C2", "0", "--E2", "0.5", "--rA", "0", "--rC", "0", "--rE", "0"]
    cmd += ["--assort1", "0", "--output-pedigree", str(p)]
    subprocess.run(cmd, check=True, env=os.environ | {"OMP_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1"})
    rows = pedigree_rows(pl.read_parquet(p), 1)
    p.unlink()
    return rows


with ThreadPoolExecutor(6) as pool:
    runs = list(pool.map(one, range(first, first + n_seeds)))

print(f"simACE, C = {C}, r = 0, N = {N}, seeds {first}..{first + n_seeds - 1}")
print("| g | stat | mean | SE | expected | z |")
print("|---|---|---|---|---|---|")
for g in range(2, G + 1):
    for s in STATS:
        v = np.array([next(r for r in run if r["generation"] == g)[s] for run in runs])
        exp = predict1(0.5, 0.5, 0.0, g, c=C)[s]
        se = v.std(ddof=1) / np.sqrt(len(v))
        print(f"| {g} | {s} | {v.mean():.4f} | {se:.4f} | {exp:.4f} | {(v.mean() - exp) / se:+.2f} |")
