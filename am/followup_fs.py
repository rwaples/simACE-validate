"""Follow-up of a pilot signal: simACE full-sib liability correlation under random mating.

The 40-seed pilot (am/pilot/summary.md) gave simACE r_FS_L at r = 0, g = 3
of 0.2418 (z -4.04 against 0.25). This reruns simACE alone at the same
configuration on fresh pilot seeds and reports, by g, the full-sib
correlation of L, A and E, with z against 0.25, 0.5 and 0.

Usage: followup_fs.py FIRST_SEED N_SEEDS WORKDIR
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stats import _sib_pairs, _sym_corr

first, n_seeds = int(sys.argv[1]), int(sys.argv[2])
work = Path(sys.argv[3])
work.mkdir(parents=True, exist_ok=True)
N, G = 10000, 6
out = {}
for seed in range(first, first + n_seeds):
    p = work / f"fs_{seed}.parquet"
    cmd = ["simace", "simulate", "-q", "--seed", str(seed), "--N", str(N), "--G-ped", str(G), "--G-sim", str(G)]
    cmd += ["--mating-lambda", "0.5", "--p-mztwin", "0", "--A1", "0.5", "--C1", "0", "--E1", "0.5"]
    cmd += ["--A2", "0.5", "--C2", "0", "--E2", "0.5", "--rA", "0", "--rC", "0", "--rE", "0"]
    cmd += ["--assort1", "0", "--output-pedigree", str(p)]
    subprocess.run(cmd, check=True)
    ped = pl.read_parquet(p)
    p.unlink()
    for k in range(1, G):
        cur = ped.filter(pl.col("generation") == k)
        key = cur["mother"].to_numpy().astype(np.int64) * (2**31) + cur["father"].to_numpy()
        _, couple = np.unique(key, return_inverse=True)
        i, j = _sib_pairs(couple)
        for col in ["liability1", "A1", "E1"]:
            out.setdefault((k + 1, col), []).append(_sym_corr(cur[col].to_numpy().astype(np.float64), i, j))

print(f"simACE, r = 0, N = {N}, seeds {first}..{first + n_seeds - 1}")
print("| g | stat | mean | SE | expected | z |")
print("|---|---|---|---|---|---|")
for (g, col), vals in sorted(out.items()):
    v = np.array(vals)
    exp = {"liability1": 0.25, "A1": 0.5, "E1": 0.0}[col]
    se = v.std(ddof=1) / np.sqrt(len(v))
    print(f"| {g} | r_FS_{col} | {v.mean():.4f} | {se:.4f} | {exp} | {(v.mean() - exp) / se:+.2f} |")
