"""Probe for rwaples/simACE-validate#5, open question 2: simACE's household C when a mother has several partners.

Runs `simace simulate` with C transmission under random mating and reports, for the last
recorded generation, the C correlation and the share of a household id among maternal
half-sibs (MHS) and paternal half-sibs (PHS), against the random-mating law
MHS 1 - b_f^2, PHS b_f^2 (simACE `_SHARED_C`). Sharing the household id across a mother's
partners is what the xftsim construction (the mother's private draw) reproduces.

Usage (default env, repo root): python -m am.ctm.probe_halfsib SEED N G BM BF C OUTDIR
"""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl

seed, n, g = (int(a) for a in sys.argv[1:4])
bm, bf, c = (float(a) for a in sys.argv[4:7])
out = Path(sys.argv[7]) / f"halfsib_s{seed}.parquet"
out.parent.mkdir(parents=True, exist_ok=True)
args = {
    "seed": seed, "N": n, "G-ped": g, "G-sim": g, "mating-lambda": 0.5, "p-mztwin": 0,
    "A1": 0.5, "C1": c, "E1": 0.5 - c, "A2": 0.5, "C2": 0, "E2": 0.5, "rA": 0, "rC": 0, "rE": 0,
    "Ctm-maternal1": bm, "Ctm-paternal1": bf,
}  # fmt: skip
cmd = [
    "simace",
    "simulate",
    "-q",
    *[x for k, v in args.items() for x in (f"--{k}", str(v))],
    "--output-pedigree",
    str(out),
]
subprocess.run(cmd, check=True)
full = pl.read_parquet(out)
out.unlink()
ped = full.filter(pl.col("generation") == pl.col("generation").max())
parent_c = dict(zip(full["id"].to_list(), full["C1"].to_list(), strict=True))
mother, father = ped["mother"].to_numpy(), ped["father"].to_numpy()
c1, hh = ped["C1"].to_numpy(), ped["household_id"].to_numpy()

rows = []
for kind, same, other in (("MHS", mother, father), ("PHS", father, mother), ("FS", mother, father)):
    order = np.argsort(same, kind="stable")
    pi, pj = [], []
    s = same[order]
    starts = np.flatnonzero(np.r_[True, s[1:] != s[:-1]])
    for a, b in zip(starts, np.r_[starts[1:], len(s)], strict=True):
        idx = order[a:b]
        ii, jj = np.triu_indices(len(idx), 1)
        pi.append(idx[ii])
        pj.append(idx[jj])
    pi, pj = np.concatenate(pi), np.concatenate(pj)
    keep = (other[pi] == other[pj]) if kind == "FS" else (other[pi] != other[pj])
    pi, pj = pi[keep], pj[keep]
    x, y = np.r_[c1[pi], c1[pj]], np.r_[c1[pj], c1[pi]]
    rows.append({
        "pair": kind,
        "n_pairs": len(pi),
        "r_C": round(float(np.corrcoef(x, y)[0, 1]), 4),
        "share_household": round(float(np.mean(hh[pi] == hh[pj])), 4),
    })  # fmt: skip
for kind, par in (("MO", mother), ("FO", father)):
    pc = np.array([parent_c[k] for k in par])
    rows.append({"pair": kind, "n_pairs": len(par), "r_C": round(float(np.corrcoef(pc, c1)[0, 1]), 4)})
expected = {"MHS": 1 - bf**2, "PHS": bf**2, "FS": 1.0, "MO": bm, "FO": bf}
print(json.dumps({"seed": seed, "n": n, "bm": bm, "bf": bf, "C": c, "var_C": round(float(np.var(c1, ddof=1)), 4)}))
for r in rows:
    print(json.dumps({**r, "expected_r_C": expected[r["pair"]]}))
