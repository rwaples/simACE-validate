"""Summary of the unequal-r two-trait pilot (rank-1 mate correlation; am/xft_rank1.py).

Reads the simACE and xftsim JSONL of `CELLS=rank1 bash am/pilot.sh` (committed as
am/pilot/rank1_simace.jsonl and rank1_xft.jsonl), keeps the rank-1 records (r is
the string "R1/R2"; every kept record must share r, g and h2) and prints, per statistic and generation, the simACE and xftsim means, the
Welch p between them, and each tool's z against a two-trait Bulmer recursion
with a general mate correlation matrix R (bulmer2 in am/stats.py takes one r).

Usage: python -m am.rank1_summary SIMACE_JSONL XFT_JSONL
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

COUNTS = ("n", "n_couples", "n_fs_pairs")  # family structure differs by design (differences.md D4)


def bulmer2_general(a0: float, e: float, rmat: np.ndarray, g: int) -> list[dict[str, float]]:
    """am/stats.py bulmer2 with the liability mate correlation R = rmat (female rows, male columns)."""
    sa = np.eye(2) * a0
    out = [{}]
    for _ in range(g):
        sy = sa + e * np.eye(2)
        d = np.sqrt(np.diag(sy))
        k = sa @ np.linalg.solve(sy, np.eye(2))
        c = k @ (np.outer(d, d) * rmat) @ k.T
        sda = np.sqrt(np.diag(sa))
        mate = {
            "r_mate_A": c[0, 0] / sa[0, 0],
            "r_mate_A_t2": c[1, 1] / sa[1, 1],
            "r_mate_A12": c[0, 1] / (sda[0] * sda[1]),
            "r_mate_A21": c[1, 0] / (sda[0] * sda[1]),
            "r_mate_L": rmat[0, 0],
            "r_mate_L_t2": rmat[1, 1],
            "r_mate_L12": rmat[0, 1],
            "r_mate_L21": rmat[1, 0],
            "r_w_L_parents": sy[0, 1] / (d[0] * d[1]),
        }
        sa = sa / 2 + (c + c.T) / 4 + np.eye(2) * a0 / 2
        sy = sa + e * np.eye(2)
        out.append(
            {
                "var_A": sa[0, 0],
                "var_A_t2": sa[1, 1],
                "var_L": sy[0, 0],
                "var_L_t2": sy[1, 1],
                "r_w_A": sa[0, 1] / np.sqrt(sa[0, 0] * sa[1, 1]),
                "r_w_L": sy[0, 1] / np.sqrt(sy[0, 0] * sy[1, 1]),
                **mate,
            }
        )
    return out


def load(path: Path) -> tuple[tuple, dict[tuple[str, int], np.ndarray]]:
    """The shared (r, g, h2) of the rank-1 records, and per (stat, g) their per-seed values."""
    runs = [json.loads(line) for line in path.read_text().splitlines() if line.strip()]
    runs = [run for run in runs if run["traits"] == 2 and "/" in str(run["r"])]
    spec = {(run["r"], run["g"], run["h2"]) for run in runs}
    if len(spec) != 1:
        raise SystemExit(f"{path}: need one rank-1 (r, g, h2), got {sorted(spec)}")
    vals: dict[tuple[str, int], list[float]] = {}
    for run in runs:
        for row in run["rows"]:
            for k, v in row.items():
                if k != "generation" and row["generation"] > 0:
                    vals.setdefault((k, row["generation"]), []).append(v)
    return spec.pop(), {k: np.array(v) for k, v in vals.items()}


def z(v: np.ndarray, exp: float) -> float:
    """Return the z of the mean against exp."""
    return (v.mean() - exp) / (v.std(ddof=1) / np.sqrt(len(v)))


spec, sim = load(Path(sys.argv[1]))
spec_xft, xft = load(Path(sys.argv[2]))
if spec != spec_xft:
    raise SystemExit(f"simACE {spec} and xftsim {spec_xft} differ")
(R1, R2), G, h2 = (float(v) for v in spec[0].split("/")), spec[1], spec[2]
x = np.sqrt(R1 * R2)
pred = bulmer2_general(h2, 1 - h2, np.array([[R1, x], [x, R2]]), G)
keys = sorted(set(sim) & set(xft), key=lambda k: (k[1], k[0]))
for tool, only in (("simACE", set(sim) - set(xft)), ("xftsim", set(xft) - set(sim))):
    if only:
        names = ", ".join(f"{k} g{g}" for k, g in sorted(only, key=lambda k: (k[1], k[0])))
        print(f"{tool} only, not compared: {names}")
n_sim, n_xft = len(sim[keys[0]]), len(xft[keys[0]])
print(
    f"rank-1 two-trait AM, G {G}, h2 {h2}, r1 {R1}, r2 {R2}, cross {x:.4f}; simACE n = {n_sim}, xftsim n = {n_xft} seeds"
)
print()
print("| g | stat | simACE | xftsim | Welch p | Bulmer | z simACE | z xftsim |")
print("|---|---|---|---|---|---|---|---|")
ps, zs = [], {"simACE": [], "xftsim": []}
for stat, g in keys:
    a, b = sim[(stat, g)], xft[(stat, g)]
    p = stats.ttest_ind(a, b, equal_var=False).pvalue
    if stat.removesuffix("_t2") not in COUNTS:
        ps.append((p, stat, g))
    exp = pred[g].get(stat)
    cells = ["", "", ""] if exp is None else [f"{exp:.4f}", f"{z(a, exp):+.2f}", f"{z(b, exp):+.2f}"]
    if exp is not None:
        zs["simACE"].append((abs(z(a, exp)), stat, g))
        zs["xftsim"].append((abs(z(b, exp)), stat, g))
    print(f"| {g} | {stat} | {a.mean():.4f} | {b.mean():.4f} | {p:.3g} | " + " | ".join(cells) + " |")
print()
p, stat, g = min(ps)
print(
    f"Welch rows other than {', '.join(COUNTS)}: {len(ps)}; smallest p {p:.3g} ({stat}, g {g}); Bonferroni 0.05/{len(ps)} = {0.05 / len(ps):.3g}"
)
for tool, v in zs.items():
    m, stat, g = max(v)
    print(f"{tool} vs Bulmer: {len(v)} rows; largest |z| {m:.2f} ({stat}, g {g})")
