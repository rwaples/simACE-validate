"""Summarise the C-transmission pilot (rwaples/simACE-validate#5); prints Markdown.

Usage: python -m am.ctm.summarize PILOT_DIR   (reads simace.jsonl and xft.jsonl from `CELLS=ctm bash am/pilot.sh`)

For each cell (r 0: random mating with half-sibs; r 0.5: LAM) and each statistic and generation:
simACE and xftsim mean (SE), the random-mating law (`am.bulmer.predict_ctm`, r 0 only), the
Welch p, and the minimum mean gap a Welch t at R = 100 per side detects with 80% power at
alpha/K = 0.05/1461 (the K this case brings the harness to). Then check 3: Var(C) by
generation on each side against C, under AM.
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats as st

from am.bulmer import predict_ctm

BM, BF, C = 0.5, 0.3, 0.2
ALPHA_K = 0.05 / 1461
R_TARGET = 100
STATS = [
    "var_A",
    "var_L",
    "var_C",
    "r_FS_L",
    "r_FS_A",
    "r_PO_L",
    "r_PO_A",
    "r_MO_L",
    "r_FO_L",
    "r_MO_C",
    "r_FO_C",
    "r_MHS_C",
    "r_PHS_C",
    "r_MHS_L",
    "r_PHS_L",
    "b_MP_L",
    "r_mate_L",
    "r_mate_A",
]


def load(path: Path) -> dict[float, list[dict]]:
    """Records by r, each as {stat_g: value} (generation 0 dropped)."""
    out = defaultdict(list)
    for line in path.read_text().splitlines():
        if line.strip():
            rec = json.loads(line)
            flat = {
                f"{k}_g{row['generation']}": v for row in rec["rows"] if row["generation"] > 0 for k, v in row.items()
            }
            out[float(rec["r"])].append(flat)
    return out


def mdd(sa: float, sb: float) -> float:
    """Normal-approximation minimum detectable mean gap at R_TARGET per side, 80% power, alpha/K."""
    z = st.norm.isf(ALPHA_K / 2) + st.norm.isf(0.2)
    return z * np.sqrt((sa**2 + sb**2) / R_TARGET)


def main() -> None:
    """Print the tables."""
    root = Path(sys.argv[1])
    sim, xft = load(root / "simace.jsonl"), load(root / "xft.jsonl")
    law = predict_ctm(0.5, C, 0.5 - C, BM, BF)
    for r in sorted(sim):
        a_runs, b_runs = sim[r], xft[r]
        print(f"## r = {r} (simACE {len(a_runs)} seeds, xftsim {len(b_runs)} seeds)\n")
        print("| stat | g | simACE | xftsim | law | Welch p | MDD (R 100) |")
        print("|---|---|---|---|---|---|---|")
        ps = []
        for s in STATS:
            for g in range(1, 7):
                key = f"{s}_g{g}"
                a = np.array([x[key] for x in a_runs if key in x], float)
                b = np.array([x[key] for x in b_runs if key in x], float)
                if len(a) < 3 or len(b) < 3 or np.isnan(a).any() or np.isnan(b).any():
                    continue
                p = st.ttest_ind(a, b, equal_var=False).pvalue
                ps.append((p, key))
                pred = f"{law[s]:.4f}" if r == 0 and s in law else ""
                print(
                    f"| {s} | {g} | {a.mean():.4f} ({a.std(ddof=1) / np.sqrt(len(a)):.4f}) "
                    f"| {b.mean():.4f} ({b.std(ddof=1) / np.sqrt(len(b)):.4f}) | {pred} | {p:.3g} "
                    f"| {mdd(a.std(ddof=1), b.std(ddof=1)):.4f} |"
                )
        ps.sort()
        print(f"\n{len(ps)} rows; smallest Welch p: {', '.join(f'{k} {p:.3g}' for p, k in ps[:3])}")
        print(f"rows with p < 0.05: {sum(p < 0.05 for p, _ in ps)} (expected {0.05 * len(ps):.1f})\n")
        if r == 0:
            continue
        print("Check 3: Var(C) minus C by generation, mean (95% CI)\n")
        print("| g | simACE (realised top-up) | xftsim (fixed residual) |")
        print("|---|---|---|")
        for g in range(1, 7):
            cells = []
            for runs in (a_runs, b_runs):
                d = np.array([x[f"var_C_g{g}"] for x in runs]) - C
                h = st.t.ppf(0.975, len(d) - 1) * d.std(ddof=1) / np.sqrt(len(d))
                cells.append(f"{d.mean():+.4f} ({d.mean() - h:+.4f}, {d.mean() + h:+.4f})")
            print(f"| {g} | {cells[0]} | {cells[1]} |")
        print()


if __name__ == "__main__":
    main()
