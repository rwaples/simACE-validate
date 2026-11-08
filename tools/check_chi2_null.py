"""Check the exact chi-square reference of calibration.py (G4) and print every number status.md and ADR 0001 quote.

1. calibration.chi2_null_tail against brute-force enumeration in integer arithmetic.
2. The size (P(p < alpha), the harness's rule) of the exact and the asymptotic reference, from the exact pmf.
3. The simulated rejection rate of both references over --draws null count vectors per R.
4. calibrate() on uniform PITs at conventional alphas.
5. The chi-square p-values of results/results.tsv against those of --base (the asymptotic run).

Run: pixi run python -P tools/check_chi2_null.py [--draws 100000000] [--base ec94688]
"""

from __future__ import annotations

import argparse
import csv
import io
import itertools
import math
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from calibration import CHI2_BINS, calibrate, chi2_null_tail  # noqa: E402  (after the sys.path entry)

ALPHAS = {"0.05/587": 0.05 / 587, "0.05/1364": 0.05 / 1364}
R_VALUES = (30, 50, 100)  # every R of a calibrated case in cases.yaml


def enumerated_tail(R: int, bins: int) -> dict[int, float]:
    """Return P(S >= s) for each reachable s, summed over every count vector in integer arithmetic."""
    ways: defaultdict[int, int] = defaultdict(int)
    for cuts in itertools.combinations(range(R + bins - 1), bins - 1):  # stars and bars
        edges = (-1, *cuts, R + bins - 1)
        c = [edges[i + 1] - edges[i] - 1 for i in range(bins)]
        ways[sum(x * x for x in c)] += math.factorial(R) // math.prod(math.factorial(x) for x in c)
    total = bins**R
    tail, acc = {}, 0
    for s in sorted(ways, reverse=True):
        acc += ways[s]
        tail[s] = acc / total
    return tail


def check_enumeration() -> None:
    """Compare the DP with enumeration, including the production 10 bins."""
    for R, bins in ((8, 4), (12, 3), (6, 10), (9, 10)):
        exact = enumerated_tail(R, bins)
        dp = chi2_null_tail(R, bins)
        worst = max(abs(dp[s] - p) / p for s, p in exact.items())
        print(f"enumeration R={R} bins={bins}: {len(exact)} values of S, max relative |DP - enumeration| = {worst:.3g}")


def sizes(R: int, alpha: float) -> tuple[float, float]:
    """Return P(p < alpha) under the null for the exact and the asymptotic p, from the exact pmf."""
    tail = chi2_null_tail(R, CHI2_BINS)
    pmf = tail - np.append(tail[1:], 0.0)
    s = np.arange(len(tail))
    asym = stats.chi2.sf(CHI2_BINS / R * s - R, CHI2_BINS - 1)
    return float(pmf[tail < alpha].sum()), float(pmf[asym < alpha].sum())


def check_sizes() -> None:
    """Print the exact size of both references at each alpha and R."""
    for name, alpha in ALPHAS.items():
        for R in R_VALUES:
            exact, asym = sizes(R, alpha)
            print(f"size at alpha {name} = {alpha:.4g}, R={R}: exact {exact:.4g}, asymptotic {asym:.4g}")


def simulate(R: int, draws: int, seed: int) -> None:
    """Print the simulated rejection rate of both references at alpha 0.05/1364."""
    alpha = ALPHAS["0.05/1364"]
    rng = np.random.default_rng(seed)
    tail = chi2_null_tail(R, CHI2_BINS)
    hits = {"exact": 0, "asymptotic": 0}
    chunk = 5_000_000
    for start in range(0, draws, chunk):
        counts = rng.multinomial(R, np.full(CHI2_BINS, 1 / CHI2_BINS), size=min(chunk, draws - start))
        s = np.einsum("ij,ij->i", counts, counts)
        hits["exact"] += int(np.sum(tail[s] < alpha))
        hits["asymptotic"] += int(np.sum(stats.chi2.sf(CHI2_BINS / R * s - R, CHI2_BINS - 1) < alpha))
    for name, h in hits.items():
        ci = stats.binomtest(h, draws).proportion_ci(0.95, method="exact")
        print(f"simulated R={R} {name}: {h}/{draws} = {h / draws:.4g} (95% CI {ci.low:.4g} to {ci.high:.4g})")


def end_to_end(R: int, samples: int, seed: int) -> None:
    """Print calibrate()'s chi-square rejection rate on uniform PITs."""
    rng = np.random.default_rng(seed)
    p = np.array([calibrate(rng.uniform(size=R), np.zeros(R)).p_values["chi2"] for _ in range(samples)])
    print(
        f"calibrate R={R}, {samples} uniform samples: rejection at 0.05 {np.mean(p < 0.05):.4g}, at 0.01 {np.mean(p < 0.01):.4g}"
    )


def chi2_p(path_text: str) -> dict[tuple[str, str], tuple[str, float]]:
    """Return (status, chi2 p) of each row of a results.tsv text that has a chi-square."""
    out = {}
    for r in csv.DictReader(io.StringIO(path_text), delimiter="\t"):
        p = dict(part.split(" ", 1) for part in r["p_values"].split(", ") if " " in part)
        if "chi2" in p:
            out[(r["case"], r["stat"])] = (r["status"], float(p["chi2"]))
    return out


def compare(base: str) -> None:
    """Print how the chi-square p-values of results/results.tsv moved from those of base."""
    old_text = subprocess.run(
        ["git", "-C", str(ROOT), "show", f"{base}:results/results.tsv"], capture_output=True, text=True, check=True
    ).stdout
    old, new = chi2_p(old_text), chi2_p((ROOT / "results" / "results.tsv").read_text())
    keys = sorted(old.keys() & new.keys())
    ratio = np.array([new[k][1] / old[k][1] for k in keys])
    moved = [k for k in keys if old[k][0] != new[k][0]]
    print(
        f"results.tsv against {base}: {len(keys)} chi-square rows ({len(old.keys() ^ new.keys())} unmatched), "
        f"{len(moved)} status changes; new/old p min {ratio.min():.4g}, median {np.median(ratio):.4g}, max {ratio.max():.4g}"
    )
    scored = [k for k in keys if new[k][0] != "INFO"]
    k = min(scored, key=lambda k: new[k][1])
    print(
        f"smallest scored chi-square p: {k[0]} {k[1]} {new[k][1]:.3g} (was {old[k][1]:.3g}), of {len(scored)} scored rows"
    )


def main() -> None:
    """Run every check."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--draws", type=int, default=100_000_000, help="null count vectors per R in the simulation")
    parser.add_argument(
        "--base", default="ec94688", help="git revision whose results.tsv used the asymptotic reference"
    )
    args = parser.parse_args()
    check_enumeration()
    check_sizes()
    for R in (50, 100):
        simulate(R, args.draws, seed=20261108 + R)
    end_to_end(50, 20_000, seed=1)
    compare(args.base)


if __name__ == "__main__":
    main()
