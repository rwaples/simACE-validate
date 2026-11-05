"""Compare a simACE-ts Case's Gaussian-sensitive INFO statistics across architectures (differences.md D3).

For each INFO statistic that reads a Gaussian liability (``tetrachoric_*``,
``affected_corr_*``, ``h2_observed_*``) and ``prevalence``, and each
architecture other than the reference, prints the paired per-seed
difference from the reference: its mean, SE, t and n. Then it prints how
many mean differences are negative, the largest |t|, and the two-sided
Bonferroni cutoff over the computable comparisons. The comparisons share
the reference and the seeds' phenotypes, so they are correlated; the
cutoff is conservative, and a sign count is not a test.

With ``--mde-k K`` it also prints each scored row's minimum detectable mean
gap at alpha / K, by the formula of run_harness.minimum_detectable, for the
K of a full run.

Usage::

    pixi run python simace_ts/sweep_info.py [--case simace_ts_p2] [--reference c1m_a05] [--mde-k 846]
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "simace_ts"))
from score import case_rows, load_parts  # noqa: E402

GAUSSIAN_SENSITIVE = ("tetrachoric_", "affected_corr_", "h2_observed_", "prevalence")
ALPHA = 0.05


def contrasts(parts: list[dict], reference: str) -> list[tuple[str, str, float, float, float, int]]:
    """Return (stat, architecture, mean, SE, t, n) of each architecture's paired difference from ``reference``."""
    archs = [a for a in parts[0]["architectures"] if a != reference]
    names = [k for k in parts[0]["architectures"][reference]["info"] if k.startswith(GAUSSIAN_SENSITIVE)]
    out = []
    for stat in names:
        for arch in archs:
            d = np.array(
                [
                    p["architectures"][arch]["info"][stat] - p["architectures"][reference]["info"][stat]
                    for p in parts
                    if p["architectures"][arch]["info"].get(stat) is not None
                    and p["architectures"][reference]["info"].get(stat) is not None
                ]
            )
            if len(d) < 2:
                continue
            se = d.std(ddof=1) / np.sqrt(len(d))
            out.append((stat, arch, float(d.mean()), float(se), float(d.mean() / se), len(d)))
    return out


def main() -> int:
    """Print the contrast table and its summary for one case."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--case", default="simace_ts_p2")
    ap.add_argument("--reference", default="c1m_a05")
    ap.add_argument("--mde-k", type=int, help="also print scored rows' minimum detectable gap at alpha / K")
    args = ap.parse_args()
    case = {c["id"]: c for c in yaml.safe_load((ROOT / "cases.yaml").read_text())["cases"]}[args.case]
    parts = load_parts(case)

    rows = contrasts(parts, args.reference)
    print(f"| Statistic | Architecture | Mean difference from {args.reference} | SE | t | n |")
    print("|---|---|---|---|---|---|")
    for stat, arch, mean, se, t, n in rows:
        print(f"| {stat} | {arch} | {mean:+.5f} | {se:.5f} | {t:+.2f} | {n} |")
    worst = max(rows, key=lambda r: abs(r[4]))
    largest = max(rows, key=lambda r: abs(r[2]))
    n = worst[5]
    cutoff = stats.t.ppf(1 - ALPHA / len(rows) / 2, n - 1)
    print()
    print(f"comparisons: {len(rows)}; negative mean differences: {sum(r[2] < 0 for r in rows)}")
    for arch in dict.fromkeys(r[1] for r in rows):
        mine = [r for r in rows if r[1] == arch]
        print(f"  {arch}: {sum(r[2] < 0 for r in mine)} of {len(mine)} negative")
    print(f"largest |t|: {abs(worst[4]):.2f} ({worst[0]}, {worst[1]}, mean {worst[2]:+.5f})")
    print(f"largest |mean difference|: {abs(largest[2]):.5f} ({largest[0]}, {largest[1]})")
    print(f"two-sided Bonferroni cutoff at {ALPHA} over {len(rows)}, df {n - 1}: {cutoff:.2f}")

    if args.mde_k:
        threshold = ALPHA / args.mde_k
        print()
        print(f"minimum detectable gap at alpha/K = {threshold:.3g} (K = {args.mde_k}):")
        for stat, gaps in case_rows(case, parts).scored.items():
            m = len(gaps)
            factor = (stats.t.ppf(1 - threshold / 2, m - 1) + stats.t.ppf(0.8, m - 1)) / np.sqrt(m)
            print(f"  {stat}: {factor * gaps.std(ddof=1):.4g} (n {m})")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
