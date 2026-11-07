"""Compare a simACE-ts Case's Gaussian-sensitive INFO statistics across architectures (differences.md D3).

For each INFO statistic that reads a Gaussian liability (``tetrachoric_*``,
``affected_corr_*``, ``h2_observed_*``) and ``prevalence``, and each
architecture other than the reference, prints the paired per-seed
difference from the reference: its mean, SE, t and n. Then it prints how
many mean differences are negative, the largest |t|, and the two-sided
Bonferroni cutoff over the computable comparisons. The comparisons share
the reference and the seeds' phenotypes, so they are correlated; the
cutoff is conservative, and a sign count is not a test.

Then the joint test (``joint_test``), once with each architecture as the
reference. Per seed and statistic, the contrast is the reference's value
minus the mean of the other architectures; the statistic is the sum over
statistics of each contrast's squared one-sample t across seeds. Its null is
that the architectures are exchangeable within a seed, so the permutation
draws, independently per seed, which architecture plays the reference. That
keeps each seed's dependence across statistics and architectures, which the
pairwise comparisons above ignore. Statistics empty in any seed or
architecture are left out. Two limits: a rejection means the architectures
are not exchangeable, and one outlier also shifts every other reference's
contrast (by -1/(m-1) of its gap), so read the line for every reference; and
averaging over the others dilutes a gap between the reference and a single
architecture, which only the pairwise table shows.

With ``--mde-k K`` it also prints each scored row's minimum detectable mean
gap at alpha / K, by the formula of run_harness.minimum_detectable, for the
K of a full run.

Usage::

    pixi run python simace_ts/sweep_info.py [--case simace_ts_p2] [--reference c1m_a05] [--mde-k 1364]

Parts built under another simace (e.g. 1f8bd5d's, at 0.2.0) are read with
``--parts-dir DIR --accept-part-simace``; every other staleness check still applies.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "simace_ts"))
from score import PARTS, case_rows, load_parts  # noqa: E402

GAUSSIAN_SENSITIVE = ("tetrachoric_", "affected_corr_", "h2_observed_", "prevalence")
ALPHA = 0.05


def gaussian_names(parts: list[dict], reference: str) -> list[str]:
    """Return the reference's Gaussian-sensitive INFO statistics, in Part order."""
    return [k for k in parts[0]["architectures"][reference]["info"] if k.startswith(GAUSSIAN_SENSITIVE)]


def one_sample_t(d: np.ndarray) -> np.ndarray:
    """Return the one-sample t of ``d`` along its second-to-last axis (seeds)."""
    n = d.shape[-2]
    with np.errstate(divide="ignore", invalid="ignore"):  # a zero SD gives inf/NaN, which the caller rejects
        return d.mean(axis=-2) / (d.std(axis=-2, ddof=1) / np.sqrt(n))


def contrasts(parts: list[dict], reference: str) -> list[tuple[str, str, float, float, float, int]]:
    """Return (stat, architecture, mean, SE, t, n) of each architecture's paired difference from ``reference``."""
    archs = [a for a in parts[0]["architectures"] if a != reference]
    out = []
    for stat in gaussian_names(parts, reference):
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


def joint_test(parts: list[dict], reference: str, draws: int, seed: int = 0) -> tuple[float, float, list[str]]:
    """Return (statistic, permutation p, statistics used) of the joint test with ``reference`` as the reference.

    Uses only the statistics present in every seed and architecture. Raises ValueError when there are
    fewer than two architectures or a t is not finite (a contrast with no variance across seeds), since
    a NaN statistic would compare False against every draw and report the smallest possible p.
    """
    archs = list(parts[0]["architectures"])
    m = len(archs)
    if m < 2:
        raise ValueError(f"the joint test needs two or more architectures, got {archs}")
    names = [
        k
        for k in gaussian_names(parts, reference)
        if all(p["architectures"][a]["info"].get(k) is not None for p in parts for a in archs)
    ]
    x = np.array([[[p["architectures"][a]["info"][k] for k in names] for a in archs] for p in parts])
    contrast = (x * m - x.sum(axis=1, keepdims=True)) / (m - 1)  # [seed, arch, stat]: arch minus the mean of the rest

    def statistic(d: np.ndarray) -> np.ndarray:  # d: [..., seed, stat]
        t = one_sample_t(d)
        if not np.isfinite(t).all():
            bad = sorted({names[i] for i in np.nonzero(~np.isfinite(t))[-1]})
            raise ValueError(f"joint test: t is not finite for {bad} (no variance across seeds)")
        return (t**2).sum(axis=-1)

    seeds = np.arange(len(parts))
    observed = float(statistic(contrast[:, archs.index(reference)]))
    rng = np.random.default_rng(seed)
    hits = 0
    for start in range(0, draws, 1000):
        pick = rng.integers(m, size=(min(1000, draws - start), len(parts)))
        hits += int((statistic(contrast[seeds, pick]) >= observed).sum())
    return observed, (1 + hits) / (1 + draws), names


def main() -> int:
    """Print the contrast table and its summary for one case."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--case", default="simace_ts_p2")
    ap.add_argument("--reference", default="c1m_a05")
    ap.add_argument("--mde-k", type=int, help="also print scored rows' minimum detectable gap at alpha / K")
    ap.add_argument("--draws", type=int, default=20000, help="permutation draws of the joint test")
    ap.add_argument("--parts-dir", type=Path, default=PARTS, help="directory holding <case>/seed*.json")
    ap.add_argument(
        "--accept-part-simace",
        action="store_true",
        help="accept the simace version the Parts record instead of the installed one (INFO use only)",
    )
    args = ap.parse_args()
    case = {c["id"]: c for c in yaml.safe_load((ROOT / "cases.yaml").read_text())["cases"]}[args.case]
    simace = None
    if args.accept_part_simace:
        first = args.parts_dir / case["id"] / f"seed{case['seed0']}.json"
        simace = json.loads(first.read_text())["manifest"]["simace"]
    parts = load_parts(case, args.parts_dir, simace)

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

    archs = list(parts[0]["architectures"])
    print()
    print(
        f"joint test, each architecture against the mean of the others ({len(parts)} seeds, {args.draws} draws;"
        " a small p for one reference can come from any architecture that departs, so read the column):"
    )
    for ref in [args.reference, *(a for a in archs if a != args.reference)]:
        observed, p, names = joint_test(parts, ref, draws=args.draws)
        print(f"  {ref}: sum of squared t over {len(names)} statistics = {observed:.2f}, permutation p = {p:.4f}")

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
