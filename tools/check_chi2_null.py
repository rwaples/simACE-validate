"""Check calibration.chi2_exact_sf: against brute-force enumeration, and its simulated size at alpha/K (G4).

Run: pixi run python -P tools/check_chi2_null.py (2 x 1e8 null draws).
"""

import itertools
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from calibration import CHI2_BINS, _sum_sq_tail, calibrate, chi2_exact_sf

ALPHA_K = 0.05 / 1364


def brute(R: int, k: int) -> None:
    """Enumerate every count vector and compare each exact tail with the DP."""
    worst = 0.0
    vecs = [c for c in itertools.product(range(R + 1), repeat=k) if sum(c) == R]
    prob = {c: math.factorial(R) / math.prod(math.factorial(x) for x in c) / k**R for c in vecs}
    for c in vecs:
        s = sum(x * x for x in c)
        p_bf = sum(q for d, q in prob.items() if sum(x * x for x in d) >= s)
        worst = max(worst, abs(p_bf - chi2_exact_sf(np.array(c))))
    print(f"brute force R={R} bins={k}: {len(vecs)} vectors, max |DP - enumeration| = {worst:.3g}")


def simulate(R: int, n_sim: int, seed: int) -> None:
    """Simulated rejection rate of the exact and asymptotic p at ALPHA_K under uniform PITs."""
    rng = np.random.default_rng(seed)
    tail = _sum_sq_tail(R, CHI2_BINS)
    hits_exact = hits_asym = 0
    chunk = 5_000_000
    for start in range(0, n_sim, chunk):
        m = min(chunk, n_sim - start)
        counts = rng.multinomial(R, np.full(CHI2_BINS, 1 / CHI2_BINS), size=m)
        s = np.sum(counts.astype(np.int64) ** 2, axis=1)
        hits_exact += int(np.sum(tail[s] <= ALPHA_K))
        hits_asym += int(np.sum(stats.chi2.sf(CHI2_BINS / R * s - R, CHI2_BINS - 1) <= ALPHA_K))
    exact_size = float(np.sum(np.diff(-tail, append=0.0)[tail <= ALPHA_K]))  # P(S = s) summed where p <= alpha
    for name, h in (("exact", hits_exact), ("asymptotic", hits_asym)):
        lo, hi = stats.binomtest(h, n_sim).proportion_ci(0.95, method="exact")
        print(f"R={R} {name}: {h}/{n_sim} = {h / n_sim:.4g} (95% CI {lo:.4g} to {hi:.4g}); nominal {ALPHA_K:.4g}")
    print(f"R={R} exact size from the pmf: {exact_size:.4g}")


def end_to_end(R: int, n_sim: int, seed: int) -> None:
    """calibrate() on uniform PITs: the chi2 p is valid at the conventional 0.05 too."""
    rng = np.random.default_rng(seed)
    p = np.array([calibrate(rng.uniform(size=R), np.zeros(R)).p_values["chi2"] for _ in range(n_sim)])
    print(
        f"calibrate R={R}: {n_sim} uniform samples, rejection at 0.05 = {np.mean(p <= 0.05):.4g}, at 0.01 = {np.mean(p <= 0.01):.4g}"
    )


def main() -> None:
    """Run the enumeration checks, then the size simulation at R = 50 and 100."""
    brute(8, 4)
    brute(12, 3)
    for R in (50, 100):
        simulate(R, 100_000_000, seed=20261108 + R)
    end_to_end(50, 20_000, seed=1)


if __name__ == "__main__":
    main()
