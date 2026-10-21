"""Scores for the harness's rows (ADR 0001).

Layer 1 (conditional calibration). For one pedigree, an oracle returns how
many of its M draws fall below (``n_less``) and equal to (``n_equal``)
simACE's value. If simACE samples from the oracle's law, the randomised PIT

    u = (n_less + V * (n_equal + 1)) / (M + 1),   V ~ U(0, 1)

is exactly U(0, 1), ties included. Across seeds the u values are i.i.d., and
two tests score them: Kolmogorov-Smirnov (exact) and a chi-square on equal
bins.

Layer 2 (one gene-drop per pedigree). Each seed gives one difference whose
expectation is 0 under the null, and a one-sample t-test scores the seeds.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats

CHI2_BINS = 10


@dataclass(frozen=True)
class Score:
    """The tests of one row. Every p-value counts toward K."""

    n: int
    center: float  # mean PIT (calibration) or mean difference (t-test)
    effect: float  # mean z of simACE against the draws, or the t statistic
    p_values: dict[str, float]

    @property
    def p(self) -> float:
        """Return the smallest p-value of the row's tests."""
        return min(self.p_values.values())


def randomised_pit(n_less: np.ndarray, n_equal: np.ndarray, M: int, rng: np.random.Generator) -> np.ndarray:
    """Return the randomised PIT of each rank against M oracle draws."""
    v = rng.uniform(size=len(n_less))
    return (n_less + v * (n_equal + 1)) / (M + 1)


def calibrate(u: np.ndarray, z: np.ndarray) -> Score:
    """Score PIT values u over seeds; z is simACE's per-seed value in oracle SDs."""
    counts = np.histogram(u, bins=CHI2_BINS, range=(0.0, 1.0))[0]
    return Score(
        n=len(u),
        center=float(u.mean()),
        effect=float(z.mean()),
        p_values={
            "KS": float(stats.kstest(u, "uniform", method="exact").pvalue),
            "chi2": float(stats.chisquare(counts).pvalue),
        },
    )


def mean_zero(d: np.ndarray) -> Score:
    """Score per-seed differences whose expectation is 0 by a one-sample t-test."""
    test = stats.ttest_1samp(d, 0.0)
    return Score(n=len(d), center=float(d.mean()), effect=float(test.statistic), p_values={"t": float(test.pvalue)})


def two_sample(a: np.ndarray, b: np.ndarray) -> Score:
    """Score simACE's per-pedigree values a against a comparator's b: two-sample KS and Welch t."""
    if np.ptp(a) == 0 and np.ptp(b) == 0:
        p = 1.0 if a[0] == b[0] else 0.0  # both constant: they agree or they do not
        return Score(n=len(a), center=float(a.mean() - b.mean()), effect=np.nan, p_values={"KS": p})
    t = stats.ttest_ind(a, b, equal_var=False)
    return Score(
        n=len(a),
        center=float(a.mean() - b.mean()),
        effect=float(t.statistic),
        p_values={"KS": float(stats.ks_2samp(a, b).pvalue), "t": float(t.pvalue)},
    )
