"""Scores for the harness's rows (ADR 0001).

Layer 1 (conditional calibration). For one pedigree, an oracle returns how
many of its M draws fall below (``n_less``) and equal to (``n_equal``)
simACE's value. If simACE samples from the oracle's law, the randomised PIT

    u = (n_less + V * (n_equal + 1)) / (M + 1),   V ~ U(0, 1)

is exactly U(0, 1), ties included. Across seeds the u values are i.i.d., and
two tests score them: Kolmogorov-Smirnov (exact) and a Pearson chi-square on
equal bins, referred to its exact multinomial null (``chi2_exact_sf``), not
the asymptotic chi2(bins - 1), which is anti-conservative at alpha/K.

Layer 2 (one gene-drop per pedigree). Each seed gives one difference whose
expectation is 0 under the null, and a one-sample t-test scores the seeds.

Layer 3b (std_mating). A discrete statistic x with an exact null gets the
randomised PIT F(x - 1) + V * P(x), also exactly U(0, 1); the per-generation
PITs are combined per pedigree in structure/mating.py and scored like layer 1.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from functools import cache

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


def randomised_pit_discrete(
    cdf_below: np.ndarray | float, pmf_at: np.ndarray | float, rng: np.random.Generator
) -> np.ndarray:
    """Return the randomised PIT F(x - 1) + V * P(x) of a discrete statistic x under its exact null."""
    return cdf_below + rng.uniform(size=np.shape(pmf_at)) * pmf_at


@cache
def chi2_null_tail(R: int, bins: int) -> np.ndarray:
    """Return P(S >= s) for s = 0 .. R^2, read-only, where S is the sum of squared counts of multinomial(R, 1/bins).

    A DP over bins in probabilities: g[m, s] is the probability that m of the
    R draws are left for the later bins and the earlier bins' squared counts
    sum to s. A bin with b bins still open takes c ~ Binomial(m, 1/b), and the
    last bin takes the remaining m. Every entry is a probability, at least
    bins^-R, so it stays a normal float while R * log10(bins) <= 300.
    """
    if R < 1 or R * math.log10(bins) > 300:
        raise ValueError(f"exact chi-square null needs 1 <= R and R * log10(bins) <= 300, got R = {R}, bins = {bins}")
    s_max = R * R
    m = np.arange(R + 1)
    g = np.zeros((R + 1, s_max + 1))
    g[R, 0] = 1.0
    for open_bins in range(bins, 1, -1):
        new = np.zeros_like(g)
        for c in range(R + 1):
            pmf = stats.binom.pmf(c, m[c:], 1.0 / open_bins)  # P(c of the m left fall in this bin)
            new[: R + 1 - c, c * c :] += pmf[:, None] * g[c:, : s_max + 1 - c * c]
        g = new
    pmf_s = np.zeros(s_max + 1)
    for left in range(R + 1):
        pmf_s[left * left :] += g[left, : s_max + 1 - left * left]
    tail = np.minimum(np.cumsum(pmf_s[::-1])[::-1], 1.0)
    tail.flags.writeable = False
    return tail


def chi2_exact_sf(counts: np.ndarray) -> float:
    """Return the exact p-value of the Pearson equal-bin chi-square of counts under multinomial(sum, 1/bins).

    The statistic is (bins / R) * sum(counts^2) - R, a function of the sum of
    squares alone, so its upper tail is that of the sum of squares. With no
    counts (R = 0) the statistic is undefined and the p-value is NaN, which
    makes the row an ERROR.
    """
    counts = np.asarray(counts, dtype=np.int64)
    R = int(counts.sum())
    if R == 0:
        return float("nan")
    s = int(np.dot(counts, counts))
    if counts.min() < 0 or not R * R // len(counts) <= s <= R * R:
        raise ValueError(f"not a vector of bin counts: {counts}")
    return float(chi2_null_tail(R, len(counts))[s])


def calibrate(u: np.ndarray, z: np.ndarray) -> Score:
    """Score PIT values u over seeds; z is simACE's per-seed value in oracle SDs."""
    counts = np.histogram(u, bins=CHI2_BINS, range=(0.0, 1.0))[0]
    return Score(
        n=len(u),
        center=float(u.mean()),
        effect=float(z.mean()),
        p_values={
            "KS": float(stats.kstest(u, "uniform", method="exact").pvalue),
            "chi2": chi2_exact_sf(counts),
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


def welch(a: np.ndarray, b: np.ndarray) -> Score:
    """Score simACE's per-seed values a against a reference sample's b by a Welch t alone (layer 5, ADR 0003).

    The two simulators' spreads differ by design, so a distribution test (KS) would reject with no bug.
    """
    t = stats.ttest_ind(a, b, equal_var=False)
    return Score(
        n=len(a), center=float(a.mean() - b.mean()), effect=float(t.statistic), p_values={"t": float(t.pvalue)}
    )
