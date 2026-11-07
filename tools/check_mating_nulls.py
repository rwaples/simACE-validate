"""Check structure/mating.py's null samplers and invariants.

    pixi run python -P tools/check_mating_nulls.py

1. surjection_draws against the exact law, enumerated on small cases
   (equal weights, unequal weights with ties, Dirichlet-multinomial), and
   at a simACE generation's K ~ 1,000 couples against the pre-#6 rejection
   sampler, an independent method, on three summaries.
2. Each Part A null against simACE's own `_mating_standard` on one parent
   generation, read through `generations` (so its invariants hold on
   simACE's output): the per-generation PITs must be U(0, 1). Run at the
   case's p_mz and at p_mz = 0.9, where couples carry several MZ pairs and
   the N // 2 cap on t binds. This is a consistency check of the
   derivations, not evidence about simACE (the harness is).
3. The MZ and order invariants raise on hand-broken pedigrees and accept
   valid ones, including a couple with two MZ pairs.
Exits 1 if a test p-value is below 0.001 or an invariant check is wrong.
"""

from __future__ import annotations

import itertools
import math
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import optimize, stats
from simace.simulation.simulate import _mating_standard

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from structure.mating import (
    DIRICHLET_SHAPE,
    allocation_draws,
    binomial_pit,
    cell_pmf,
    generations,
    mz_rate,
    n_multi,
    observed_allocation,
    partner_draws,
    rank_pit,
    split_sums,
    surjection_draws,
    ztp,
)

P_MIN = 1e-3


def exact_surjection(n: int, w: np.ndarray, shape: float | None = None) -> dict[tuple[int, ...], float]:
    """Return P(c) of Multinomial(n, w/sum w), or with shape Dirichlet(shape)-multinomial, conditioned on every count >= 1."""
    probs = {}
    for c in itertools.product(range(1, n + 1), repeat=len(w)):
        if sum(c) == n:
            if shape is None:
                probs[c] = math.prod(wi**ci / math.factorial(ci) for wi, ci in zip(w, c, strict=True))
            else:
                probs[c] = math.prod(math.gamma(ci + shape) / math.factorial(ci) for ci in c)
    total = sum(probs.values())
    return {c: p / total for c, p in probs.items()}


SURJECTION_CASES = [  # (n, weights, dirichlet shape): K = 5 makes uneven halves and a third level; ties share a law
    (7, np.ones(3), None),
    (7, np.array([1.0, 2.0, 3.0]), None),
    (10, np.ones(5), None),
    (10, np.array([2.0, 1.0, 3.0, 1.0, 2.0]), None),
    (9, np.ones(4), 2.0),
    (9, np.ones(4), DIRICHLET_SHAPE),
]


def check_surjection(rng: np.random.Generator) -> list[tuple[str, float]]:
    """Compare sampled count vectors with the enumerated law by chi-square."""
    out = []
    for n, w, shape in SURJECTION_CASES:
        exact = exact_surjection(n, w, shape)
        draws = surjection_draws(rng, n, w, 100_000, dirichlet_shape=shape)
        keys = list(exact)
        seen = dict.fromkeys(keys, 0)
        for row in map(tuple, draws):
            seen[row] += 1
        obs = np.array([seen[k] for k in keys])
        exp = np.array([exact[k] for k in keys]) * len(draws)
        label = f"surjection n={n} w={w.tolist()}" + (f" Dirichlet({shape:g})" if shape else "")
        out.append(
            (
                f"{label} (chi-square, {len(keys)} cells, min expected {exp.min():.0f})",
                float(stats.chisquare(obs, exp).pvalue),
            )
        )
    return out


def rejection_surjection(rng: np.random.Generator, n: int, w: np.ndarray, m: int) -> np.ndarray:
    """Return m draws of surjection_draws' law by an independent method: K - 1 ZTP counts, the last cell's pmf accepts.

    The counts are ZTP(mu * w_k); the last is n minus the rest, kept with
    probability p(last) / max p for its Poisson pmf p. The pre-#6 sampler.
    """
    k = len(w)
    order = np.argsort(w, kind="stable")
    lam = optimize.brentq(lambda s: np.sum(s * w / -np.expm1(-s * w)) - n, 1e-12, 1e6) * w[order]
    mode = max(1.0, np.floor(lam[-1]))
    log_max = max(stats.poisson.logpmf(mode, lam[-1]), stats.poisson.logpmf(mode + 1, lam[-1]))
    batch = max(64, (1 << 22) // k)
    got, n_got = [], 0
    while n_got < m:
        x = ztp(rng, np.broadcast_to(lam[:-1], (batch, k - 1)), (batch, k - 1))
        last = n - x.sum(axis=1)
        log_p = stats.poisson.logpmf(np.maximum(last, 1), lam[-1])
        keep = (last >= 1) & (np.log(rng.uniform(size=batch)) < log_p - log_max)
        got.append(np.column_stack([x[keep], last[keep]]))
        n_got += int(keep.sum())
    draws = np.empty((m, k), dtype=np.int64)
    draws[:, order] = np.concatenate(got)[:m]
    return draws


def two_sample_chisquare(x: np.ndarray, y: np.ndarray, bins: int = 20) -> float:
    """Return the chi-square p-value of x and y counted in bins cut at the pooled sample's quantiles.

    The cut points are values, so tied (integer) data stay in one bin and the
    test keeps its nominal size, unlike a two-sample KS on ties.
    """
    cuts = np.unique(np.quantile(np.concatenate([x, y]), np.linspace(0, 1, bins + 1)[1:-1]))
    table = [np.bincount(np.searchsorted(cuts, v, side="right"), minlength=len(cuts) + 1) for v in (x, y)]
    table = np.array(table)
    return float(stats.chi2_contingency(table[:, table.sum(axis=0) > 0]).pvalue)


def check_surjection_large(rng: np.random.Generator, m: int = 4000) -> list[tuple[str, float]]:
    """Two-sample chi-square of surjection_draws against rejection_surjection on a simACE generation at N = 2000."""
    (gen,) = generations(simace_pedigree(rng, 2000, 0.02))
    k, n = len(gen.counts), gen.n_kids - gen.n_mz
    multi = gen.mother_couples >= 2
    out = []
    for label, w in (("equal weights", np.ones(k)), ("weights 1/(mother's couples)", 1.0 / gen.mother_couples)):
        tree, ref = surjection_draws(rng, n, w, m), rejection_surjection(rng, n, w, m)
        for stat, f in (
            ("sum c^2", lambda d: (d**2).sum(axis=1)),
            ("couples with c = 1", lambda d: (d == 1).sum(axis=1)),
            ("mean c over multi-couple mothers", lambda d: d[:, multi].mean(axis=1)),
        ):
            p = two_sample_chisquare(f(tree), f(ref))
            out.append((f"surjection K={k} n={n} {label}: {stat} vs rejection (binned chi-square, {m} + {m})", p))
    return out


def simace_pedigree(rng: np.random.Generator, n: int, p_mz: float) -> pl.DataFrame:
    """Run simACE's standard mating and reproduction sex rule on one parent generation; return both generations."""
    sex = rng.binomial(1, 0.5, size=n)
    parents, twins, _ = _mating_standard(rng, sex, 0.5, p_mz)
    kid_sex = rng.binomial(1, 0.5, size=n)
    kid_sex[twins[:, 1]] = kid_sex[twins[:, 0]]
    twin = np.full(n, -1)
    twin[twins[:, 0]], twin[twins[:, 1]] = twins[:, 1] + n, twins[:, 0] + n
    none = np.full(n, -1)
    return pl.DataFrame(
        {
            "generation": np.repeat([0, 1], n),
            "sex": np.concatenate([sex, kid_sex]),
            "mother": np.concatenate([none, parents[:, 0]]),
            "father": np.concatenate([none, parents[:, 1]]),
            "twin": np.concatenate([none, twin]),
        }
    )


def check_nulls(rng: np.random.Generator, n: int, p_mz: float, reps: int = 300) -> list[tuple[str, float]]:
    """KS-test the per-generation PIT of every Part A row on simACE's own mating."""
    pits: dict[str, list[float]] = {}
    for _ in range(reps):
        (gen,) = generations(simace_pedigree(rng, n, p_mz))
        pits.setdefault("sex_males", []).append(binomial_pit(gen.males, n - gen.n_mz, 0.5, rng)[0])
        pits.setdefault("mz_pairs", []).append(binomial_pit(gen.n_mz, n, mz_rate(p_mz), rng, cap=n // 2)[0])
        sumsq, mean_multi = allocation_draws(rng, gen, 199, None)
        obs_sumsq, obs_multi = observed_allocation(gen)
        pits.setdefault("alloc_sumsq", []).append(rank_pit(obs_sumsq, sumsq, rng)[0])
        if not np.isnan(obs_multi):
            pits.setdefault("alloc_multi_mother", []).append(rank_pit(obs_multi, mean_multi, rng)[0])
        partners = partner_draws(rng, gen, 0.5, 199)
        observed = (len(gen.counts), n_multi(gen.couple_mother), n_multi(gen.couple_father))
        for j, row in enumerate(("n_couples", "n_mothers_multi", "n_fathers_multi")):
            pits.setdefault(row, []).append(rank_pit(observed[j], partners[:, j], rng)[0])
    return [
        (f"null {row} vs simACE at N={n}, p_mz={p_mz} (KS of {len(u)} PITs)", float(stats.kstest(u, "uniform").pvalue))
        for row, u in pits.items()
    ]


# generation 0: ids 0-3 (0, 1 female; 2, 3 male); generation 1: ids 4-7, with one valid MZ pair (4, 5)
VALID = {
    "generation": [0, 0, 0, 0, 1, 1, 1, 1],
    "sex": [0, 0, 1, 1, 0, 0, 1, 1],
    "mother": [-1, -1, -1, -1, 0, 0, 1, 1],
    "father": [-1, -1, -1, -1, 2, 2, 3, 3],
    "twin": [-1, -1, -1, -1, 5, 4, -1, -1],
}


# generation 1 of each: ids 4-7, all children of couple (0, 2)
ONE_COUPLE = {"mother": [-1, -1, -1, -1, 0, 0, 0, 0], "father": [-1, -1, -1, -1, 2, 2, 2, 2]}


def valid_pedigrees() -> dict[str, pl.DataFrame]:
    """Return small pedigrees that every invariant must accept."""
    cases = {
        "one MZ pair": {},
        "a couple with two MZ pairs": ONE_COUPLE | {"twin": [-1, -1, -1, -1, 5, 4, 7, 6]},
        "an MZ pair after a singleton": ONE_COUPLE
        | {"twin": [-1, -1, -1, -1, -1, 6, 5, -1], "sex": [0, 0, 1, 1, 0, 1, 1, 0]},
    }
    return {name: pl.DataFrame(VALID | change) for name, change in cases.items()}


def broken_pedigrees() -> dict[str, pl.DataFrame]:
    """Return small pedigrees that each break one parent, MZ or order invariant."""
    cases = {
        "generation 1: no children": {"generation": [0, 0, 0, 0, 2, 2, 2, 2]},
        "an MZ pair differs in sex": {"twin": [-1, -1, -1, -1, -1, -1, 7, 6], "sex": [0, 0, 1, 1, 0, 0, 0, 1]},
        "an MZ pair has different parents": {"twin": [-1, -1, -1, -1, -1, 6, 5, -1], "sex": [0, 0, 1, 1, 0, 1, 1, 1]},
        "co-twins do not have consecutive ids": ONE_COUPLE
        | {"twin": [-1, -1, -1, -1, 6, -1, 4, -1], "sex": [0, 0, 1, 1, 0, 1, 0, 1]},
        "mothers are not in id order": {"mother": [-1, -1, -1, -1, 1, 1, 0, 0], "father": [-1, -1, -1, -1, 3, 3, 2, 2]},
        "a couple's children are not one block of consecutive ids": {"father": [-1, -1, -1, -1, 2, 3, 2, 3]}
        | {"mother": [-1, -1, -1, -1, 0, 0, 0, 0], "twin": [-1] * 8},
        "not reciprocal": {"twin": [-1, -1, -1, -1, 5, 6, -1, -1]},
        "a mother is not a generation-0 individual of that sex": {"mother": [-1, -1, -1, -1, 0, 0, 2, 2]},
        "a father is not a generation-0 individual of that sex": {"father": [-1, -1, -1, -1, 2, 2, 5, 5]},
        "generation 0: an MZ pair differs in sex": {"twin": [-1, 2, 1, -1, 5, 4, -1, -1]},
    }
    return {msg: pl.DataFrame(VALID | change) for msg, change in cases.items()}


def check_surjection_edges(rng: np.random.Generator) -> list[tuple[str, bool]]:
    """Return (name, passed) for edge inputs that broke the first tree sampler (#6 review)."""
    tiny = surjection_draws(rng, 2049, np.ones(2048), 500, dirichlet_shape=1e-4)
    left, right = cell_pmf(3, 2.0**-535, None), cell_pmf(3, 2.0**-540, None)
    split = split_sums(rng, np.full(100_000, 4), left, right)
    # exact P(left = 2 | 4) for ZTP(2^-535) and ZTP(2^-540) halves: weights 2^-1080 / 6, 2^-1075 / 4, 2^-1070 / 6
    w = np.array([2.0**-10 / 6, 2.0**-5 / 4, 1 / 6])
    p_two = float(stats.binomtest(int(np.sum(split == 2)), len(split), w[1] / w.sum()).pvalue)
    sums, halves = np.arange(2, 400), cell_pmf(399, 2.0, None)
    # widths 1..204, so max_cells = 1000 puts rows of different widths in one chunk and leaves a short last chunk
    whole, *chunked = (
        split_sums(np.random.default_rng(5), sums, halves, halves, max_cells=c) for c in (1 << 20, 1000, 7)
    )
    # the capped z over 10^6 simulated min(Bin(101, q), 50) at p_mz = 0.9, where the cap binds (SE of the mean 0.001)
    q = mz_rate(0.9)
    z = np.array([binomial_pit(v, 101, q, np.random.default_rng(0), cap=50)[1] for v in range(51)])
    z_sim = z[np.minimum(rng.binomial(101, q, size=1_000_000), 50)]
    z_mean, z_sd = float(z_sim.mean()), float(z_sim.std())
    return [
        *(
            (f"split_sums in chunks of {c} products draws what one chunk does", bool(np.array_equal(whole, got)))
            for c, got in zip((1000, 7), chunked, strict=True)
        ),
        (
            f"mz_pairs z over simulated capped t: mean {z_mean:.2g}, SD {z_sd:.4f}",
            abs(z_mean) < 0.005 and abs(z_sd - 1) < 0.005,
        ),
        ("n = K gives all ones", bool(np.all(surjection_draws(rng, 3, np.ones(3), 2) == 1))),
        *(
            (
                f"weights at scale {scale:g} sum to n",
                bool(np.all(surjection_draws(rng, 3, np.full(2, scale), 5).sum(axis=1) == 3)),
            )
            for scale in (1e-12, 1e14)
        ),
        *(  # a long leaf support must be trimmed, and mu's bracket must not cancel at n = K + 1
            (f"n = {n:,}, K = {k:,} sums to n", bool(np.all(surjection_draws(rng, n, np.ones(k), 2).sum(axis=1) == n)))
            for n, k in ((40_000, 4), (50_001, 50_000), (100_001, 100_000))
        ),
        ("Dirichlet(1e-4), n = K + 1: one couple with 2", bool(np.all(np.sort(tiny, axis=1)[:, -2:] == [1, 2]))),
        (f"split of subnormal halves: P(left = 2) binomial p = {p_two:.3g}", p_two >= P_MIN),
    ]


def main() -> int:
    """Run every check and print one line each."""
    rng = np.random.default_rng(20261022)
    failed = False
    for name, ped in valid_pedigrees().items():
        try:
            generations(ped)
        except ValueError as exc:
            print(f"FAIL valid pedigree '{name}' raises: {exc}")
            failed = True
        else:
            print(f"ok   valid pedigree '{name}' passes")
    checks = [
        *check_surjection(rng),
        *check_surjection_large(rng),
        *check_nulls(rng, 500, 0.02),
        *check_nulls(rng, 101, 0.9),
    ]
    for name, ok in check_surjection_edges(rng):
        failed |= not ok
        print(f"{'ok  ' if ok else 'FAIL'} {name}")
    for name, p in checks:
        failed |= p < P_MIN
        print(f"{'FAIL' if p < P_MIN else 'ok  '} {name}: p = {p:.3g}")
    for msg, ped in broken_pedigrees().items():
        try:
            generations(ped)
        except ValueError as exc:
            ok = msg in str(exc)
            print(f"{'ok  ' if ok else 'FAIL'} invariant '{msg}' raises: {exc}")
            failed |= not ok
        else:
            print(f"FAIL invariant '{msg}' did not raise")
            failed = True
    return int(failed)


if __name__ == "__main__":
    sys.exit(main())
