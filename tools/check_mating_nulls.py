"""Check structure/mating.py's null samplers and invariants.

    pixi run python -P tools/check_mating_nulls.py

1. surjection_draws against the exact law, enumerated on a small case.
2. Each Part A null against simACE's own `_mating_standard` on one parent
   generation: the per-generation PITs must be U(0, 1). This is a consistency
   check of the derivations, not evidence about simACE (the harness is).
3. The MZ invariants raise on hand-broken pedigrees.
Exits 1 if a test p-value is below 0.001 or an invariant does not raise.
"""

from __future__ import annotations

import itertools
import math
import sys
from pathlib import Path

import numpy as np
import polars as pl
from scipy import stats
from simace.simulation.simulate import _mating_standard

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from structure.mating import (
    Generation,
    allocation_draws,
    binomial_pit,
    generations,
    n_multi,
    observed_allocation,
    occurrences,
    partner_draws,
    rank_pit,
    surjection_draws,
)

P_MIN = 1e-3


def exact_surjection(n: int, w: np.ndarray) -> dict[tuple[int, ...], float]:
    """Return P(c) of Multinomial(n, w/sum w) conditioned on every count >= 1, by enumeration."""
    probs = {}
    for c in itertools.product(range(1, n + 1), repeat=len(w)):
        if sum(c) == n:
            probs[c] = math.factorial(n) * math.prod(wi**ci / math.factorial(ci) for wi, ci in zip(w, c, strict=True))
    total = sum(probs.values())
    return {c: p / total for c, p in probs.items()}


def check_surjection(rng: np.random.Generator) -> list[tuple[str, float]]:
    """Compare sampled count vectors with the enumerated law by chi-square."""
    out = []
    for w in (np.ones(3), np.array([1.0, 2.0, 3.0])):
        exact = exact_surjection(7, w)
        draws = surjection_draws(rng, 7, w, 100_000)
        keys = list(exact)
        seen = dict.fromkeys(keys, 0)
        for row in map(tuple, draws):
            seen[row] += 1
        obs = np.array([seen[k] for k in keys])
        exp = np.array([exact[k] for k in keys]) * len(draws)
        out.append((f"surjection n=7 w={w.tolist()}", float(stats.chisquare(obs, exp).pvalue)))
    return out


def simace_generation(rng: np.random.Generator, n: int, p_mz: float) -> Generation:
    """Run simACE's standard mating on one parent generation and read it as a Generation."""
    sex = rng.binomial(1, 0.5, size=n)
    parents, twins, _ = _mating_standard(rng, sex, 0.5, p_mz)
    keys, counts = np.unique(parents.T, axis=1, return_counts=True)
    kid_sex = rng.binomial(1, 0.5, size=n)
    kid_sex[twins[:, 1]] = kid_sex[twins[:, 0]]
    return Generation(
        n_female=int(np.sum(sex == 0)),
        n_male=int(np.sum(sex == 1)),
        n_kids=n,
        n_mz=len(twins),
        males=int(np.sum(kid_sex == 1) - np.sum(kid_sex[twins[:, 1]] == 1)),
        counts=counts,
        couple_mother=keys[0],
        couple_father=keys[1],
        mother_couples=occurrences(keys[0]),
    )


def check_nulls(rng: np.random.Generator, reps: int = 300, n: int = 500) -> list[tuple[str, float]]:
    """KS-test the per-generation PIT of every Part A row on simACE's own mating."""
    pits: dict[str, list[float]] = {}
    for _ in range(reps):
        gen = simace_generation(rng, n, 0.05)
        pits.setdefault("sex_males", []).append(binomial_pit(gen.males, n - gen.n_mz, 0.5, rng)[0])
        pits.setdefault("mz_pairs", []).append(binomial_pit(gen.n_mz, int(np.sum(gen.counts >= 2)), 0.05, rng)[0])
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
        (f"null {row} vs simACE (KS of {reps} PITs)", float(stats.kstest(u, "uniform").pvalue))
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


def broken_pedigrees() -> dict[str, pl.DataFrame]:
    """Return small pedigrees that each break one parent or MZ invariant."""
    cases = {
        "an MZ pair differs in sex": {"twin": [-1, -1, -1, -1, -1, -1, 7, 6], "sex": [0, 0, 1, 1, 0, 0, 0, 1]},
        "an MZ pair has different parents": {"twin": [-1, -1, -1, -1, -1, 6, 5, -1], "sex": [0, 0, 1, 1, 0, 1, 1, 1]},
        "a couple has two MZ pairs": {"twin": [-1, -1, -1, -1, 5, 4, 7, 6], "mother": [-1, -1, -1, -1, 0, 0, 0, 0]}
        | {"father": [-1, -1, -1, -1, 2, 2, 2, 2], "sex": [0, 0, 1, 1, 0, 0, 1, 1]},
        "not reciprocal": {"twin": [-1, -1, -1, -1, 5, 6, -1, -1]},
        "a mother is not a generation-0 individual of that sex": {"mother": [-1, -1, -1, -1, 0, 0, 2, 2]},
        "a father is not a generation-0 individual of that sex": {"father": [-1, -1, -1, -1, 2, 2, 5, 5]},
        "generation 0: an MZ pair differs in sex": {"twin": [-1, 2, 1, -1, 5, 4, -1, -1]},
    }
    return {msg: pl.DataFrame(VALID | change) for msg, change in cases.items()}


def main() -> int:
    """Run every check and print one line each."""
    rng = np.random.default_rng(20261022)
    failed = False
    generations(pl.DataFrame(VALID))  # raises if a valid pedigree is rejected
    for name, p in [*check_surjection(rng), *check_nulls(rng)]:
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
