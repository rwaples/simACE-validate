"""Count duplicate couples left after simACE's swap passes, with simACE's own functions.

    pixi run python -P tools/residual_duplicates.py --N 2000 --gens 200000

A duplicate that survives `pair_partners`' 5 swap passes merges two matings
into one visible couple, which breaks the exactness of std_mating's
partner-count, allocation and twin nulls (structure/mating.py). Parent sexes
are Bernoulli(1/2), as in simACE.
"""

from __future__ import annotations

import argparse

import numpy as np
from simace.simulation.simulate import (
    _find_duplicate_pairs,
    balance_mating_slots,
    draw_mating_counts,
    pair_partners,
)


def main() -> None:
    """Simulate one mating step per generation and count the duplicates left."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--N", type=int, default=2000)
    parser.add_argument("--gens", type=int, default=200_000)
    parser.add_argument("--lam", type=float, default=0.5)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    left = 0
    gens_with = 0
    for _ in range(args.gens):
        n_male = rng.binomial(args.N, 0.5)
        males, females = np.arange(n_male), np.arange(n_male, args.N)
        m, f = balance_mating_slots(
            rng, draw_mating_counts(rng, len(males), args.lam), draw_mating_counts(rng, len(females), args.lam)
        )
        k = int(_find_duplicate_pairs(pair_partners(rng, males, m, females, f)).sum())
        left += k
        gens_with += k > 0
    bound = 3 / args.gens if gens_with == 0 else float("nan")
    print(
        f"N={args.N} lambda={args.lam} generations={args.gens}: residual duplicates {left}, "
        f"generations with one {gens_with}"
        + (f", 95% upper bound per generation {bound:.2g}" if gens_with == 0 else "")
    )


if __name__ == "__main__":
    main()
