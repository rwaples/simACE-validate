"""Layer-5 statistics of one simACE pedigree, by generation (default env).

simACE's recorded generation k follows k + 1 rounds of mating (ADR 0001 §3),
so it is reported as g = k + 1, matching xftsim's generation g (the offspring
of g matings). simACE does not record founders, so g = 0 has no row.
"""

from __future__ import annotations

import numpy as np
import polars as pl

from am.stats import cross_trait_stats, generation_stats


def _values(frame: pl.DataFrame, traits: int) -> list[tuple[np.ndarray, np.ndarray]]:
    """(A, liability) of each scored trait."""
    return [(frame[f"A{t}"].to_numpy(), frame[f"liability{t}"].to_numpy()) for t in range(1, traits + 1)]


def _suffixed(d: dict, t: int) -> dict:
    """Trait 1 keys unchanged; trait 2 keys get `_t2`, as in am/run_xft.py."""
    return d if t == 0 else {f"{k}_t{t + 1}": v for k, v in d.items()}


def pedigree_rows(ped: pl.DataFrame, traits: int) -> list[dict]:
    """Per-generation statistics of one pedigree, one dict per g (with key `generation`)."""
    ped = ped.sort("id")
    rows = []
    prev = None
    for k in range(int(ped["generation"].max()) + 1):
        cur = ped.filter(pl.col("generation") == k)
        v = _values(cur, traits)
        row: dict = {"generation": k + 1}
        if prev is None:
            for t, (a, y) in enumerate(v):
                row.update(_suffixed({"var_A": float(np.var(a, ddof=1)), "var_L": float(np.var(y, ddof=1))}, t))
        else:
            pid = prev["id"].to_numpy()
            mother, father = cur["mother"].to_numpy(), cur["father"].to_numpy()
            pos, fos = np.searchsorted(pid, mother), np.searchsorted(pid, father)
            if not ((pid[pos] == mother).all() and (pid[fos] == father).all()):
                raise ValueError(f"generation {k}: a parent is not in generation {k - 1}")
            pv = _values(prev, traits)
            for t in range(traits):
                row.update(_suffixed(generation_stats(v[t][0], v[t][1], pos, fos, pv[t][0], pv[t][1]), t))
            if traits == 2:
                row.update(cross_trait_stats(v, pv, pos, fos))
        rows.append(row)
        prev = cur
    return rows
