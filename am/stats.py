"""Per-generation statistics for the assortative-mating pilot (issue #4).

One function computes every statistic on simACE's and on xftsim's values, so
the two simulators are compared on the same definitions. numpy only: it runs
in both the default env (Python 3.14) and the xftsim env (Python 3.9).
"""

from __future__ import annotations

import numpy as np


def _corr(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def _sib_pairs(couple: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Every unordered pair of offspring that share a couple id."""
    order = np.argsort(couple, kind="stable")
    c = couple[order]
    starts = np.flatnonzero(np.r_[True, c[1:] != c[:-1]])
    sizes = np.diff(np.r_[starts, len(c)])
    i_list, j_list = [], []
    for s, k in zip(starts[sizes >= 2], sizes[sizes >= 2]):  # noqa: B905  (strict= needs Python 3.10; this runs on 3.9)
        idx = order[s : s + k]
        ii, jj = np.triu_indices(k, 1)
        i_list.append(idx[ii])
        j_list.append(idx[jj])
    if not i_list:
        return np.empty(0, np.int64), np.empty(0, np.int64)
    return np.concatenate(i_list), np.concatenate(j_list)


def _sym_corr(x: np.ndarray, i: np.ndarray, j: np.ndarray) -> float:
    """Pearson correlation over pairs entered in both orders (intraclass-like)."""
    return _corr(np.r_[x[i], x[j]], np.r_[x[j], x[i]])


def generation_stats(
    a: np.ndarray,
    y: np.ndarray,
    mother: np.ndarray,
    father: np.ndarray,
    parent_a: np.ndarray,
    parent_y: np.ndarray,
) -> dict[str, float]:
    """Statistics of one offspring generation.

    a, y: offspring additive value and phenotype (liability).
    mother, father: each offspring's parents, as row indices into parent_a
        and parent_y (the previous generation).
    """
    a = np.asarray(a, np.float64)
    y = np.asarray(y, np.float64)
    pa = np.asarray(parent_a, np.float64)
    py = np.asarray(parent_y, np.float64)
    mother = np.asarray(mother, np.int64)
    father = np.asarray(father, np.int64)

    couple_key = mother * (len(pa) + 1) + father
    uniq, first = np.unique(couple_key, return_index=True)
    cm, cf = mother[first], father[first]
    _, couple = np.unique(couple_key, return_inverse=True)
    fi, fj = _sib_pairs(couple)

    return {
        "n": float(len(a)),
        "n_couples": float(len(uniq)),
        "n_fs_pairs": float(len(fi)),
        "var_A": float(np.var(a, ddof=1)),
        "var_L": float(np.var(y, ddof=1)),
        # mate correlations over distinct couples, in the parents' generation
        "r_mate_L": _corr(py[cm], py[cf]),
        "r_mate_A": _corr(pa[cm], pa[cf]),
        "r_FS_A": _sym_corr(a, fi, fj),
        "r_FS_L": _sym_corr(y, fi, fj),
        "r_MO_L": _corr(py[mother], y),
        "r_FO_L": _corr(py[father], y),
        "r_PO_L": _corr(np.r_[py[mother], py[father]], np.r_[y, y]),
        "r_PO_A": _corr(np.r_[pa[mother], pa[father]], np.r_[a, a]),
    }


def bulmer(v0: float, s: float, r: float, g: int) -> list[float]:
    """Var(A) after each of g rounds of mating on liability with mate correlation r.

    Infinitesimal model, no inbreeding: V' = V/2 * (1 + r V / (V + S)) + v0/2,
    with S the non-additive liability variance (C + E). Mates' A correlation is
    r V / V_L. Element t is Var(A) after t matings (element 0 is v0).
    """
    out = [v0]
    v = v0
    for _ in range(g):
        v = v / 2 * (1 + r * v / (v + s)) + v0 / 2
        out.append(v)
    return out
