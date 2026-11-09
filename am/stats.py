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


def _half_sib_pairs(shared: np.ndarray, other: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Every unordered pair of offspring with the same `shared` parent and different `other` parents."""
    i, j = _sib_pairs(np.asarray(shared, np.int64))
    keep = other[i] != other[j]
    return i[keep], j[keep]


def ctm_stats(
    c: np.ndarray,
    y: np.ndarray,
    mother: np.ndarray | None = None,
    father: np.ndarray | None = None,
    parent_c: np.ndarray | None = None,
    parent_y: np.ndarray | None = None,
) -> dict[str, float]:
    """C-transmission statistics of one offspring generation (rwaples/simACE-validate#5).

    c, y: offspring C and phenotype (liability). Without parents (simACE's first recorded
    generation) only Var(C) is defined. MHS share the mother and not the father, PHS the
    reverse; b_MP_L is the slope of offspring L on midparent L, (A + (b_m + b_f) C) / V_P
    under random mating.
    """
    c = np.asarray(c, np.float64)
    out = {"var_C": float(np.var(c, ddof=1))}
    if mother is None:
        return out
    y = np.asarray(y, np.float64)
    pc, py = np.asarray(parent_c, np.float64), np.asarray(parent_y, np.float64)
    mother, father = np.asarray(mother, np.int64), np.asarray(father, np.int64)
    mi, mj = _half_sib_pairs(mother, father)
    pi, pj = _half_sib_pairs(father, mother)
    mid = (py[mother] + py[father]) / 2
    out |= {
        "r_MO_C": _corr(pc[mother], c),
        "r_FO_C": _corr(pc[father], c),
        "r_MHS_C": _sym_corr(c, mi, mj),
        "r_PHS_C": _sym_corr(c, pi, pj),
        "r_MHS_L": _sym_corr(y, mi, mj),
        "r_PHS_L": _sym_corr(y, pi, pj),
        "n_mhs_pairs": float(len(mi)),
        "n_phs_pairs": float(len(pi)),
        "b_MP_L": float(np.cov(mid, y)[0, 1] / np.var(mid, ddof=1)),
    }
    return out


def cross_trait_stats(
    v: list[tuple[np.ndarray, np.ndarray]],
    pv: list[tuple[np.ndarray, np.ndarray]],
    mother: np.ndarray,
    father: np.ndarray,
) -> dict[str, float]:
    """Cross-trait statistics of one offspring generation, for two-trait AM.

    v, pv: [(A, L) of trait 1, (A, L) of trait 2] for the offspring and for
        their parents' generation; mother, father index into pv's rows.
    Mate correlations are over distinct couples, mother's trait first: r_mate_L12
    pairs the mother's L1 with the father's L2. r_w_* is the within-person
    cross-trait correlation in the offspring generation.
    """
    mother = np.asarray(mother, np.int64)
    father = np.asarray(father, np.int64)
    (a1, y1), (a2, y2) = ([np.asarray(x, np.float64) for x in t] for t in v)
    (pa1, py1), (pa2, py2) = ([np.asarray(x, np.float64) for x in t] for t in pv)
    key = mother * (len(py1) + 1) + father
    _, first = np.unique(key, return_index=True)
    cm, cf = mother[first], father[first]
    return {
        "r_w_A": _corr(a1, a2),
        "r_w_L": _corr(y1, y2),
        "r_w_L_parents": _corr(py1, py2),
        "r_mate_L12": _corr(py1[cm], py2[cf]),
        "r_mate_L21": _corr(py2[cm], py1[cf]),
        "r_mate_A12": _corr(pa1[cm], pa2[cf]),
        "r_mate_A21": _corr(pa2[cm], pa1[cf]),
    }


def flatten(rows: list[dict]) -> dict[str, float]:
    """One run's per-generation rows as {stat_g<g>: value}; generation 0 (founders) is dropped."""
    return {
        f"{k}_g{r['generation']}": v for r in rows if r["generation"] > 0 for k, v in r.items() if k != "generation"
    }
