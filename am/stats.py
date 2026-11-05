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


def bulmer2(a0: float, e: float, r: float, g: int) -> list[dict[str, float]]:
    """Two-trait Bulmer recursion: both traits A = a0, E = e, independent, every mate correlation r.

    Infinitesimal model, no inbreeding, C = 0. Per generation, with Sigma_A the
    2x2 A covariance and Sigma_y = Sigma_A + e*I, mates' liabilities have
    cross-covariance R_y = D R D (D the liability SDs, R = r * ones), and
    mates' A cross-covariance is C = Sigma_A Sigma_y^-1 R_y Sigma_y^-1 Sigma_A.
    Offspring: Sigma_A' = Sigma_A / 2 + (C + C') / 4 + (a0 / 2) I.
    Element t holds the moments after t matings, plus the mate moments of the
    parents of generation t (t >= 1).
    """
    sa = np.eye(2) * a0
    out = [{"var_A": a0, "var_A_t2": a0, "r_w_A": 0.0, "r_w_L": 0.0}]
    for _ in range(g):
        sy = sa + e * np.eye(2)
        d = np.sqrt(np.diag(sy))
        ry = np.outer(d, d) * r
        k = sa @ np.linalg.solve(sy, np.eye(2))
        c = k @ ry @ k.T
        sda = np.sqrt(np.diag(sa))
        mate = {
            "r_mate_A": c[0, 0] / sa[0, 0],
            "r_mate_A12": c[0, 1] / (sda[0] * sda[1]),
            "r_w_L_parents": sy[0, 1] / (d[0] * d[1]),
        }
        sa = sa / 2 + (c + c.T) / 4 + np.eye(2) * a0 / 2
        sy = sa + e * np.eye(2)
        out.append(
            {
                "var_A": sa[0, 0],
                "var_A_t2": sa[1, 1],
                "var_L": sy[0, 0],
                "r_w_A": sa[0, 1] / np.sqrt(sa[0, 0] * sa[1, 1]),
                "r_w_L": sy[0, 1] / np.sqrt(sy[0, 0] * sy[1, 1]),
                **mate,
            }
        )
    return out


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


def predict1(v0: float, s: float, r: float, g: int, c: float = 0.0) -> dict[str, float]:
    """Single-trait Bulmer expectation of every statistic after g >= 1 matings at mate correlation r.

    v0 is the founders' Var(A), s the non-additive liability variance (C + E),
    and c the household C that full sibs share (household = mother; a child's
    C is independent of its parents' liabilities, so it adds to no other row).
    """
    v = bulmer(v0, s, r, g)
    vt, vp = v[g], v[g - 1]
    lt, lp = vt + s, vp + s
    mid = vp / 2 * (1 + r * vp / lp)  # Var(midparent A) = Cov(A) of full sibs
    po_l = vp * (1 + r) / 2 / np.sqrt(lp * lt)
    return {
        "var_A": vt,
        "var_L": lt,
        "r_mate_L": r,
        "r_mate_A": r * vp / lp,
        "r_FS_A": mid / vt,
        "r_FS_L": (mid + c) / lt,
        "r_PO_A": mid / np.sqrt(vp * vt),
        "r_PO_L": po_l,
        "r_MO_L": po_l,
        "r_FO_L": po_l,
    }
