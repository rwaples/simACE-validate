"""Infinitesimal (Bulmer) expectations under assortative mating (issue #4, ADR 0003).

They feed the pilot summaries and the harness's INFO `bulmer:` rows only; no
reference sample is built from them, so they sit outside the reference
manifest's hashes (am/manifest.py). numpy only.
"""

from __future__ import annotations

import numpy as np


def bulmer2(a0: float, e: float, r, g: int) -> list[dict[str, float]]:
    """Two-trait Bulmer recursion: both traits A = a0, E = e, independent founders.

    r is the liability mate correlation: a float (every mate correlation r) or a
    2x2 matrix R (female trait rows, male trait columns). Infinitesimal model,
    no inbreeding, C = 0. Per generation, with Sigma_A the 2x2 A covariance and
    Sigma_y = Sigma_A + e*I, mates' liabilities have cross-covariance
    R_y = D R D (D the liability SDs), and mates' A cross-covariance is
    C = Sigma_A Sigma_y^-1 R_y Sigma_y^-1 Sigma_A.
    Offspring: Sigma_A' = Sigma_A / 2 + (C + C') / 4 + (a0 / 2) I.
    Element t holds the moments after t matings, plus the mate moments of the
    parents of generation t (t >= 1).
    """
    rmat = np.full((2, 2), float(r)) if np.ndim(r) == 0 else np.asarray(r, dtype=float)
    sa = np.eye(2) * a0
    out = [{"var_A": a0, "var_A_t2": a0, "r_w_A": 0.0, "r_w_L": 0.0}]
    for _ in range(g):
        sy = sa + e * np.eye(2)
        d = np.sqrt(np.diag(sy))
        ry = np.outer(d, d) * rmat
        k = sa @ np.linalg.solve(sy, np.eye(2))
        c = k @ ry @ k.T
        sda = np.sqrt(np.diag(sa))
        mate = {
            "r_mate_A": c[0, 0] / sa[0, 0],
            "r_mate_A_t2": c[1, 1] / sa[1, 1],
            "r_mate_A12": c[0, 1] / (sda[0] * sda[1]),
            "r_mate_A21": c[1, 0] / (sda[0] * sda[1]),
            "r_mate_L": rmat[0, 0],
            "r_mate_L_t2": rmat[1, 1],
            "r_mate_L12": rmat[0, 1],
            "r_mate_L21": rmat[1, 0],
            "r_w_L_parents": sy[0, 1] / (d[0] * d[1]),
        }
        sa = sa / 2 + (c + c.T) / 4 + np.eye(2) * a0 / 2
        sy = sa + e * np.eye(2)
        out.append(
            {
                "var_A": sa[0, 0],
                "var_A_t2": sa[1, 1],
                "var_L": sy[0, 0],
                "var_L_t2": sy[1, 1],
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
