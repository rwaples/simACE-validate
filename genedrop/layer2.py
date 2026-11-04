"""Layer-2 statistics: one gene-drop per simACE pedigree.

Per seed, two kinds of difference with expectation 0:

- ``ibd_<class>``: mean over the class's pairs of (realised IBD fraction -
  2 phi), with phi from kinship2 (oracles/pairs.R). This tests that simACE's
  recorded links (mother, father, twin) form the pedigree their kinship
  says they do.
- ``innov_*``: a statistic of the Mendelian-sampling term X - midparent,
  computed on simACE's A1 and on the gene-drop's finite-locus G, as simACE
  minus gene-drop. Linkage makes realised sharing vary at any number of
  loci, so only the means are compared. The spread ratio is INFO
  (``spread_*``).

``innov_var_g*`` uses only children whose recorded parents both have F = 0.
There simACE's A/2 and the drop's (A/2)(1 - (F_m + F_f)/2) agree, so the
null needs no formula. The drop treats recorded generation 0 as unrelated,
so recorded F is exactly the drop's F; simACE's innovations ignore F, so
generation 0's hidden sibships cannot leak in. At F = 0 a parent is
heterozygous at a locus with probability 1/2 whatever the drift elsewhere
in the pedigree, and drop.py scales effects by the expected (not realised)
founder variance, so the drop's innovation variance there is A/2 in
expectation. Children of inbred parents are reported as INFO rows
``d1_innov_var_g*``, gene-drop evidence on D1. The sib correlations
``innov_r_*`` keep every child: independent meioses leave sibs'
innovations uncorrelated in both engines at any F.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
import polars as pl

from genedrop.drop import drop

if TYPE_CHECKING:
    from genedrop.drop import Genome

SIB_CLASSES = ("FS", "MHS", "PHS")


def innovations(x: np.ndarray, ped: pl.DataFrame) -> np.ndarray:
    """Return x - midparent for rows in generation >= 1 (NaN for generation 0)."""
    mo = ped["mother"].to_numpy()
    fa = ped["father"].to_numpy()
    out = np.full(len(x), np.nan)
    nf = mo >= 0
    out[nf] = x[nf] - (x[mo[nf]] + x[fa[nf]]) / 2
    return out


def pair_r(d: np.ndarray, i: np.ndarray, j: np.ndarray) -> float:
    """Return the double-entry pair correlation of d over pairs (i, j)."""
    a = np.concatenate([d[i], d[j]])
    b = np.concatenate([d[j], d[i]])
    return float(np.corrcoef(a, b)[0, 1])


def var_over(d: np.ndarray, mask: np.ndarray, A: float) -> float:
    """Return the variance of d over mask in units of A; NaN (absent) with fewer than two rows."""
    return float(np.var(d[mask], ddof=1) / A) if mask.sum() >= 2 else np.nan


def innovation_stats(
    x: np.ndarray, ped: pl.DataFrame, pairs: pl.DataFrame, F: np.ndarray, A: float
) -> dict[str, float]:
    """Return x's innovation variance (in units of A) by generation and parental F, and its sib correlations."""
    d = innovations(x, ped)
    gen = ped["generation"].to_numpy()
    twin = ped["twin"].to_numpy()
    own = (twin < 0) | (twin > np.arange(len(x)))  # one row per MZ pair
    mo = ped["mother"].to_numpy()
    fa = ped["father"].to_numpy()
    outbred = np.zeros(len(x), dtype=bool)
    nf = mo >= 0
    outbred[nf] = (F[mo[nf]] == 0) & (F[fa[nf]] == 0)  # kinship2's phi is an exact dyadic sum
    out = {}
    for g in range(1, int(gen.max()) + 1):
        rows = (gen == g) & own
        out[f"innov_var_g{g}"] = var_over(d, rows & outbred, A)
        out[f"d1_innov_var_g{g}"] = var_over(d, rows & ~outbred, A)
    for cl in SIB_CLASSES:
        p = pairs.filter(pl.col("class") == cl)
        out[f"innov_r_{cl}"] = pair_r(d, p["i"].to_numpy(), p["j"].to_numpy())
    return out


def seed_stats(
    ped: pl.DataFrame, pairs: pl.DataFrame, F: np.ndarray, genome: Genome, A: float, seed: int
) -> dict[str, float]:
    """Gene-drop once through one pedigree and return each statistic's difference; F is recorded inbreeding."""
    G, ibd = drop(ped, pairs.select("i", "j").to_numpy(), genome, A, seed)
    pairs = pairs.with_columns(gap=ibd - pairs["phi2"].to_numpy())
    out = {f"ibd_{cl}": float(g["gap"].mean()) for (cl,), g in pairs.group_by("class", maintain_order=True)}
    sim = innovation_stats(ped["A1"].to_numpy().astype(np.float64), ped, pairs, F, A)
    gd = innovation_stats(G, ped, pairs, F, A)
    out |= {k: sim[k] - gd[k] for k in sim}
    out |= {f"spread_{k}": (sim[k], gd[k]) for k in sim}
    return out
