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


def innovation_stats(x: np.ndarray, ped: pl.DataFrame, pairs: pl.DataFrame, A: float) -> dict[str, float]:
    """Return the variance (in units of A) by generation and sib correlations of x's innovations."""
    d = innovations(x, ped)
    gen = ped["generation"].to_numpy()
    twin = ped["twin"].to_numpy()
    own = (twin < 0) | (twin > np.arange(len(x)))  # one row per MZ pair
    out = {f"innov_var_g{g}": float(np.var(d[(gen == g) & own], ddof=1) / A) for g in range(1, int(gen.max()) + 1)}
    for cl in SIB_CLASSES:
        p = pairs.filter(pl.col("class") == cl)
        out[f"innov_r_{cl}"] = pair_r(d, p["i"].to_numpy(), p["j"].to_numpy())
    return out


def seed_stats(ped: pl.DataFrame, pairs: pl.DataFrame, genome: Genome, A: float, seed: int) -> dict[str, float]:
    """Gene-drop once through one pedigree and return each statistic's difference."""
    G, ibd = drop(ped, pairs.select("i", "j").to_numpy(), genome, A, seed)
    pairs = pairs.with_columns(gap=ibd - pairs["phi2"].to_numpy())
    out = {f"ibd_{cl}": float(g["gap"].mean()) for (cl,), g in pairs.group_by("class", maintain_order=True)}
    sim = innovation_stats(ped["A1"].to_numpy().astype(np.float64), ped, pairs, A)
    gd = innovation_stats(G, ped, pairs, A)
    out |= {k: sim[k] - gd[k] for k in sim}
    out |= {f"spread_{k}": (sim[k], gd[k]) for k in sim}
    return out
