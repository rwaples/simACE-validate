"""Layer-3 statistics of one recorded pedigree, computed the same way for simACE and SLiM.

Each statistic is one number per pedigree and generation g >= 1, so the unit
of the two-sample tests is the pedigree (ADR 0001). Pair counts come from the
parent links directly, not from pedigree-graph.
"""

from __future__ import annotations

import subprocess
from typing import TYPE_CHECKING

import numpy as np
import polars as pl

if TYPE_CHECKING:
    from pathlib import Path


def pairs_within(counts: np.ndarray) -> int:
    """Return the number of unordered pairs within groups of the given sizes."""
    return int(np.sum(counts * (counts - 1) // 2))


def structure_stats(ped: pl.DataFrame) -> dict[str, float]:
    """Return per-generation counts and offspring-number variances of one pedigree."""
    gen = ped["generation"].to_numpy()
    sex = ped["sex"].to_numpy()
    mother = ped["mother"].to_numpy()
    father = ped["father"].to_numpy()
    out: dict[str, float] = {}
    for g in range(1, int(gen.max()) + 1):
        kids = gen == g
        mo, fa = mother[kids], father[kids]
        prev = np.flatnonzero(gen == g - 1)
        out[f"frac_female_g{g}"] = float(np.mean(sex[kids] == 0))
        out[f"n_mothers_g{g}"] = len(np.unique(mo))
        out[f"n_fathers_g{g}"] = len(np.unique(fa))
        for label, parents, want in (("mother", mo, 0), ("father", fa, 1)):
            pool = prev[sex[prev] == want]
            idx = np.searchsorted(pool, parents)
            if not np.array_equal(pool[np.minimum(idx, len(pool) - 1)], parents):
                raise ValueError(f"generation {g}: a {label} is not a previous-generation individual of that sex")
            n_kids = np.bincount(idx, minlength=len(pool))
            out[f"var_offspring_per_{label}_g{g}"] = float(np.var(n_kids, ddof=1))
        _, by_couple = np.unique(np.stack([mo, fa]), axis=1, return_counts=True)
        fs = pairs_within(by_couple)
        out[f"n_FS_g{g}"] = fs
        out[f"n_MHS_g{g}"] = pairs_within(np.unique(mo, return_counts=True)[1]) - fs
        out[f"n_PHS_g{g}"] = pairs_within(np.unique(fa, return_counts=True)[1]) - fs
        if g >= 2:
            mo_mo, mo_fa = mother[mo], father[mo]
            fa_mo, fa_fa = mother[fa], father[fa]
            out[f"n_sib_mated_g{g}"] = int(np.sum((mo_mo == fa_mo) | (mo_fa == fa_fa)))
    return out


def run_slim(slim_script: Path, out: Path, N: int, G: int, seed: int) -> pl.DataFrame:
    """Run slim/wf.slim and return its pedigree with rows in id order and 0-based links."""
    cmd = ["slim", "-d", f"N={N}", "-d", f"G={G}", "-d", f"SEED={seed}", "-d", f"OUT='{out}'", str(slim_script)]
    subprocess.run(cmd, check=True, capture_output=True)
    raw = pl.read_csv(out, separator="\t").sort("generation", "id")
    index = {pid: i for i, pid in enumerate(raw["id"].to_list())}

    def remap(col: str) -> np.ndarray:
        return np.array([index.get(p, -1) for p in raw[col].to_list()], dtype=np.int64)

    return pl.DataFrame(
        {
            "generation": raw["generation"],
            "sex": raw["sex"],
            "mother": remap("mother"),
            "father": remap("father"),
        }
    )
