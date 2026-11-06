"""Compute one simACE-ts seed's Part from its outputs (ADR 0002 items 2 to 5 and 9).

Per architecture, from that architecture's augmented pedigree, cohort and
``report.yaml``:

- ``scored``: simACE's value and its realised truth, for fields whose
  expectation given the supplied A follows from the configured C and E law.
  A pair correlation over pairs (i, j) has truth
  ``(cov(A_i, A_j) + C1 * share) / sqrt((var A_i + C1 + E1) (var A_j + C1 + E1))``
  and a midparent slope ``(cov(A_child, A_mid) + C1 * share_cm) / var_mid`` with
  ``var_mid = var A_mid + (C1 (1 + share_mf) + E1) / 2``. The moments of A are
  population moments of the realised A over exactly simACE's pairs; each
  ``share`` is the realised fraction of pairs in one household (the
  pedigree's ``household_id``). The C and E variances are the configured
  ones, not the realised, so a wrong C or E variance in simACE shows.
- ``exact``: fields that are functions of the supplied values alone,
  recomputed here with simACE's arithmetic; they must match to ``tol``.
- ``info``: reported, never scored (Gaussian-liability quantities and the
  pipeline diagnostic, by class and generation).

A value simACE leaves empty (too few pairs) is recorded as ``None``, never
dropped, so the harness can make the Case ERROR (ADR 0002 item 9).

Pairs are enumerated here from the mother and father links, each once as
(lower id, higher id), as simACE's moments table orders them. Only classes
defined by direct links are used (FS, MHS, PHS, MO, FO).

Usage::

    python simace_ts/part.py <seed directory>     # writes part.json; exits 1 on an exact mismatch
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import polars as pl
import yaml
from run_seed import validate_commit
from score import PAIR_CLASSES, PART_VERSION, po_generations

TWO_PHI = {"FS": 0.5, "MHS": 0.25, "PHS": 0.25, "MO": 0.5, "FO": 0.5}  # outbred 2 * kinship, for the INFO diagnostic
TOL = 1e-10  # simACE's own arithmetic, so only summation order differs: up to 4.7e-12 on a 360k-row smoke seed


def class_pairs(ped: pl.DataFrame, classes: tuple[str, ...] = PAIR_CLASSES) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Return each requested class's (i, j) row arrays, i < j, among the rows of ``ped`` (ids equal row positions)."""
    kids = ped.filter(pl.col("mother") >= 0).select("id", "mother", "father")
    ids = ped.select("id")
    out = {}
    if {"FS", "MHS"} & set(classes):
        by_mother = kids.join(kids, on="mother", suffix="_r").filter(pl.col("id") < pl.col("id_r"))
        same_father = pl.col("father") == pl.col("father_r")
        out["FS"] = by_mother.filter(same_father)
        out["MHS"] = by_mother.filter(~same_father)
    if "PHS" in classes:
        out["PHS"] = kids.join(kids, on="father", suffix="_r").filter(
            (pl.col("id") < pl.col("id_r")) & (pl.col("mother") != pl.col("mother_r"))
        )
    if "MO" in classes:
        out["MO"] = kids.join(ids, left_on="mother", right_on="id", how="semi").select(id="mother", id_r="id")
    if "FO" in classes:
        out["FO"] = kids.join(ids, left_on="father", right_on="id", how="semi").select(id="father", id_r="id")
    return {k: (out[k]["id"].to_numpy(), out[k]["id_r"].to_numpy()) for k in classes}


def cov0(x: np.ndarray, y: np.ndarray) -> float:
    """Population covariance."""
    return float(np.mean((x - x.mean()) * (y - y.mean())))


def pearson(x: np.ndarray, y: np.ndarray) -> float:
    """Pearson correlation; NaN when either side is constant."""
    denom = math.sqrt(float(np.var(x)) * float(np.var(y)))
    return cov0(x, y) / denom if denom > 0 else math.nan


def truth_r(a: np.ndarray, i: np.ndarray, j: np.ndarray, share: float, c1: float, e1: float) -> float:
    """Expected liability correlation over pairs (i, j) given the realised A, with ``share`` of pairs in one household."""
    return (cov0(a[i], a[j]) + share * c1) / math.sqrt((np.var(a[i]) + c1 + e1) * (np.var(a[j]) + c1 + e1))


def midparent(v: np.ndarray, child: np.ndarray, mother: np.ndarray, father: np.ndarray) -> tuple[float, float]:
    """Return the OLS slope and correlation of v[child] on the midparent of v."""
    x = (v[mother] + v[father]) / 2
    y = v[child]
    var_x = float(np.var(x))
    return (cov0(x, y) / var_x if var_x > 0 else math.nan), pearson(x, y)


def midparent_truth(
    a: np.ndarray, hh: np.ndarray, child: np.ndarray, mother: np.ndarray, father: np.ndarray, c1: float, e1: float
) -> tuple[float, float]:
    """Return the expected midparent slope and correlation of liability given the realised A and households."""
    x = (a[mother] + a[father]) / 2
    y = a[child]
    share_mf = float(np.mean(hh[mother] == hh[father]))  # parents from one household share C
    share_cm = (float(np.mean(hh[child] == hh[mother])) + float(np.mean(hh[child] == hh[father]))) / 2
    cov = cov0(x, y) + c1 * share_cm
    var_mid = float(np.var(x)) + (c1 * (1 + share_mf) + e1) / 2
    return cov / var_mid, cov / math.sqrt(var_mid * (float(np.var(y)) + c1 + e1))


def num(value: object) -> float | None:
    """Return a finite float, or None for a value simACE left empty or that is not finite."""
    if value is None:
        return None
    f = float(value)
    return f if math.isfinite(f) else None


def architecture_part(d: Path, c1: float, e1: float, gens: list[int]) -> dict:
    """Return the scored, exact and INFO statistics of one architecture's outputs in ``d``."""
    ped = pl.read_parquet(d / "pedigree.full.tstrait.parquet")
    if not np.array_equal(ped["id"].to_numpy(), np.arange(ped.height)):
        raise ValueError(f"{d}: pedigree ids must equal row order")
    cohort = pl.read_parquet(d / "cohort.parquet")
    report = yaml.safe_load((d / "report.yaml").read_text())
    obs = report["observed"]["analysis_sample"]
    est = report["estimators"]["heritability"]["liability_scale"]["trait1"]

    a = ped["A1"].to_numpy().astype(np.float64)
    # simACE reads the stored liability for the analysis sample and sums the float32 parts for the recorded pedigree.
    liab = ped["liability1"].to_numpy().astype(np.float64)
    liab_parts = (ped["A1"].to_numpy() + ped["C1"].to_numpy() + ped["E1"].to_numpy()).astype(np.float64)
    hh = ped["household_id"].to_numpy()
    gen = ped["generation"].to_numpy()
    in_sample = np.zeros(ped.height, dtype=bool)
    in_sample[cohort.filter(pl.col("affected1").is_not_null())["id"].to_numpy()] = True

    scored: dict[str, dict[str, float | None]] = {}
    exact: dict[str, dict[str, float | None]] = {}
    info: dict[str, float | None] = {}

    def add_exact(name: str, simace: object, harness: float, tol: float = TOL) -> None:
        exact[name] = {"simace": num(simace), "harness": num(harness), "tol": tol}

    def add_scored(name: str, simace: object, truth: float) -> None:
        scored[name] = {"simace": num(simace), "truth": num(truth)}

    # Analysis sample: liability correlations by class.
    for k, (i, j) in class_pairs(ped.filter(pl.Series(in_sample))).items():
        add_exact(f"n_pairs_{k}", obs["relationship_pair_counts"][k], len(i), tol=0)
        simace_r = obs["liability_correlations"]["trait1"][k]
        add_exact(f"lc_{k}_value", simace_r, pearson(liab[i], liab[j]))
        add_scored(f"lc_{k}", simace_r, truth_r(a, i, j, float(np.mean(hh[i] == hh[j])), c1, e1))
        info[f"tetrachoric_{k}"] = num(obs["tetrachoric"]["trait1"][k]["r"])
        info[f"affected_corr_{k}"] = num(obs["affected_correlations"]["trait1"][k])
        info[f"gd_rA_{k}_minus_2phi"] = num(pearson(a[i], a[j]) - TWO_PHI[k])
        for g in np.unique(gen[j]).tolist():  # j is the younger member (the child for PO)
            m = gen[j] == g
            info[f"gd_rA_{k}_g{g}_minus_2phi"] = num(pearson(a[i[m]], a[j[m]]) - TWO_PHI[k])

    # Analysis sample: midparent regression per child generation.
    kids = ped.filter(pl.col("mother") >= 0)
    child, mother, father = (kids[c].to_numpy() for c in ("id", "mother", "father"))
    for g in gens:
        k = (gen[child] == g) & in_sample[child] & in_sample[mother] & in_sample[father]
        entry = obs["parent_offspring_corr"]["trait1"][f"gen{g}"]
        add_exact(f"po_g{g}_n", entry["n_pairs"], int(k.sum()), tol=0)
        slope, r = midparent(liab, child[k], mother[k], father[k])
        add_exact(f"po_g{g}_slope_value", entry["slope"], slope)
        add_exact(f"po_g{g}_r_value", entry["r"], r)
        t_slope, t_r = midparent_truth(a, hh, child[k], mother[k], father[k], c1, e1)
        add_scored(f"po_g{g}_slope", entry["slope"], t_slope)
        add_scored(f"po_g{g}_r", entry["r"], t_r)

    # Recorded pedigree: full-sib correlations and the pooled midparent regression.
    i, j = class_pairs(ped, ("FS",))["FS"]
    add_exact("dz_A_corr", est["dz_sibling_A_corr"], pearson(a[i], a[j]))
    add_exact("dz_L_corr_value", est["dz_sibling_liability_corr"], pearson(liab_parts[i], liab_parts[j]))
    add_scored("dz_L", est["dz_sibling_liability_corr"], truth_r(a, i, j, float(np.mean(hh[i] == hh[j])), c1, e1))
    a_slope, a_r = midparent(a, child, mother, father)
    add_exact("po_A_slope", est["parent_offspring_A_slope"], a_slope)
    add_exact("po_A_r2", est["parent_offspring_A_r2"], a_r**2)
    add_exact(
        "po_L_slope_value", est["parent_offspring_liability_slope"], midparent(liab_parts, child, mother, father)[0]
    )
    add_scored(
        "po_L_slope", est["parent_offspring_liability_slope"], midparent_truth(a, hh, child, mother, father, c1, e1)[0]
    )
    # simACE's realised var_A is generation 0's, with ddof=1 (simace quality/statistical.py, _check_variance).
    var_a = report["truth"]["recorded_pedigree"]["traits"]["trait1"]["realized"]["var_A"]
    add_exact("var_A_g0", var_a, np.var(a[gen == 0], ddof=1))

    # INFO: the Gaussian-liability quantities (D3) and the pre-rescale gene-drop variance.
    info["prevalence"] = num(obs["prevalence"]["trait1"])
    for name, value in report["estimators"]["heritability"]["observed_scale"]["trait1"].items():
        info[f"h2_observed_{name}"] = num(value)
    raw = json.loads((d / "augment_meta.json").read_text())["rescale_by_generation"]
    for row in raw:
        info[f"gd_raw_var_g{row['generation']}_over_g0"] = row["realized_var"] / raw[0]["realized_var"]
    return {"scored": scored, "exact": exact, "info": info}


def mismatches(part: dict) -> list[str]:
    """Return ``architecture/statistic`` for every exact row that does not match, or has a missing side."""
    return [
        f"{arch}/{name}"
        for arch, p in part["architectures"].items()
        for name, e in p["exact"].items()
        if e["simace"] is None or e["harness"] is None or abs(e["simace"] - e["harness"]) > e["tol"]
    ]


def seed_part(out: Path) -> dict:
    """Return the Part of the seed whose outputs are in ``out``."""
    manifest = json.loads((out / "manifest.json").read_text())
    case = manifest["case"]
    c1, e1 = float(case["simulate"]["C1"]), float(case["simulate"]["E1"])
    gens = po_generations(case)
    return {
        "part_version": PART_VERSION,
        "manifest": manifest,
        "part_by": validate_commit(),
        "architectures": {a["name"]: architecture_part(out / a["name"], c1, e1, gens) for a in case["architectures"]},
    }


def main(argv: list[str] | None = None) -> int:
    """Write ``part.json`` into a seed directory; return 1 when an exact row does not match (ADR 0002 item 3)."""
    out = Path((argv or sys.argv[1:])[0])
    part = seed_part(out)
    tmp = out / "part.json.tmp"
    tmp.write_text(json.dumps(part, indent=1) + "\n")
    tmp.replace(out / "part.json")
    print(out / "part.json")
    bad = mismatches(part)
    if bad:
        print("exact recomputation mismatches: " + ", ".join(bad), file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
