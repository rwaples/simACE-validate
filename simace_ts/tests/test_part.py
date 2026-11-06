"""simace_ts/part.py against a small real simACE seed: every exact row must match simACE's report."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import polars as pl
import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import part
import run_seed

CASE = {
    "id": "tiny",
    "simulate": {
        "N": 400, "G-ped": 5, "G-sim": 5, "mating-model": "standard", "mating-lambda": 0.5, "p-mztwin": 0,
        "A1": 0.4, "C1": 0.1, "E1": 0.5, "A2": 0.4, "C2": 0.1, "E2": 0.5, "rA": 0, "rC": 0, "rE": 0,
    },
    "cohort": {
        "G-pheno": 3, "phenotype-model1": "frailty", "phenotype-params1": "{distribution: weibull, scale: 2160, rho: 0.8}",
        "phenotype-model2": "frailty", "phenotype-params2": "{distribution: weibull, scale: 2160, rho: 0.8}",
        "N-sample": 0, "censor-age": 80,
    },
    "architectures": [{"name": "gaussian"}],
}  # fmt: skip


def _simace(*args: str) -> None:
    subprocess.run([sys.executable, "-m", "simace", *args], check=True, capture_output=True)


@pytest.fixture(scope="module")
def seed_dir(tmp_path_factory) -> Path:
    """A seed directory laid out as run_seed.py leaves it, with simACE's own A1 in place of a gene-drop A."""
    out = tmp_path_factory.mktemp("seed")
    d = out / "gaussian"
    d.mkdir()
    ped = d / "pedigree.full.tstrait.parquet"
    _simace("simulate", "--seed", "11", *run_seed.flags(CASE["simulate"]), "--output-pedigree", str(ped))
    params = run_seed.write_params(CASE, 11, out / "params.yaml")
    cohort, phenotyped = d / "cohort.parquet", d / "phenotyped_population.yaml"
    _simace(
        "cohort", "--pedigree", str(ped), "--output-cohort", str(cohort), "--output-phenotyped-population",
        str(phenotyped), "--seed", "12", *run_seed.flags(CASE["cohort"]),
    )  # fmt: skip
    _simace(
        "analyze", "--pedigree", str(ped), "--params", str(params), "--cohort", str(cohort),
        "--phenotyped-population", str(phenotyped), "--report-output", str(d / "report.yaml"),
        "--plot-payload-output", str(d / "plot_payload.yaml"), "--samples-output", str(d / "samples.parquet"),
        "--seed", "13", "--censor-age", "80",
    )  # fmt: skip
    gens = pl.read_parquet(ped)["generation"].unique().sort().to_list()
    meta = {"rescale_by_generation": [{"generation": g, "realized_var": 1.0} for g in gens]}
    (d / "augment_meta.json").write_text(json.dumps(meta))
    (out / "manifest.json").write_text(json.dumps({"case": CASE}))
    return out


@pytest.fixture(scope="module")
def tiny_part(seed_dir) -> dict:
    return part.seed_part(seed_dir)["architectures"]["gaussian"]


def test_every_exact_row_matches_simace(tiny_part):
    bad = {k: e for k, e in tiny_part["exact"].items() if not abs(e["simace"] - e["harness"]) <= e["tol"]}
    assert not bad
    assert {f"n_pairs_{k}" for k in part.PAIR_CLASSES} <= set(tiny_part["exact"])


def test_scored_rows_are_finite_and_near_their_truth(tiny_part):
    expected = {f"lc_{k}" for k in part.PAIR_CLASSES} | {"dz_L", "po_L_slope", "po_g3_slope", "po_g3_r"}
    expected |= {"po_g4_slope", "po_g4_r"}
    assert set(tiny_part["scored"]) == expected
    for name, row in tiny_part["scored"].items():
        assert np.isfinite(row["simace"]), name
        assert np.isfinite(row["truth"]), name
        assert abs(row["simace"] - row["truth"]) < 0.25, name  # N = 400: only a gross error lands this far


def test_truth_of_pairs_identical_within_each_pair():
    a = np.array([1.0, 1.0, -1.0, -1.0, 0.0, 0.0])
    i, j = np.array([0, 2, 4]), np.array([1, 3, 5])
    # cov(A_i, A_j) = var A = 2/3 on both sides, so r = (2/3 + share * 0.1) / (2/3 + 0.6).
    assert part.truth_r(a, i, j, 1.0, 0.1, 0.5) == pytest.approx((2 / 3 + 0.1) / (2 / 3 + 0.6))
    assert part.truth_r(a, i, j, 0.0, 0.1, 0.5) == pytest.approx((2 / 3) / (2 / 3 + 0.6))


def test_midparent_truth_counts_parents_from_one_household():
    rng = np.random.default_rng(3)
    n = 2000
    a = rng.normal(size=3 * n)
    child, mother, father = np.arange(n), np.arange(n, 2 * n), np.arange(2 * n, 3 * n)
    hh = np.arange(3 * n)
    apart = part.midparent_truth(a, hh, child, mother, father, 0.1, 0.5)[0]
    hh[father] = hh[mother]  # every mate pair from one household: Var((C_m + C_f) / 2) = C1, not C1 / 2
    together = part.midparent_truth(a, hh, child, mother, father, 0.1, 0.5)[0]
    x = (a[mother] + a[father]) / 2
    cov = part.cov0(x, a[child])
    assert apart == pytest.approx(cov / (np.var(x) + (0.1 + 0.5) / 2))
    assert together == pytest.approx(cov / (np.var(x) + (2 * 0.1 + 0.5) / 2))


def test_empty_simace_value_is_none_and_a_mismatch():
    p = {"architectures": {"x": {"exact": {"a": {"simace": None, "harness": 1.0, "tol": 0}}}}}
    assert part.num(None) is None
    assert part.num(float("nan")) is None
    assert part.mismatches(p) == ["x/a"]


def test_class_pairs_on_a_hand_pedigree():
    #   0 x 1 -> 4, 5 (FS); 0 x 2 -> 6 (MHS of 4, 5); 3 x 1 -> 7 (PHS of 4, 5)
    ped = pl.DataFrame({"id": range(8), "mother": [-1, -1, -1, -1, 0, 0, 0, 3], "father": [-1, -1, -1, -1, 1, 1, 2, 1]})
    got = {k: sorted(zip(i.tolist(), j.tolist(), strict=True)) for k, (i, j) in part.class_pairs(ped).items()}
    assert got["FS"] == [(4, 5)]
    assert got["MHS"] == [(4, 6), (5, 6)]
    assert got["PHS"] == [(4, 7), (5, 7)]
    assert got["MO"] == [(0, 4), (0, 5), (0, 6), (3, 7)]
    assert got["FO"] == [(1, 4), (1, 5), (1, 7), (2, 6)]


def test_main_writes_part_json(seed_dir, capsys):
    assert part.main([str(seed_dir)]) == 0
    assert json.loads((seed_dir / "part.json").read_text())["architectures"]["gaussian"]["scored"]
    assert not list(seed_dir.glob("*.tmp"))
    assert "mismatch" not in capsys.readouterr().err


def test_main_exits_1_on_an_exact_mismatch(seed_dir, tmp_path):
    copy = tmp_path / "seed"
    shutil.copytree(seed_dir, copy)
    report = copy / "gaussian" / "report.yaml"
    data = yaml.safe_load(report.read_text())
    data["observed"]["analysis_sample"]["relationship_pair_counts"]["FS"] += 1
    report.write_text(yaml.safe_dump(data))
    assert part.main([str(copy)]) == 1


def test_missing_report_fails(seed_dir, tmp_path):
    copy = tmp_path / "seed"
    shutil.copytree(seed_dir, copy)
    (copy / "gaussian" / "report.yaml").unlink()
    with pytest.raises(FileNotFoundError):
        part.seed_part(copy)
