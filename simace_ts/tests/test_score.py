"""simace_ts/score.py and its harness runner: Parts that are missing, stale or empty make the Case ERROR."""

from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
import run_harness
from simace_ts import score

ROOT = Path(__file__).resolve().parents[2]
SIMACE = {"version": "0.2.0", "direct_url": {"url": "https://github.com/rwaples/simACE.git"}}
DATA = {"trees_zip_sha256": "aa", "simhumanity_commit": "bb", "preprocessed_fingerprint": "cc"}
CASE = {
    "id": "ts_test",
    "layer": 4,
    "oracle": "simace_ts",
    "simulate": {"A1": 0.4, "G-ped": 5},
    "cohort": {"G-pheno": 3},
    "architectures": [{"name": "first"}, {"name": "second"}],
    "data": DATA,
    "seed0": 100,
    "R": 20,  # the harness's MIN_SEEDS, below which a row is INFO
}
SCHEMA = score.scored_schema(CASE)
GAPS = np.random.default_rng(5).normal(0, 0.01, size=20)


def _part(seed: int, gap: float) -> dict:
    arch = {
        "scored": {stat: {"simace": 0.3 + gap, "truth": 0.3} for stat in SCHEMA},
        "exact": {"n_pairs_FS": {"simace": 10.0, "harness": 10.0, "tol": 0}},
        "info": {"prevalence": 0.1, "empty": None},
    }
    manifest = {
        "case": CASE,
        "seed": seed,
        "chromosomes": list(score.CHROMOSOMES),
        "simace": SIMACE,
        "data": DATA,
        "runs": [{"simace_validate": "37166f987aa7"}],
    }
    return {
        "part_version": score.PART_VERSION,
        "part_by": "37166f987aa7",
        "manifest": manifest,
        "architectures": {"first": copy.deepcopy(arch), "second": copy.deepcopy(arch)},
    }


def _write(parts_dir: Path, parts: list[dict], seed: int | None = None) -> None:
    d = parts_dir / CASE["id"]
    d.mkdir(parents=True, exist_ok=True)
    for p in parts:
        (d / f"seed{p['manifest']['seed'] if seed is None else seed}.json").write_text(json.dumps(p))


@pytest.fixture
def parts_dir(tmp_path: Path) -> Path:
    _write(tmp_path, [_part(100 + k, float(g)) for k, g in enumerate(GAPS)])
    return tmp_path


@pytest.fixture
def harness_parts(parts_dir, monkeypatch) -> Path:
    monkeypatch.setattr(run_harness, "load_parts", lambda case: score.load_parts(case, parts_dir, SIMACE))
    return parts_dir


def _rows(parts_dir: Path) -> score.CaseRows:
    return score.case_rows(CASE, score.load_parts(CASE, parts_dir, SIMACE))


def test_schema_is_fixed_by_the_case():
    assert [
        *(f"lc_{k}" for k in score.PAIR_CLASSES),
        *("po_g3_slope", "po_g3_r", "po_g4_slope", "po_g4_r"),
        *("dz_L", "po_L_slope"),
    ] == SCHEMA


def test_valid_parts_give_gaps_per_seed(parts_dir):
    rows = _rows(parts_dir)
    np.testing.assert_allclose(rows.scored["lc_FS"], GAPS, atol=1e-12)
    assert set(rows.info_gaps) == {f"second/{s}" for s in SCHEMA}
    assert rows.exact["n_pairs_FS"] == (0, 0.0)
    assert rows.exact_sweep["n_pairs_FS"] == (0, 0.0)
    assert np.isnan(rows.info_values["empty"]).all()


def test_missing_seed_is_refused(parts_dir):
    (parts_dir / CASE["id"] / "seed101.json").unlink()
    with pytest.raises(FileNotFoundError, match="1 of 20 Parts missing"):
        score.load_parts(CASE, parts_dir, SIMACE)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("simace", {"version": "0.3.0", "direct_url": None}, "simace 0.3.0"),
        ("data", {**DATA, "preprocessed_fingerprint": "zz"}, "data pins"),
        ("chromosomes", [22], "chromosomes"),
        ("seed", 999, "seed 999"),
        ("runs", [{"simace_validate": "37166f987aa7-dirty"}], "uncommitted"),
    ],
)
def test_stale_manifest_is_refused(parts_dir, field, value, message):
    p = _part(101, 0.0)
    p["manifest"][field] = value
    _write(parts_dir, [p], seed=101)
    with pytest.raises(ValueError, match=message):
        score.load_parts(CASE, parts_dir, SIMACE)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [("part_version", 0, "part_version 0"), ("part_by", "abc-dirty", "uncommitted")],
)
def test_part_from_other_or_dirty_code_is_refused(parts_dir, field, value, message):
    p = _part(101, 0.0)
    p[field] = value
    _write(parts_dir, [p])
    with pytest.raises(ValueError, match=message):
        score.load_parts(CASE, parts_dir, SIMACE)


def test_changed_output_key_is_refused_but_r_and_note_are_not(parts_dir):
    with pytest.raises(ValueError, match="simulate"):
        score.load_parts(CASE | {"simulate": {"A1": 0.5, "G-ped": 5}}, parts_dir, SIMACE)
    assert len(score.load_parts(CASE | {"note": "edited", "R": 19}, parts_dir, SIMACE)) == 19


def test_part_missing_a_schema_row_is_refused(parts_dir):
    p = _part(105, 0.0)
    del p["architectures"]["first"]["scored"]["po_g4_r"]
    _write(parts_dir, [p])
    with pytest.raises(ValueError, match="not the schema"):
        _rows(parts_dir)


def test_part_with_different_exact_rows_is_refused(parts_dir):
    p = _part(105, 0.0)
    p["architectures"]["second"]["exact"]["extra"] = {"simace": 1.0, "harness": 1.0, "tol": 0}
    _write(parts_dir, [p])
    with pytest.raises(ValueError, match="exact statistics differ"):
        _rows(parts_dir)


def test_empty_scored_value_is_refused(parts_dir):
    p = _part(102, 0.0)
    p["architectures"]["first"]["scored"]["lc_FS"]["simace"] = None
    _write(parts_dir, [p])
    with pytest.raises(ValueError, match="first/lc_FS is empty"):
        _rows(parts_dir)


def test_missing_info_value_is_nan_not_an_error(parts_dir):
    p = _part(103, 0.0)
    del p["architectures"]["first"]["info"]["prevalence"]
    _write(parts_dir, [p])
    v = _rows(parts_dir).info_values["prevalence"]
    assert np.isnan(v[3])
    assert np.isfinite(np.delete(v, 3)).all()


def test_exact_mismatch_counts_scored_and_sweep_separately(parts_dir):
    p = _part(102, 0.0)
    p["architectures"]["second"]["exact"]["n_pairs_FS"]["harness"] = 11.0
    q = _part(101, 0.0)
    q["architectures"]["first"]["exact"]["n_pairs_FS"]["simace"] = None
    _write(parts_dir, [p, q])
    rows = _rows(parts_dir)
    assert rows.exact["n_pairs_FS"] == (1, 0.0)
    assert rows.exact_sweep["n_pairs_FS"] == (1, 1.0)


def test_harness_runner_rows_and_status(harness_parts):
    rows = run_harness.run_simace_ts(CASE)
    by_stat = {r["stat"]: r for r in rows}
    threshold = run_harness.assign_status(rows)
    assert threshold == pytest.approx(0.05 / len(SCHEMA))
    assert by_stat["lc_FS"]["status"] == "PASS"
    assert by_stat["exact/n_pairs_FS"]["status"] == "PASS"
    assert by_stat["exact_sweep/n_pairs_FS"]["status"] == "INFO"
    assert by_stat["second/lc_FS"]["status"] == "INFO"
    assert by_stat["info/prevalence"]["status"] == "INFO"
    md = run_harness.minimum_detectable(rows, threshold)
    assert any(line.startswith("| ts_test | lc_FS | 20 |") for line in md)
    assert not any("second/" in line for line in md)


def test_exact_mismatch_in_scored_architecture_fails_outside_k(harness_parts):
    p = _part(100, float(GAPS[0]))
    p["architectures"]["first"]["exact"]["n_pairs_FS"]["harness"] = 9.0
    _write(harness_parts, [p])
    rows = run_harness.run_simace_ts(CASE)
    assert run_harness.assign_status(rows) == pytest.approx(0.05 / len(SCHEMA))
    exact = {r["stat"]: r for r in rows}["exact/n_pairs_FS"]
    assert exact["status"] == "FAIL"
    assert exact["score"].effect == 1


def test_exact_mismatch_in_sweep_is_info_only(harness_parts):
    p = _part(100, float(GAPS[0]))
    p["architectures"]["second"]["exact"]["n_pairs_FS"]["harness"] = 9.0
    _write(harness_parts, [p])
    rows = run_harness.run_simace_ts(CASE)
    run_harness.assign_status(rows)
    status = {r["stat"]: r["status"] for r in rows}
    assert status["exact/n_pairs_FS"] == "PASS"
    assert status["exact_sweep/n_pairs_FS"] == "INFO"


def test_run_case_does_not_simulate_a_simace_ts_case(harness_parts, monkeypatch):
    monkeypatch.setattr(run_harness, "simulate", lambda *a: pytest.fail("simace_ts must not simulate"))
    assert run_harness.run_case(CASE)


def test_replicate_run_skips_layer_4(monkeypatch, capsys):
    seen = []
    monkeypatch.setattr(sys, "argv", ["run_harness.py", "--only", "simace_ts_p2,rm_base", "--replicate", "1"])
    monkeypatch.setattr(run_harness, "run_case", lambda case: seen.append(case["id"]) or [])
    monkeypatch.setattr(run_harness, "write_results", lambda *a: None)
    assert run_harness.main() == 0
    assert seen == ["rm_base~1"]
    assert "--replicate skips simace_ts_p2" in capsys.readouterr().out


def test_cases_yaml_case_has_pins_and_a_fixed_r():
    case = {c["id"]: c for c in run_harness.yaml.safe_load((ROOT / "cases.yaml").read_text())["cases"]}["simace_ts_p2"]
    assert set(case["data"]) == {"trees_zip_sha256", "simhumanity_commit", "preprocessed_fingerprint"}
    assert 0 < case["R"] < 100
