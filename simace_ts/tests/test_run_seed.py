"""simace_ts/run_seed.py: seed namespaces and the case config (ADR 0002 item 8)."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_seed


def test_stage_seeds_are_distinct_across_stages_chromosomes_and_seeds():
    seeds = range(990000, 990003)
    keys = [(s, run_seed.STAGE_SIMULATE, 0) for s in seeds]
    keys += [(s, run_seed.STAGE_DROP, c) for s in seeds for c in run_seed.CHROMOSOMES]
    keys += [(s, stage, i) for s in seeds for stage in (2, 3, 4) for i in range(6)]
    values = [run_seed.stage_seed(*k) for k in keys]
    assert len(set(values)) == len(values)
    assert all(1 <= v < 2**31 - 1 for v in values)


def test_stage_seed_is_deterministic():
    assert run_seed.stage_seed(9000, run_seed.STAGE_DROP, 22) == run_seed.stage_seed(9000, run_seed.STAGE_DROP, 22)


def test_consecutive_seeds_do_not_share_drop_seeds():
    # The old offsets s + 100 * chrom collided between reps; namespaced seeds must not.
    a = {run_seed.stage_seed(9000 + r, run_seed.STAGE_DROP, c) for r in range(150) for c in run_seed.CHROMOSOMES}
    assert len(a) == 150 * len(run_seed.CHROMOSOMES)


def test_case_config_keeps_the_adr_restrictions():
    case = yaml.safe_load(run_seed.CASE.read_text())
    assert case["simulate"]["p-mztwin"] == 0
    assert case["cohort"]["N-sample"] == 0
    assert "analyze" not in case  # analyze takes censor-age from cohort
    assert not any(k.startswith("Ctm") for k in case["simulate"])
    assert len({a["name"] for a in case["architectures"]}) == len(case["architectures"])


def test_flags():
    assert run_seed.flags({"G-ped": 10, "alpha": -0.5}) == ["--G-ped", "10", "--alpha", "-0.5"]


def test_stage_with_a_slash_in_its_name_runs_and_is_skipped_on_rerun(tmp_path):
    stages = run_seed.Stages(tmp_path)
    flag = tmp_path / "ran"
    cmd = [sys.executable, "-c", f"open({str(flag)!r}, 'a').write('x')"]

    stages.run("c1k_a05/effects", cmd)
    stages.run("c1k_a05/effects", cmd)

    assert flag.read_text() == "x"
    assert stages.marker("c1k_a05/effects").exists()
    rows = (tmp_path / "timing.tsv").read_text().splitlines()
    assert rows[0] == "stage\twall_s\tmax_rss_mb\tcgroup_peak_mb"
    assert [r.split("\t")[0] for r in rows[1:]] == ["c1k_a05/effects"]


def test_failed_stage_raises_and_leaves_no_marker(tmp_path):
    stages = run_seed.Stages(tmp_path)
    with pytest.raises(RuntimeError, match="exited 3"):
        stages.run("bad", [sys.executable, "-c", "raise SystemExit(3)"])
    assert not stages.marker("bad").exists()


def test_rerun_with_a_changed_case_is_refused(tmp_path):
    first = {"case": {"A1": 0.4}, "seed": 1, "runs": [{"simace_validate": "a"}]}
    run_seed.check_manifest(tmp_path, first)
    run_seed.check_manifest(tmp_path, {"case": {"A1": 0.4}, "seed": 1, "runs": [{"simace_validate": "b"}]})
    assert [r["simace_validate"] for r in json.loads((tmp_path / "manifest.json").read_text())["runs"]] == ["a", "b"]

    with pytest.raises(SystemExit, match="case changed"):
        run_seed.check_manifest(tmp_path, {"case": {"A1": 0.5}, "seed": 1, "runs": [{"simace_validate": "c"}]})


def test_chromosome_subset_gets_its_own_directory(tmp_path):
    case = {"id": "x"}
    full = run_seed.seed_dir(tmp_path, case, 7, run_seed.CHROMOSOMES)
    subset = run_seed.seed_dir(tmp_path, case, 7, (22,))
    assert full.name == "seed7"
    assert subset.name == "seed7_chr22"


def test_manifest_write_leaves_no_temporary_file(tmp_path):
    run_seed.check_manifest(tmp_path, {"case": {}, "runs": []})
    assert [p.name for p in tmp_path.iterdir()] == ["manifest.json"]
