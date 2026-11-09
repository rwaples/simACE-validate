"""inheritance.power: mutations of the scored parameters, cell configs and the table row."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from inheritance.power import BASE, CELLS, DELTAS, MUTATIONS, merged, parse_span, study_seed, summarize

PARAMS = {"A1": 0.5, "Ctm_paternal1": 0.3, "Ctm_maternal1": 0.5, "C1": 0.2}


def test_mutations_shift_one_trait1_parameter() -> None:
    assert len(MUTATIONS) == 1 + 3 * 2 * len(DELTAS)
    assert MUTATIONS["null"].apply(PARAMS) == PARAMS
    assert MUTATIONS["bf_minus_010"].apply(PARAMS) == PARAMS | {"Ctm_paternal1": pytest.approx(0.2)}
    assert MUTATIONS["bm_plus_050"].apply(PARAMS) == PARAMS | {"Ctm_maternal1": pytest.approx(1.0)}
    # A is off by a share: the variance cell then expects (1 + d) * A / 2.
    assert MUTATIONS["a_minus_020"].apply(PARAMS) == PARAMS | {"A1": pytest.approx(0.4)}
    assert PARAMS["A1"] == 0.5  # apply copies


def test_cells_change_base_in_one_respect() -> None:
    ctm = merged(BASE, CELLS["ctm_base"][1])
    assert ctm["pedigree"]["trait1"] == {
        "A": 0.5,
        "C": 0.2,
        "E": 0.3,
        "C_transmission": {"maternal": 0.5, "paternal": 0.3},
    }
    assert ctm["pedigree"]["rA"] == BASE["pedigree"]["rA"]
    assert "C_transmission" not in BASE["pedigree"]["trait1"]
    assert merged(BASE, CELLS["rm_n20000"][1])["N"] == 20000


def test_study_seeds_leave_room_for_the_reps() -> None:
    seeds = [study_seed(cell, j) for cell in CELLS for j in range(1, 101)]
    assert len(set(seeds)) == len(seeds)
    assert max(seeds) + 19 < 2**31


def test_parse_span() -> None:
    assert parse_span("3-5") == [3, 4, 5]
    assert parse_span("7") == [7]


def test_summarize_attributes_failures_to_the_target_test() -> None:
    outcomes = [
        {"verdict": "fail", "failed": ["mendelian_residual/trait1_variance/mean_z", "total_variance/trait1/mean_z"]},
        {"verdict": "fail", "failed": ["total_variance/trait1/sd_z"]},
        {"verdict": "pass", "failed": []},
        {"verdict": "insufficient", "failed": []},
    ]
    row = summarize("rm_base", "a_plus_005", outcomes).split(" | ")
    assert row[3] == "2/4 (1 insufficient)"
    assert row[5] == "1/4"
    assert row[6] == "1 / 0 / 0"
    assert row[7] == "total_variance (2) |"
