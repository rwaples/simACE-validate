"""inheritance.power: mutations of the scored parameters, cell configs and the table row."""

from __future__ import annotations

import sys
from itertools import pairwise
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from inheritance.power import (
    BASE,
    CELLS,
    DELTAS,
    MUTATIONS,
    REPS_PER_VERDICT,
    cell_mutations,
    in_domain,
    merged,
    parse_span,
    study_seed,
    summarize,
)

PARAMS = {"A1": 0.5, "Ctm_paternal1": 0.3, "Ctm_maternal1": 0.5, "C1": 0.2}


def test_mutations_shift_one_trait1_parameter() -> None:
    assert len(MUTATIONS) == 1 + 3 * 2 * len(DELTAS)
    assert MUTATIONS["null"].apply(PARAMS) == PARAMS
    assert MUTATIONS["bf_minus_010"].apply(PARAMS) == PARAMS | {"Ctm_paternal1": pytest.approx(0.2)}
    assert MUTATIONS["bm_plus_050"].apply(PARAMS) == PARAMS | {"Ctm_maternal1": pytest.approx(1.0)}
    # A is off by a share: the variance cell then expects (1 + d) * A / 2.
    assert MUTATIONS["a_minus_020"].apply(PARAMS) == PARAMS | {"A1": pytest.approx(0.4)}
    assert PARAMS["A1"] == 0.5  # apply copies


def test_out_of_domain_transmission_is_skipped() -> None:
    # No true transmission: b minus d is negative.
    assert not in_domain("rm_base", "bf_minus_005")
    assert in_domain("rm_base", "bm_plus_050")
    # ctm_base's trait 1 is b_m 0.5, b_f 0.3: b_m 1.0 breaks b_m² + b_f² ≤ 1, b_m 0.8 does not.
    assert not in_domain("ctm_base", "bm_plus_050")
    assert in_domain("ctm_base", "bm_plus_030")
    assert in_domain("ctm_base", "bf_minus_030")
    assert all(in_domain(cell, n) for cell in CELLS for n in MUTATIONS if n == "null" or n.startswith("a_"))
    assert cell_mutations("rm_base")[0] == "null"
    assert len(cell_mutations("rm_base")) == 1 + 2 * len(DELTAS) + 2 * len(DELTAS)


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


def test_study_replicates_never_share_a_rep_seed() -> None:
    starts = sorted(study_seed(cell, j) for cell in CELLS for j in range(1, 101))
    assert all(b - a >= REPS_PER_VERDICT for a, b in pairwise(starts))
    assert starts[-1] + REPS_PER_VERDICT - 1 < 2**31


def test_parse_span() -> None:
    assert parse_span("3-5") == [3, 4, 5]
    assert parse_span("7") == [7]
    with pytest.raises(ValueError, match="not a span"):
        parse_span("5-3")


def test_summarize_attributes_failures_to_the_target_test() -> None:
    outcomes = [
        {
            "verdict": "fail",
            "n_tests": 30,
            "failed": ["mendelian_residual/trait1_variance/mean_z", "total_variance/trait1/mean_z"],
        },
        {"verdict": "fail", "n_tests": 30, "failed": ["total_variance/trait1/sd_z"]},
        {"verdict": "pass", "n_tests": 30, "failed": []},
        {"verdict": "insufficient", "n_tests": 22, "failed": []},
    ]
    row = summarize("rm_base", "a_plus_005", outcomes).split(" | ")
    assert row[3] == "2/3 (+1 insufficient)"
    assert row[5] == "22/30"
    assert row[6] == "1/3"
    assert row[7] == "1 / 0 / 0"
    assert row[8] == "total_variance (2) |"
