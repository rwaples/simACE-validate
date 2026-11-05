"""simace_ts/sweep_info.py: paired differences from the reference architecture."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sweep_info import contrasts


def _part(ref: float, other: float, h2: float | None) -> dict:
    return {
        "architectures": {
            "ref": {"info": {"tetrachoric_FS": ref, "h2_observed_falconer": h2, "gd_rA_FS_minus_2phi": 0.0}},
            "other": {"info": {"tetrachoric_FS": other, "h2_observed_falconer": None, "gd_rA_FS_minus_2phi": 9.0}},
        }
    }


def test_contrasts_are_paired_and_skip_empty_and_unlisted_stats():
    parts = [_part(0.5, 0.4, None), _part(0.6, 0.4, None), _part(0.7, 0.6, 0.3)]
    [(stat, arch, mean, se, t, n)] = contrasts(parts, "ref")
    assert (stat, arch, n) == ("tetrachoric_FS", "other", 3)
    assert mean == pytest.approx(-0.4 / 3)
    assert se == pytest.approx(0.0577350 / 3**0.5, rel=1e-5)
    assert t == pytest.approx(mean / se)
