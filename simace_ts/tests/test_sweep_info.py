"""simace_ts/sweep_info.py: paired differences from the reference architecture."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from sweep_info import contrasts, joint_test


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


def _parts(rng, n, offsets):
    archs = ["a", "b", "c", "ref"]
    return [
        {
            "architectures": {
                a: {
                    "info": {
                        f"tetrachoric_{k}": float(rng.normal()) + common + offsets.get(a, 0.0)
                        for k, common in zip("XY", rng.normal(size=2), strict=True)
                    }
                    | {"h2_observed_falconer": None}
                }
                for a in archs
            }
        }
        for _ in range(n)
    ]


def test_joint_test_finds_a_shifted_reference_and_not_an_exchangeable_one():
    rng = np.random.default_rng(1)
    _, p_shift, names = joint_test(_parts(rng, 30, {"ref": 1.0}), "ref", draws=2000)
    assert names == ["tetrachoric_X", "tetrachoric_Y"]
    assert p_shift < 0.01
    _, p_null, _ = joint_test(_parts(rng, 30, {}), "ref", draws=2000)
    assert p_null > 0.01


def test_joint_test_flags_the_outlier_when_it_plays_the_reference():
    parts = _parts(np.random.default_rng(2), 30, {"a": 1.0})
    stat_a, p_a, _ = joint_test(parts, "a", draws=2000)
    stat_ref, _, _ = joint_test(parts, "ref", draws=2000)
    assert p_a < 0.01
    assert stat_a > stat_ref


def test_joint_test_drops_a_statistic_missing_in_one_seed():
    parts = _parts(np.random.default_rng(3), 10, {})
    parts[4]["architectures"]["b"]["info"]["tetrachoric_Y"] = None
    _, _, names = joint_test(parts, "ref", draws=100)
    assert names == ["tetrachoric_X"]


def test_joint_test_refuses_a_contrast_with_no_variance():
    parts = _parts(np.random.default_rng(4), 10, {})
    for p in parts:
        for info in (a["info"] for a in p["architectures"].values()):
            info["tetrachoric_X"] = 0.3
    with pytest.raises(ValueError, match="tetrachoric_X"):
        joint_test(parts, "ref", draws=100)


def test_joint_test_refuses_a_single_architecture():
    parts = [{"architectures": {"ref": {"info": {"tetrachoric_X": float(v)}}}} for v in range(5)]
    with pytest.raises(ValueError, match="two or more"):
        joint_test(parts, "ref", draws=100)
