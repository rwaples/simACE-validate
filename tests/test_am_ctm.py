"""C-transmission pieces of layer 5 (rwaples/simACE-validate#5): statistics, schema, law."""

import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from am.bulmer import predict_ctm
from am.manifest import CTM, CTM_HALFSIB, PER_TRAIT_INFO, defaults, schema
from am.stats import _half_sib_pairs, ctm_stats


def test_half_sib_pairs_split_by_the_other_parent():
    # mother 0 has children by fathers 10 and 11; father 11 also has a child by mother 1
    mother = np.array([0, 0, 0, 1])
    father = np.array([10, 10, 11, 11])
    mhs = {tuple(sorted(p)) for p in zip(*_half_sib_pairs(mother, father), strict=True)}
    phs = {tuple(sorted(p)) for p in zip(*_half_sib_pairs(father, mother), strict=True)}
    assert mhs == {(0, 2), (1, 2)}
    assert phs == {(2, 3)}


def test_ctm_stats_recovers_transmission_and_midparent_slope():
    rng = np.random.default_rng(1)
    n_par, n = 4000, 40000
    pc, py = rng.normal(size=n_par), rng.normal(size=n_par)
    mother = rng.integers(0, n_par // 2, n)
    father = rng.integers(n_par // 2, n_par, n)
    c = 0.5 * pc[mother] + 0.3 * pc[father] + np.sqrt(1 - 0.34) * rng.normal(size=n)
    y = 0.4 * (py[mother] + py[father]) + rng.normal(size=n)
    out = ctm_stats(c, y, mother, father, pc, py)
    assert out["r_MO_C"] == pytest.approx(0.5, abs=0.02)
    assert out["r_FO_C"] == pytest.approx(0.3, abs=0.02)
    assert out["b_MP_L"] == pytest.approx(0.8, abs=0.05)
    assert out["r_MHS_C"] == pytest.approx(0.25, abs=0.03)  # MHS share only b_m c_mother here
    assert out["r_PHS_C"] == pytest.approx(0.09, abs=0.03)


def test_first_generation_has_only_var_c():
    assert set(ctm_stats(np.arange(5.0), np.arange(5.0))) == {"var_C"}


def test_schema_scores_mo_fo_and_half_sibs_only_under_c_transmission():
    plain_scored, plain_info = schema(1, 6)
    scored, info = schema(1, 6, ctm=True, halfsib=True)
    lam_scored, _ = schema(1, 6, ctm=True)
    assert {f"{s}_g2" for s in PER_TRAIT_INFO} <= set(plain_info)
    assert not info
    assert len(scored) - len(plain_scored) == len(PER_TRAIT_INFO) * 5 + 6 + (len(CTM) + len(CTM_HALFSIB)) * 5
    assert not {f"{s}_g2" for s in CTM_HALFSIB} & set(lam_scored)


def test_law_matches_simace_shared_c():
    a, c, e, bm, bf = 0.5, 0.2, 0.3, 0.5, 0.3
    law = predict_ctm(a, c, e, bm, bf)
    # CLAUDE.md of simACE: gamma MO b_m, FO b_f, MHS 1 - b_f^2, PHS b_f^2; r = 2 phi A + gamma C over V_P
    assert law["r_MO_L"] == pytest.approx(0.25 + bm * c)
    assert law["r_PHS_L"] == pytest.approx(0.125 + bf**2 * c)
    assert law["r_MHS_L"] == pytest.approx(0.125 + (1 - bf**2) * c)
    assert law["b_MP_L"] == pytest.approx(a + (bm + bf) * c)


def test_defaults_leave_existing_specs_without_transmission():
    spec = defaults({"N": 10, "M": 10, "G": 6, "r": 0.5, "seed0": 1, "R": 2})
    assert (spec["bm"], spec["bf"], spec["pairing"]) == (0, 0, "lam")
