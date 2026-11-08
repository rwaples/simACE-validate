"""simace_ts/g6.py and score.load_g6: the layer-2 cross-check (ADR 0002 item 11)."""

from __future__ import annotations

import copy
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[2]
sys.path[:0] = [str(ROOT / "simace_ts"), str(ROOT / "tskit_pipeline" / "scripts")]
import g6  # noqa: E402
import polars as pl  # noqa: E402
import score  # noqa: E402
import tskit  # noqa: E402
from build_pedigree_tables import pedigree_to_msprime  # noqa: E402
from part import class_pairs, pearson  # noqa: E402
from pedigree_graph import PedigreeGraph  # noqa: E402


def _pedigree(rng: np.random.Generator, n: int = 40, generations: int = 3) -> pd.DataFrame:
    rows, gen_ids = [], []
    for g in range(generations):
        ids = list(range(len(rows), len(rows) + n))
        if g == 0:
            rows += [(i, -1, -1, g, i % 2) for i in ids]
        else:
            prev = gen_ids[-1]
            mothers, fathers = [i for i in prev if i % 2 == 0], [i for i in prev if i % 2 == 1]
            couples = [(mothers[k], fathers[rng.integers(len(fathers))]) for k in range(len(mothers))]
            for i in ids:
                m, f = couples[rng.integers(len(couples))]
                rows.append((i, m, f, g, i % 2))
        gen_ids.append(ids)
    return pd.DataFrame(rows, columns=["id", "mother", "father", "generation", "sex"])


def test_tskit_spans_give_half_for_parent_and_child_of_unrelated_mates(tmp_path):
    ped = _pedigree(np.random.default_rng(1))
    tables = tmp_path / "pedigree_tables.trees"
    pedigree_to_msprime(ped, 3, 3, 40).finalise(sequence_length=1).dump(str(tables))
    kids = ped[ped["mother"] >= 0]
    pairs = np.column_stack([np.r_[kids["mother"], kids["father"]], np.r_[kids["id"], kids["id"]]])
    graph = PedigreeGraph.from_arrays(ids=ped["id"].to_numpy(), mother_ids=ped["mother"].to_numpy(),
                                      father_ids=ped["father"].to_numpy())  # fmt: skip
    two_phi = 2 * np.asarray(graph.pair_kinship(pairs[:, 0], pairs[:, 1]), dtype=np.float64)
    length = 10_000_000
    ibd = g6.chromosome_spans(tskit.TableCollection.load(str(tables)), length, 1e-8, 7, pairs) / (2 * length)
    outbred = two_phi == 0.5
    assert outbred.sum() > 20
    assert (~outbred).sum() > 0  # the third generation has some related mates
    np.testing.assert_array_equal(ibd[outbred], 0.5)
    assert np.all(ibd[~outbred] >= 0.5)


def test_sample_pairs_keep_part_order_so_the_tskit_correlation_matches_the_part(tmp_path):
    ped = pl.from_pandas(_pedigree(np.random.default_rng(5), n=400, generations=4))
    ped.write_parquet(tmp_path / "pedigree.full.tstrait.parquet")
    keep = ped["id"].to_numpy()[ped["generation"].to_numpy() >= 1]
    pl.DataFrame({"id": keep, "affected1": np.ones(len(keep), dtype=bool)}).write_parquet(tmp_path / "cohort.parquet")
    a = np.random.default_rng(6).normal(size=ped.height)
    part_pairs = class_pairs(ped.filter(pl.col("generation") >= 1))  # part.py's selection
    got = g6.sample_pairs(tmp_path, g6.read_pedigree(tmp_path))
    for k, (i, j) in part_pairs.items():
        assert pearson(a[got[k][0]], a[got[k][1]]) == pytest.approx(pearson(a[i], a[j]), abs=score.G6_TOL, nan_ok=True)
        assert sorted(zip(got[k][0], got[k][1], strict=True)) == sorted(zip(i, j, strict=True))


def test_subsample_is_sorted_reproducible_and_capped():
    pairs = {"FS": (np.arange(10)[::-1], np.arange(10)[::-1] + 100), "MO": (np.arange(3), np.arange(3) + 50)}
    a = g6.subsample(pairs, 4, seed=9000)
    b = g6.subsample(pairs, 4, seed=9000)
    assert all(np.array_equal(a[k], b[k]) for k in a)
    assert a["FS"].shape == (4, 2)
    assert a["MO"].shape == (3, 2)
    assert np.all(np.diff(a["FS"][:, 0]) > 0)
    assert not np.array_equal(a["FS"], g6.subsample(pairs, 4, seed=9001)["FS"])


CASE = {"id": "ts_test", "architectures": [{"name": "first"}, {"name": "second"}]}
SIMACE = {"version": "0.3.0"}
GENOME = {"chromosomes": 22, "morgans": 1.75, "causal_per_chromosome": 45}


def _part(seed: int, rng: np.random.Generator) -> dict:
    info = {f"gd_rA_{k}_minus_2phi": float(rng.normal(0, 0.01)) for k in score.PAIR_CLASSES}
    return {"manifest": {"seed": seed, "simace": SIMACE},
            "architectures": {"first": {"info": info}, "second": {"info": copy.deepcopy(info)}}}  # fmt: skip


def _g6(part: dict, part_bytes: bytes) -> dict:
    s = part["manifest"]["seed"]
    rA = {
        arch: {k: {"l2": 0.001 + (v["info"][f"gd_rA_{k}_minus_2phi"] or 0.0), "ts": v["info"][f"gd_rA_{k}_minus_2phi"],
                   "n": 9} for k in score.PAIR_CLASSES}
        for arch, v in part["architectures"].items()
    }  # fmt: skip
    ibd = {k: {"n": 9, "two_phi_minus_nominal": 1e-4, "ibd_l2": {"mean": 0.0, "se": 1.0},
               "ibd_ts": {"mean": s * 1e-6, "se": 1.0}} for k in score.PAIR_CLASSES}  # fmt: skip
    return {"g6_version": score.G6_VERSION, "g6_by": "abc", "seed": s, "case": CASE["id"], "simace": SIMACE,
            "part_sha256": hashlib.sha256(part_bytes).hexdigest(), "pairs_per_class": 2000, "genome": GENOME,
            "ibd": ibd, "rA": rA}  # fmt: skip


@pytest.fixture
def setup(tmp_path):
    rng = np.random.default_rng(3)
    parts = [_part(s, rng) for s in (100, 101, 102)]
    parts[0]["architectures"]["second"]["info"]["gd_rA_MO_minus_2phi"] = None  # an empty class, as part.num leaves it
    d, pd_ = tmp_path / "g6" / CASE["id"], tmp_path / "parts" / CASE["id"]
    d.mkdir(parents=True)
    pd_.mkdir(parents=True)
    for p in parts:
        name = f"seed{p['manifest']['seed']}.json"
        (pd_ / name).write_text(json.dumps(p))
        (d / name).write_text(json.dumps(_g6(p, (pd_ / name).read_bytes())))
    return parts, tmp_path


def _load(parts: list[dict], root: Path) -> dict[str, np.ndarray]:
    return score.load_g6(CASE, parts, root / "g6", root / "parts")


def test_load_g6_without_a_directory_gives_nothing(tmp_path):
    assert score.load_g6(CASE, [], tmp_path) == {}


def test_load_g6_gives_values_per_seed(setup):
    parts, d = setup
    out = _load(parts, d)
    np.testing.assert_allclose(out["ibd_ts_FS"], [1e-4, 1.01e-4, 1.02e-4])
    np.testing.assert_allclose(out["rA_l2_minus_ts_MO"], 0.001, atol=1e-12)
    assert np.isnan(out["second/rA_l2_minus_ts_MO"][0])
    assert np.isfinite(out["second/rA_l2_minus_ts_MO"][1:]).all()
    assert "first/rA_l2_minus_ts_MO" not in out


@pytest.mark.parametrize(
    ("change", "message"),
    [
        (lambda g: g.update(g6_by="abc-dirty"), "uncommitted"),
        (lambda g: g.update(g6_version=0), "g6_version"),
        (lambda g: g.update(seed=7), "seed or simace"),
        (lambda g: g.update(case="other"), "seed or simace"),
        (lambda g: g.update(part_sha256="0" * 64), "part_sha256"),
        (lambda g: g["rA"]["second"]["PHS"].update(ts=g["rA"]["second"]["PHS"]["ts"] + 1e-9), "second/PHS"),
        (lambda g: g["rA"]["first"]["PHS"].update(ts=None), "first/PHS"),
        (lambda g: g.update(pairs_per_class=500), "pairs_per_class or genome"),
        (lambda g: g["genome"].update(morgans=2.0), "pairs_per_class or genome"),
    ],
)
def test_load_g6_refuses_a_stale_file(setup, change, message):
    parts, d = setup
    path = d / "g6" / CASE["id"] / "seed101.json"
    g = json.loads(path.read_text())
    change(g)
    path.write_text(json.dumps(g))
    with pytest.raises(ValueError, match=message):
        _load(parts, d)


def test_load_g6_refuses_a_missing_seed(setup):
    parts, d = setup
    (d / "g6" / CASE["id"] / "seed102.json").unlink()
    with pytest.raises(FileNotFoundError, match="seed 102"):
        _load(parts, d)
