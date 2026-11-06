"""Pipeline changes for simACE-ts (ADR 0002): per-generation rescale, and the guards on its inputs."""

from __future__ import annotations

import json

import build_pedigree_tables as build
import numpy as np
import pandas as pd
import polars as pl
import pytest
import tstrait_augment_pedigree as augment
from simace.core.cohort import read_pedigree, write_pedigree


def _pedigree(n_per_gen: int = 4, gens: int = 3) -> pd.DataFrame:
    """A small pedigree: each generation's children have parents in the one before."""
    n = n_per_gen * gens
    gen = np.repeat(np.arange(gens), n_per_gen)
    mother = np.full(n, -1)
    father = np.full(n, -1)
    for i in range(n_per_gen, n):
        base = (gen[i] - 1) * n_per_gen
        mother[i], father[i] = base + 2 * (i % 2), base + 1 + 2 * (i % 2)
    return pd.DataFrame(
        {
            "id": np.arange(n),
            "mother": mother,
            "father": father,
            "generation": gen,
            "twin": np.full(n, -1),
            "A1": np.full(n, 9.0),
            "C1": np.full(n, 0.1),
            "E1": np.full(n, 0.2),
            "liability1": np.full(n, 9.3),
        }
    )


def _run_augment(tmp_path, ped: pd.DataFrame, gv_ids: np.ndarray, *flags: str) -> tuple[pd.DataFrame, dict]:
    source = tmp_path / "pedigree.parquet"
    write_pedigree(pl.from_pandas(ped), source)
    rng = np.random.default_rng(2)
    gv_path = tmp_path / "gv_chrom_1.parquet"
    pd.DataFrame({"individual_id": gv_ids, "genetic_value": rng.normal(size=len(gv_ids)) * 7 + 3}).to_parquet(gv_path)
    out, meta = tmp_path / "out.parquet", tmp_path / "meta.json"
    augment.main(
        ["--pedigree", str(source), "--gv", str(gv_path), "--out-pedigree", str(out),
         "--meta", str(meta), "--target-var-a", "0.4", *flags]
    )  # fmt: skip
    return read_pedigree(out).to_pandas(), json.loads(meta.read_text())


def test_per_generation_rescale_sets_each_generation():
    rng = np.random.default_rng(1)
    generation = np.repeat([0, 1, 2], 50)
    gv = rng.normal(size=150) * np.repeat([1.0, 3.0, 0.5], 50) + np.repeat([5.0, -2.0, 0.0], 50)

    out, info = augment.rescale_gv_per_generation(gv, generation, 0.4)

    for g in (0, 1, 2):
        rows = generation == g
        assert np.var(out[rows]) == pytest.approx(0.4, rel=1e-12)
        assert np.mean(out[rows]) == pytest.approx(0.0, abs=1e-12)
        assert info[g]["generation"] == g
        assert info[g]["realized_var"] == pytest.approx(np.var(gv[rows]))
    # Within a generation the rescale is affine, so the order of values is kept.
    assert np.array_equal(np.argsort(out[:50]), np.argsort(gv[:50]))


def test_zero_variance_generation_is_named():
    with pytest.raises(ValueError, match="generation 1"):
        augment.rescale_gv_per_generation(np.array([1.0, 2.0, 5.0, 5.0]), np.array([0, 0, 1, 1]), 0.4)


def test_main_per_generation_overwrites_every_row(tmp_path):
    ped = _pedigree()
    got, meta = _run_augment(tmp_path, ped, ped["id"].to_numpy(), "--per-generation")

    for g in range(3):
        a = got.loc[got["generation"] == g, "A1"].to_numpy(np.float64)
        assert np.var(a) == pytest.approx(0.4, rel=1e-5)  # A1 is stored as float32
        assert np.mean(a) == pytest.approx(0.0, abs=1e-6)
    np.testing.assert_allclose(got["liability1"], got["A1"] + got["C1"] + got["E1"], rtol=1e-6)
    assert meta["rescale"] is None
    assert [r["generation"] for r in meta["rescale_by_generation"]] == [0, 1, 2]
    assert meta["n_inds_kept_parametric"] == 0


def test_main_per_generation_refuses_rows_left_parametric(tmp_path):
    ped = _pedigree()
    sampled = ped.loc[ped["generation"] >= 1, "id"].to_numpy()  # generation 0 not a sample (G_pheno < G_ped)
    with pytest.raises(ValueError, match="would keep parametric A1"):
        _run_augment(tmp_path, ped, sampled, "--per-generation")


def test_main_per_generation_refuses_twins(tmp_path):
    ped = _pedigree()
    ped.loc[[5, 7], "twin"] = [7, 5]
    with pytest.raises(ValueError, match="MZ twins"):
        _run_augment(tmp_path, ped, ped["id"].to_numpy(), "--per-generation")


def test_main_pooled_mode_is_unchanged(tmp_path):
    ped = _pedigree()
    sampled = ped.loc[ped["generation"] >= 1, "id"].to_numpy()
    got, meta = _run_augment(tmp_path, ped, sampled)

    assert np.var(got.loc[got["generation"] >= 1, "A1"].to_numpy(np.float64)) == pytest.approx(0.4, rel=1e-5)
    assert np.all(got.loc[got["generation"] == 0, "A1"] == 9.0)
    assert meta["rescale_by_generation"] is None
    assert meta["rescale"]["target_var"] == 0.4


def test_build_refuses_twins():
    ped = _pedigree()
    ped.loc[[5, 7], "twin"] = [7, 5]
    with pytest.raises(ValueError, match="MZ twins"):
        build.pedigree_to_msprime(ped, G_ped=3, G_pheno=3, Ne=4)


def test_build_refuses_ids_out_of_row_order():
    ped = _pedigree()
    ped["id"] = ped["id"] + 1
    with pytest.raises(ValueError, match="row order"):
        build.pedigree_to_msprime(ped, G_ped=3, G_pheno=3, Ne=4)


def test_build_refuses_founders_after_non_founders():
    ped = _pedigree()
    ped.loc[2, ["mother", "father"]] = [0, 1]  # row 2 becomes a child while row 3 stays a founder
    ped.loc[ped["mother"] == 2, "mother"] = 0
    with pytest.raises(ValueError, match="first rows"):
        build.pedigree_to_msprime(ped, G_ped=3, G_pheno=3, Ne=4)


def test_build_marks_every_generation_as_sample():
    tables = build.pedigree_to_msprime(_pedigree(), G_ped=3, G_pheno=3, Ne=4).finalise(sequence_length=1)
    assert np.all(tables.nodes.flags & 1)  # tskit.NODE_IS_SAMPLE
