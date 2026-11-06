r"""Overwrite pedigree A1 with rescaled tstrait GV.

Produces a sibling `pedigree.full.tstrait.parquet` where the A1 column for
sample-tagged individuals (the latest G_pheno generations) is replaced with
the tstrait genome-wide genetic value, centered to mean 0 and rescaled so
its sample variance matches the configured A1 variance. Non-sample
individuals (older ancestors not in the genotyped sample set) keep their
original parametric A1. With ``--per-generation`` each generation is
centred and rescaled on its own, so every generation has the configured
variance (simACE-validate ADR 0002); without it the sample is pooled.

`liability1 = A1 + C1 + E1` is recomputed after the overwrite so it stays
consistent. C1 and E1 come from the original pedigree unchanged — the goal
is to use simACE's standard parametric C/E components alongside a realistic
A from the genotypes.

This makes the tstrait branch comparable to standard simACE scenarios:
downstream phenotype models (frailty, threshold, etc.) can be run on the
augmented pedigree and use the simACE-style A+C+E variance composition.

``--pedigree`` comes from the ``drop_from`` scenario when one is set;
``--target-var-a`` is the scenario's ``A1``. Runs in this repo's pixi env
(``pixi run python``). Example::

    python tskit_pipeline/scripts/tstrait_augment_pedigree.py \
        --pedigree results/{folder}/{scenario}/rep{rep}/pedigree.parquet \
        --gv results/{folder}/{scenario}/rep{rep}/gv_chrom_{1..22}.parquet \
        --out-pedigree results/{folder}/{scenario}/rep{rep}/pedigree.full.tstrait.parquet \
        --meta results/{folder}/{scenario}/rep{rep}/pedigree.full.tstrait_meta.json \
        --target-var-a 0.5 \
        --log logs/{folder}/{scenario}/rep{rep}/tstrait_augment_pedigree.log
"""

import argparse
import json
import logging
import re
from pathlib import Path

import numpy as np
import pandas as pd
import polars as pl
from _script_logging import route_logging_to_file
from simace.core.cohort import read_c_transmission_clips, read_pedigree, write_pedigree

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s | %(message)s")
log = logging.getLogger("tstrait_augment_pedigree")


def _chrom_key(p: Path) -> int:
    m = re.search(r"chrom_?([0-9]+)", p.stem, re.IGNORECASE)
    if not m:
        raise ValueError(f"could not parse chromosome from {p.name}")
    return int(m.group(1))


def rescale_gv(gv: np.ndarray, target_var: float) -> tuple[np.ndarray, dict]:
    """Center GV at 0 and rescale so its sample variance == target_var exactly.

    Sample variance uses ddof=0 (population variance) so the post-rescale
    `np.var(out, ddof=0) == target_var` holds to floating-point precision.
    """
    if target_var < 0:
        raise ValueError(f"target_var must be >= 0; got {target_var}")
    realized_mean = float(gv.mean())
    realized_var = float(gv.var(ddof=0))
    if realized_var == 0:
        raise ValueError("realized GV variance is 0; cannot rescale")
    scale = float(np.sqrt(target_var / realized_var))
    out = (gv - realized_mean) * scale
    info = {
        "realized_mean": realized_mean,
        "realized_var": realized_var,
        "target_var": float(target_var),
        "scale_factor": scale,
        "rescaled_var": float(out.var(ddof=0)),
        "rescaled_mean": float(out.mean()),
    }
    return out, info


def rescale_gv_per_generation(
    gv: np.ndarray, generation: np.ndarray, target_var: float
) -> tuple[np.ndarray, list[dict]]:
    """Apply rescale_gv within each generation, so every generation has variance target_var (ADR 0002).

    Returns the rescaled values and one rescale_gv info dict per generation, with its ``generation``.
    """
    out = np.empty_like(gv, dtype=np.float64)
    info = []
    for g in np.unique(generation).tolist():
        rows = generation == g
        try:
            out[rows], row = rescale_gv(gv[rows], target_var)
        except ValueError as err:
            raise ValueError(f"generation {g} ({int(rows.sum())} individuals): {err}") from err
        info.append({"generation": g, **row})
    return out, info


def write_like(df: pd.DataFrame, path: Path, source: Path) -> None:
    """Write ``df`` as a recorded pedigree of ``source``'s results layout.

    ``simace cohort`` / ``simace analyze`` refuse a pedigree without the
    layout marker (ADR 0021), so a ``source`` without one is refused here.
    ``write_pedigree`` writes the marker, narrows dtypes as ``simace simulate``
    does, and carries ``source``'s C-transmission clip records.
    """
    read_pedigree(source, columns=["id"])
    write_pedigree(pl.from_pandas(df), path, read_c_transmission_clips(source))


def main(argv: list[str] | None = None) -> None:
    """Rescale + overwrite A1 + recompute liability1."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--pedigree", required=True)
    parser.add_argument("--gv", nargs="+", required=True)
    parser.add_argument("--out-pedigree", required=True)
    parser.add_argument("--meta", required=True)
    parser.add_argument("--target-var-a", type=float, required=True)
    parser.add_argument(
        "--per-generation",
        action="store_true",
        help="centre and rescale each generation separately instead of pooling them",
    )
    parser.add_argument("--log", default=None)
    args = parser.parse_args(argv)

    route_logging_to_file(args.log)
    pedigree_path = Path(args.pedigree)
    gv_paths = sorted([Path(p) for p in args.gv], key=_chrom_key)
    out_pedigree = Path(args.out_pedigree)
    out_meta = Path(args.meta)
    target_var_a = args.target_var_a
    out_pedigree.parent.mkdir(parents=True, exist_ok=True)

    log.info("loading pedigree %s", pedigree_path)
    ped = pd.read_parquet(pedigree_path)
    log.info("  pedigree: %d rows, columns: %s", len(ped), list(ped.columns))

    log.info("summing per-chrom GVs across %d files", len(gv_paths))
    gv_dfs = [pd.read_parquet(p) for p in gv_paths]
    gvs = pd.concat(gv_dfs, ignore_index=True).groupby("individual_id", as_index=False)["genetic_value"].sum()
    log.info("  GV summed: %d unique individuals", len(gvs))

    n_in_gv_not_in_ped = int((~gvs["individual_id"].isin(ped["id"])).sum())
    if n_in_gv_not_in_ped:
        raise ValueError(
            f"{n_in_gv_not_in_ped} GV individual_ids not present in pedigree.id — "
            "drop+graft and pedigree are out of sync"
        )
    gv = gvs["genetic_value"].to_numpy()
    rescale_by_generation = None
    if args.per_generation:
        # simACE-ts (ADR 0002): every row gets gene-drop A, and GV ids are simACE ids by row order.
        if not np.array_equal(ped["id"].to_numpy(), np.arange(len(ped))):
            raise ValueError("pedigree ids must equal row order 0..n-1")
        if "twin" in ped.columns and (ped["twin"].to_numpy() >= 0).any():
            raise ValueError("pedigree has MZ twins; simulate with p_mztwin = 0")
        if len(gvs) != len(ped):
            raise ValueError(
                f"--per-generation needs a genetic value for every row: {len(ped) - len(gvs)} of {len(ped)} "
                "would keep parametric A1 (build the tables with --g-pheno equal to --g-ped)"
            )
        generation = ped["generation"].to_numpy()[gvs["individual_id"].to_numpy()]
        rescaled_arr, rescale_by_generation = rescale_gv_per_generation(gv, generation, target_var_a)
        rescale_info = None
        for row in rescale_by_generation:
            log.info(
                "  rescale g%d: realized_var=%.4f -> target_var=%.4f (scale=%.4f)",
                row["generation"],
                row["realized_var"],
                row["target_var"],
                row["scale_factor"],
            )
    else:
        rescaled_arr, rescale_info = rescale_gv(gv, target_var_a)
        log.info(
            "  rescale: realized_var=%.4f -> target_var=%.4f (scale=%.4f)",
            rescale_info["realized_var"],
            rescale_info["target_var"],
            rescale_info["scale_factor"],
        )

    gv_map = pd.Series(rescaled_arr, index=gvs["individual_id"].to_numpy())
    mask = ped["id"].isin(gv_map.index)
    n_overwritten = int(mask.sum())

    a1_orig = ped["A1"].to_numpy().copy()
    a1_new = a1_orig.copy()
    a1_new[mask.to_numpy()] = gv_map.reindex(ped.loc[mask, "id"].to_numpy()).to_numpy()

    ped_out = ped.copy()
    ped_out["A1"] = a1_new.astype(np.float32)
    if "liability1" in ped_out.columns:
        ped_out["liability1"] = (
            ped_out["A1"].astype(np.float64) + ped_out["C1"].astype(np.float64) + ped_out["E1"].astype(np.float64)
        )

    log.info("writing %s", out_pedigree)
    write_like(ped_out, out_pedigree, pedigree_path)

    sample_mask = mask.to_numpy()
    meta = {
        "n_pedigree_rows": len(ped),
        "n_inds_overwritten": n_overwritten,
        "n_inds_kept_parametric": int((~mask).sum()),
        "rescale": rescale_info,  # pooled mode; None with --per-generation
        "rescale_by_generation": rescale_by_generation,  # --per-generation; None when pooled
        "var_A1_after_in_overwritten": float(np.var(a1_new[sample_mask], ddof=0)),
        "mean_A1_after_in_overwritten": float(np.mean(a1_new[sample_mask])),
        "var_A1_orig_in_overwritten": float(np.var(a1_orig[sample_mask], ddof=0)),
        "var_A1_orig_in_kept": (float(np.var(a1_orig[~sample_mask], ddof=0)) if (~sample_mask).any() else None),
    }
    out_meta.write_text(json.dumps(meta, indent=2, default=float))
    log.info(
        "overwrote A1 for %d / %d individuals; var(A1_new) in overwritten = %.4f",
        n_overwritten,
        len(ped),
        meta["var_A1_after_in_overwritten"],
    )


if __name__ == "__main__":
    main()
