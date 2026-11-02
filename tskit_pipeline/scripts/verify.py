r"""Verify the canonicalized concat output.

Two checks:
  1. Sidecar self-consistency: recompute the BLAKE2b fingerprint of the
     concat tables and compare to the sidecar `<trees>.fingerprint` written
     by concat_chroms.py. Detects on-disk corruption / incomplete writes.
  2. Golden expectations: read `tests/data/tskit_preprocess_fingerprint.json`
     and assert n_eligible_individuals, num_chromosomes, total_seqlen, and
     per-chrom (n_samples, n_sites, n_muts) match. Detects source-data drift.

Bootstrap: set env var `SIMACE_TSKIT_WRITE_GOLDEN=1` to write the current
summary into the golden file instead of asserting against it. Use this on
the first run for a source dataset, or after a deliberate source-data
update, and keep the golden file beside the preprocessed outputs.

Runs in this repo's pixi env (``pixi run python``). Example::

    python tskit_pipeline/scripts/verify.py \
        --trees {output_dir}/all_chroms.trees \
        --fingerprint {output_dir}/all_chroms.trees.fingerprint \
        --summary {output_dir}/preprocess_summary.json \
        --golden {output_dir}/preprocess_golden.json \
        --log logs/tskit_preprocess/verify.log
"""

import argparse
import json
import logging
import os
from pathlib import Path

import numpy as np
import tskit
from _script_logging import route_logging_to_file
from concat_chroms import fingerprint_tables

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s | %(message)s")
log = logging.getLogger("verify")


def golden_from_summary(summary: dict) -> dict:
    """Build the canonical golden dict from a preprocess summary."""
    return {
        "n_eligible_individuals": int(summary["consistency"]["n_eligible_individuals"]),
        "num_chromosomes": int(summary["n_chromosomes"]),
        "total_seqlen": float(summary["totals"]["total_seqlen"]),
        "per_chromosome": [
            {
                "name": r["name"],
                "n_samples": int(r["num_samples_post"]),
                "n_sites": int(r["num_sites_post"]),
                "n_muts": int(r["num_mutations_post"]),
            }
            for r in summary["per_chromosome"]
        ],
    }


def assert_golden(observed: dict, golden: dict) -> None:
    """Diff observed against the checked-in golden; raise if anything mismatches."""
    diffs: list[str] = [
        f"{k}: observed={observed[k]} golden={golden[k]}"
        for k in ("n_eligible_individuals", "num_chromosomes")
        if observed[k] != golden[k]
    ]
    if not np.isclose(observed["total_seqlen"], golden["total_seqlen"]):
        diffs.append(f"total_seqlen: observed={observed['total_seqlen']:.3f} golden={golden['total_seqlen']:.3f}")
    obs_pc = {r["name"]: r for r in observed["per_chromosome"]}
    gold_pc = {r["name"]: r for r in golden["per_chromosome"]}
    if set(obs_pc) != set(gold_pc):
        diffs.append(f"chromosome set differs: observed={sorted(obs_pc)} golden={sorted(gold_pc)}")
    diffs.extend(
        f"{name}.{k}: observed={obs_pc[name][k]} golden={gold_pc[name][k]}"
        for name in sorted(set(obs_pc) & set(gold_pc))
        for k in ("n_samples", "n_sites", "n_muts")
        if obs_pc[name][k] != gold_pc[name][k]
    )
    if diffs:
        raise AssertionError(
            "tskit_preprocess golden mismatch (set SIMACE_TSKIT_WRITE_GOLDEN=1 to update):\n  " + "\n  ".join(diffs)
        )


def main(argv: list[str] | None = None) -> None:
    """Golden + sidecar fingerprint checks."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--trees", required=True)
    parser.add_argument("--fingerprint", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--golden", required=True)
    parser.add_argument("--log", default=None)
    args = parser.parse_args(argv)

    route_logging_to_file(args.log)
    trees_path = Path(args.trees)
    fp_path = Path(args.fingerprint)
    summary_path = Path(args.summary)
    golden_path = Path(args.golden)

    summary = json.loads(summary_path.read_text())
    observed = golden_from_summary(summary)

    if os.environ.get("SIMACE_TSKIT_WRITE_GOLDEN") == "1":
        golden_path.parent.mkdir(parents=True, exist_ok=True)
        golden_path.write_text(json.dumps(observed, indent=2, sort_keys=True) + "\n")
        log.warning("WROTE GOLDEN to %s (SIMACE_TSKIT_WRITE_GOLDEN=1). Commit it.", golden_path)
    else:
        if not golden_path.exists():
            raise FileNotFoundError(
                f"golden file {golden_path} missing — run once with SIMACE_TSKIT_WRITE_GOLDEN=1 "
                "to bootstrap, then commit."
            )
        golden = json.loads(golden_path.read_text())
        assert_golden(observed, golden)
        log.info(
            "golden check passed (n_eligible=%d, n_chrom=%d, total_seqlen=%.1f Mb)",
            observed["n_eligible_individuals"],
            observed["num_chromosomes"],
            observed["total_seqlen"] / 1e6,
        )

    sidecar_hex = fp_path.read_text().strip()
    log.info("recomputing fingerprint from %s", trees_path)
    ts_all = tskit.load(str(trees_path))
    recomputed = fingerprint_tables(ts_all)
    if recomputed != sidecar_hex:
        raise AssertionError(f"sidecar fingerprint mismatch: sidecar={sidecar_hex} recomputed={recomputed}")
    log.info("sidecar fingerprint check passed (%s)", sidecar_hex)


if __name__ == "__main__":
    main()
