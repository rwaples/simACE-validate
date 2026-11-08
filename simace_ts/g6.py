"""G6: the layer-2 cross-check on one simACE-ts seed's pedigree (ADR 0002 item 11).

Per seed, on the pedigree simACE recorded (``<arch>/pedigree.full.tstrait.parquet``):

- **IBD of both engines against the recorded-pedigree kinship.** A sample of
  each direct-link class's analysis-sample pairs (sorted by (i, j), then drawn from
  ``stage_seed(s, STAGE_G6_PAIRS, class)``) gets its IBD fraction from two
  drops: layer 2's ``genedrop.drop`` (22 unlinked chromosomes of the mean
  SimHumanity map length, to 1/64 Morgan so the length in bp is an integer),
  and the tskit layer's own drop, rerun from the seed's
  ``pedigree_tables.trees`` with each chromosome's HapMapII map and the
  seed ``run_seed.py`` gave it (``stage_seed(s, STAGE_DROP, c)``), before
  the founder graft, so IBD is within the pedigree only. Each is compared
  with 2 phi from ``PedigreeGraph.pair_kinship``, which counts parental
  inbreeding and every path (for a parent and child,
  2 phi = (1 + F_parent) / 2 + phi(mother, father)). Expectation 0 for each
  engine.
- **A correlations, layer 2 minus tskit.** Per architecture and class, the
  Pearson correlation of layer 2's G over the architecture's analysis-sample
  pairs minus that of the architecture's A. This is a diagnostic with no
  derived null: the tskit A is centred and rescaled per generation, comes
  from a fixed founder panel with real allele frequencies and linkage, and a
  correlation is a ratio. Its tskit half is computed over the pairs in
  ``part.py``'s order and must match the Part's ``gd_rA_<class>_minus_2phi``;
  the file also records the sha256 of the seed's ``part.json``, which binds
  it to that Part.

Writes ``<seed dir>/g6.json``; ``score.load_g6`` reads the committed copies.

Usage::

    pixi run python simace_ts/g6.py <seed dir> --data ~/data/humanity_sim [--threads 4] [--pairs 2000]
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import polars as pl
import tskit

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / "tskit_pipeline" / "scripts")]
from genotype_drop_chrom import drop_one_chromosome, load_simhumanity_rate_map  # noqa: E402
from part import TWO_PHI, class_pairs, num, pearson  # noqa: E402
from pedigree_graph import PedigreeGraph  # noqa: E402
from run_seed import STAGE_DROP, STAGE_G6_DROP, STAGE_G6_PAIRS, Data, load_case, stage_seed, validate_commit  # noqa: E402
from score import CHROMOSOMES, G6_VERSION  # noqa: E402

from genedrop.drop import Genome, drop, individual_nodes, pair_spans  # noqa: E402

if TYPE_CHECKING:
    import msprime

CAUSAL_PER_CHROMOSOME = 45


def read_pedigree(d: Path) -> pl.DataFrame:
    """Return one architecture's recorded pedigree, refusing ids that are not row positions."""
    ped = pl.read_parquet(d / "pedigree.full.tstrait.parquet")
    if not np.array_equal(ped["id"].to_numpy(), np.arange(ped.height)):
        raise ValueError(f"{d}: pedigree ids must equal row order")
    return ped


def sample_pairs(d: Path, ped: pl.DataFrame) -> dict[str, tuple[np.ndarray, np.ndarray]]:
    """Return each class's analysis-sample pairs in ``d``, in part.py's order (so A correlations match it)."""
    cohort = pl.read_parquet(d / "cohort.parquet")
    in_sample = np.zeros(ped.height, dtype=bool)
    in_sample[cohort.filter(pl.col("affected1").is_not_null())["id"].to_numpy()] = True
    return class_pairs(ped.filter(pl.Series(in_sample)))


def subsample(pairs: dict[str, tuple[np.ndarray, np.ndarray]], per_class: int, seed: int) -> dict[str, np.ndarray]:
    """Return each class's subsample as rows of (i, j), sorted, drawn without replacement from its own seed namespace.

    The pairs are sorted first, since the join that found them has no guaranteed order.
    """
    out = {}
    for c, (k, (i, j)) in enumerate(pairs.items()):
        order = np.lexsort((j, i))
        rng = np.random.default_rng(stage_seed(seed, STAGE_G6_PAIRS, c))
        pick = order[np.sort(rng.choice(len(i), size=min(per_class, len(i)), replace=False))]
        out[k] = np.column_stack([i[pick], j[pick]])
    return out


def load_maps(data: Data) -> dict[int, tuple[int, msprime.RateMap]]:
    """Return each chromosome's length (bp) and HapMapII rate map."""
    base = data.simhumanity / "stdpopsim extraction" / "extracted"
    out = {}
    for c in CHROMOSOMES:
        length = int((base / f"chr{c}_length.txt").read_text().strip())
        out[c] = (length, load_simhumanity_rate_map(data.simhumanity, c, length))
    return out


def chromosome_spans(
    tables: tskit.TableCollection, length: int, rate: msprime.RateMap, seed: int, pairs: np.ndarray
) -> np.ndarray:
    """Return each pair's IBD span (bp, summed over its four node pairs) in one fixed-pedigree drop of ``tables``."""
    t = tables.copy()
    t.sequence_length = length
    ts = drop_one_chromosome(t, rate, seed)
    return pair_spans(ts, individual_nodes(ts), pairs).sum(axis=(1, 2))  # msprime individual ids follow pedigree rows


def tskit_ibd(
    tables: tskit.TableCollection, maps: dict[int, tuple[int, msprime.RateMap]], seed: int, pairs: np.ndarray,
    threads: int,
) -> np.ndarray:  # fmt: skip
    """Return each pair's IBD fraction (2 phi scale) in the tskit layer's drop, rerun before the graft."""

    def one(c: int) -> np.ndarray:
        length, rate = maps[c]
        return chromosome_spans(tables, length, rate, stage_seed(seed, STAGE_DROP, c), pairs)

    with ThreadPoolExecutor(threads) as pool:  # msprime releases the GIL; each chromosome has its own seed
        shared = np.sum(list(pool.map(one, CHROMOSOMES)), axis=0)
    return shared / (2 * sum(length for length, _ in maps.values()))


def mean_morgans(maps: dict[int, tuple[int, msprime.RateMap]]) -> float:
    """Return the mean map length of the chromosomes, rounded to 1/64 Morgan (so the length in bp is an integer)."""
    return round(float(np.mean([rate.total_mass for _, rate in maps.values()])) * 64) / 64


def gap_summary(ibd: np.ndarray, two_phi: np.ndarray) -> dict[str, float]:
    """Return the mean and SE over pairs of IBD - 2 phi."""
    gap = ibd - two_phi
    return {"mean": float(gap.mean()), "se": float(gap.std(ddof=1) / np.sqrt(len(gap)))}


def seed_g6(out: Path, data: Data, per_class: int, threads: int) -> dict:
    """Return the G6 record of one finished seed directory."""
    manifest = json.loads((out / "manifest.json").read_text())
    case = load_case(manifest["case"]["id"])
    seed = int(manifest["seed"])
    archs = [a["name"] for a in case["architectures"]]
    first = out / archs[0]
    ped = read_pedigree(first)
    first_pairs = sample_pairs(first, ped)
    sampled = subsample(first_pairs, per_class, seed)
    ids = np.concatenate(list(sampled.values()))

    graph = PedigreeGraph.from_arrays(
        ids=ped["id"].to_numpy(), mother_ids=ped["mother"].to_numpy(), father_ids=ped["father"].to_numpy()
    )
    two_phi = 2 * np.asarray(graph.pair_kinship(ids[:, 0], ids[:, 1]), dtype=np.float64)

    maps = load_maps(data)
    genome = Genome(
        chromosomes=len(CHROMOSOMES), morgans=mean_morgans(maps), causal_per_chromosome=CAUSAL_PER_CHROMOSOME
    )
    t0 = time.time()
    G, ibd_l2 = drop(ped, ids, genome, float(case["simulate"]["A1"]), stage_seed(seed, STAGE_G6_DROP))
    t_l2 = time.time() - t0
    t0 = time.time()
    tables = tskit.TableCollection.load(str(out / "pedigree_tables.trees"))
    ibd_ts = tskit_ibd(tables, maps, seed, ids, threads)
    t_ts = time.time() - t0

    classes, start = {}, 0
    for k, rows in sampled.items():
        s = slice(start, start + len(rows))
        start += len(rows)
        classes[k] = {
            "n": len(rows),
            "two_phi_minus_nominal": float(two_phi[s].mean() - TWO_PHI[k]),
            "ibd_l2": gap_summary(ibd_l2[s], two_phi[s]),
            "ibd_ts": gap_summary(ibd_ts[s], two_phi[s]),
        }
    rA = {}
    for name in archs:
        if name == archs[0]:
            arch_ped, pairs = ped, first_pairs
        else:
            arch_ped = read_pedigree(out / name)
            pairs = sample_pairs(out / name, arch_ped)
        a = arch_ped["A1"].to_numpy().astype(np.float64)
        rA[name] = {
            k: {"l2": num(pearson(G[i], G[j]) - TWO_PHI[k]), "ts": num(pearson(a[i], a[j]) - TWO_PHI[k]), "n": len(i)}
            for k, (i, j) in pairs.items()
        }
    return {
        "g6_version": G6_VERSION,
        "g6_by": validate_commit(),
        "seed": seed,
        "case": case["id"],
        "simace": manifest["simace"],
        "part_sha256": hashlib.sha256((out / "part.json").read_bytes()).hexdigest(),
        "genome": {"chromosomes": genome.chromosomes, "morgans": genome.morgans,
                   "causal_per_chromosome": genome.causal_per_chromosome},
        "pairs_per_class": per_class,
        "pairs_sha256": {k: hashlib.sha256(rows.astype("<i8").tobytes()).hexdigest() for k, rows in sampled.items()},
        "seconds": {"layer2_drop": t_l2, "tskit_drop": t_ts},
        "ibd": classes,
        "rA": rA,
    }  # fmt: skip


def main(argv: list[str] | None = None) -> int:
    """Write ``g6.json`` for one seed directory."""
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("seed_dir", type=Path)
    ap.add_argument("--data", type=Path, required=True)
    ap.add_argument("--pairs", type=int, default=2000, help="pairs per class for the IBD comparison")
    ap.add_argument("--threads", type=int, default=4, help="tskit chromosomes dropped at once")
    args = ap.parse_args(argv)
    record = seed_g6(args.seed_dir, Data(args.data), args.pairs, args.threads)
    (args.seed_dir / "g6.json").write_text(json.dumps(record, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
