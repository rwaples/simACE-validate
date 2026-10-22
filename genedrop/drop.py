"""Layer-2 gene-drop through a simACE pedigree with msprime (fixed pedigree).

Drops L unlinked chromosomes through the recorded pedigree, and returns, per
individual, a finite-locus genetic value G and, per requested pair, the
realised IBD fraction. Recorded generation 0 are the founders, with distinct
genomes and no recapitation, matching the recorded-pedigree kinship, which
also treats them as unrelated (ADR 0001 §3).

MZ co-twins are one genetic individual. Only the lower-id twin enters
msprime, the co-twin's children are re-pointed to it, and the co-twin
copies its genome. This is unlike scripts/gene_drop in simACE, which drops
co-twins as separate offspring.
"""

from __future__ import annotations

from collections.abc import Callable, Iterable, Iterator
from dataclasses import dataclass
from typing import TYPE_CHECKING

import msprime
import numpy as np

if TYPE_CHECKING:
    import polars as pl

BP_PER_MORGAN = 100_000_000


@dataclass(frozen=True)
class Genome:
    """The dropped genome: L unlinked chromosomes of one map length each."""

    chromosomes: int
    morgans: float
    causal_per_chromosome: int


def keepers(ped: pl.DataFrame) -> np.ndarray:
    """Map each row to the individual msprime simulates (the lower-id MZ twin)."""
    idx = np.arange(ped.height)
    twin = ped["twin"].to_numpy()
    co = (twin >= 0) & (twin < idx)
    idx[co] = twin[co]
    return idx


def build_pedigree(ped: pl.DataFrame, keep: np.ndarray, sequence_length: float):
    """Return msprime pedigree tables for the keeper rows, and row -> msprime individual."""
    g_max = int(ped["generation"].max())
    gen = ped["generation"].to_numpy()
    mother = ped["mother"].to_numpy()
    father = ped["father"].to_numpy()
    pb = msprime.PedigreeBuilder()
    ms_of = np.full(ped.height, -1, dtype=np.int64)
    for i in range(ped.height):
        if keep[i] != i:
            continue
        parents = None if mother[i] < 0 else [int(ms_of[keep[mother[i]]]), int(ms_of[keep[father[i]]])]
        ms_of[i] = pb.add_individual(time=g_max - int(gen[i]), parents=parents, is_sample=True)
    return pb.finalise(sequence_length=sequence_length), ms_of[keep]


def roots_at(tree, nodes: np.ndarray) -> np.ndarray:
    """Return the root (a founder node) above each node in one tree."""
    parent = tree.parent_array
    out = nodes.copy()
    while True:
        up = parent[out]
        move = up >= 0
        if not move.any():
            return out
        out[move] = up[move]


def individual_nodes(ts) -> np.ndarray:
    """Return each individual's two nodes in id order (individuals x 2), as ts.individual(k).nodes."""
    owner = ts.nodes_individual
    if not np.all(np.bincount(owner[owner >= 0], minlength=ts.num_individuals) == 2):
        raise ValueError("every individual must own exactly two nodes")
    ordered = np.argsort(owner, kind="stable")[np.count_nonzero(owner < 0) :]
    return ordered.astype(np.int32).reshape(-1, 2)


def pair_spans(ts, nodes: np.ndarray, pairs: np.ndarray) -> np.ndarray:
    """Return the IBD span of each pair's four node pairs (pairs x 2 x 2), 0 where not IBD."""
    a = nodes[pairs[:, 0]][:, :, None]
    b = nodes[pairs[:, 1]][:, None, :]
    if np.any(a[:, 0, 0] == b[:, 0, 0]):
        raise ValueError("a pair is one genetic individual (an MZ pair); its IBD is not a pair statistic")
    lo, hi = np.broadcast_arrays(np.minimum(a, b), np.maximum(a, b))
    code = lo.astype(np.int64) * ts.num_nodes + hi
    ibd = ts.ibd_segments(store_pairs=True)
    keys = ibd.pairs.astype(np.int64)
    found = np.isin(code, keys[:, 0] * ts.num_nodes + keys[:, 1])
    span = np.zeros(code.shape)
    span[found] = [ibd[(int(u), int(v))].total_span for u, v in zip(lo[found], hi[found], strict=True)]
    return span


Mapper = Callable[[Callable[[int], object], Iterable[int]], Iterator]


def drop(
    ped: pl.DataFrame, pairs: np.ndarray, genome: Genome, A: float, seed: int, map_chromosomes: Mapper = map
) -> tuple[np.ndarray, np.ndarray]:
    """Gene-drop once. Return G per row and the IBD fraction per pair (2-phi scale).

    The IBD fraction of individuals i and j is the summed IBD span over their
    four node pairs divided by twice the genome length, whose expectation is
    2 * phi_ij.

    ``map_chromosomes`` runs the per-chromosome msprime simulations. The
    builtin map draws each chromosome's seed just before its simulation; a
    thread pool's map (msprime releases the GIL) draws every seed first. The
    two give the same values only without causal loci, whose draws otherwise
    fall between the seeds, so a drop with causal loci must use the builtin map.
    """
    if genome.causal_per_chromosome and map_chromosomes is not map:
        raise ValueError("a drop with causal loci simulates its chromosomes in order (map_chromosomes=map)")
    rng = np.random.default_rng(seed)
    keep = keepers(ped)
    length = genome.morgans * BP_PER_MORGAN
    tables, ms_of = build_pedigree(ped, keep, length)
    n_causal = genome.chromosomes * genome.causal_per_chromosome
    beta = rng.normal(size=n_causal)
    scale = np.sqrt(A / (np.sum(beta**2) / 2)) if n_causal else 0.0  # Var(G) = A at allele frequency 1/2
    dosage_value = np.zeros(ped.height)
    shared = np.zeros(len(pairs))
    locus = 0
    seeds = (int(rng.integers(1, 2**31 - 1)) for _ in range(genome.chromosomes))

    def simulate(chromosome_seed: int):
        return msprime.sim_ancestry(
            initial_state=tables,
            model="fixed_pedigree",
            recombination_rate=1 / BP_PER_MORGAN,
            random_seed=chromosome_seed,
        )

    for ts in map_chromosomes(simulate, seeds):
        nodes = individual_nodes(ts)[ms_of]  # rows x 2
        span = pair_spans(ts, nodes, pairs)
        for a in range(2):  # the four node pairs in a fixed order, so the sums are reproducible
            for b in range(2):
                shared += span[:, a, b]
        allele = np.zeros(ts.num_nodes, dtype=np.int64)
        for x in np.sort(rng.uniform(0, length, size=genome.causal_per_chromosome)):
            tree = ts.at(x)
            roots = roots_at(tree, nodes.ravel()).reshape(nodes.shape)
            is_root = np.zeros(ts.num_nodes, dtype=bool)
            is_root[roots] = True
            founder_nodes = np.flatnonzero(is_root)  # np.unique(roots), without its slow hash path
            allele[founder_nodes] = rng.integers(0, 2, size=len(founder_nodes))
            dosage = allele[roots].sum(axis=1)
            dosage_value += scale * beta[locus] * (dosage - 1)
            locus += 1
    return dosage_value, shared / (2 * length * genome.chromosomes)
