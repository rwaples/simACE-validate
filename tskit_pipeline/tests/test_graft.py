"""The ancestry graft keeps every founder-haplotype mutation a descendant inherits."""

import msprime
import numpy as np
import pytest
import tskit
from genotype_drop_chrom import graft_ancestry

N_FOUNDERS = 4


def _ancestry() -> tskit.TreeSequence:
    ts = msprime.sim_ancestry(
        N_FOUNDERS, sequence_length=1e5, recombination_rate=1e-8, population_size=1e4, random_seed=1
    )
    return msprime.sim_mutations(ts, rate=1e-7, random_seed=2)


def _drop(seqlen: float) -> tskit.TreeSequence:
    """2N founder nodes at time 1 and one sample child inheriting founder node 0 everywhere."""
    tables = tskit.TableCollection(seqlen)
    for _ in range(2 * N_FOUNDERS):
        tables.nodes.add_row(flags=0, time=1)
    child = tables.nodes.add_row(flags=tskit.NODE_IS_SAMPLE, time=0)
    tables.edges.add_row(0, seqlen, 0, child)
    return tables.tree_sequence()


@pytest.fixture(scope="module")
def grafted():
    anc = _ancestry()
    out = graft_ancestry(_drop(anc.sequence_length), anc, N_FOUNDERS, 1, np.arange(N_FOUNDERS))
    return anc, out


def test_mutations_on_founder_nodes_survive(grafted):
    anc, out = grafted
    on_founders = np.isin(anc.tables.mutations.node, np.arange(2 * N_FOUNDERS)).sum()
    assert on_founders > 0
    assert out.num_mutations == anc.num_mutations


def test_child_carries_founder_node_0(grafted):
    anc, out = grafted
    child = out.samples()[-1]
    assert out.num_sites == anc.num_sites
    child_g = out.genotype_matrix(samples=[child])[:, 0]
    founder_g = anc.genotype_matrix(samples=[0])[:, 0]
    assert np.array_equal(child_g, founder_g)
