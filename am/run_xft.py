"""One xftsim run for the assortative-mating pilot (issue #4); prints one JSON line.

Runs in the `xftsim` pixi env (Python 3.9). Trait: y = A + E with Var(A) = h2
and Var(E) = 1 - h2 in the founders; M unlinked biallelic loci (xftsim's
constant map, p = 0.5 between adjacent loci). Effects are drawn N(0, 1) and
rescaled so that mean(beta^2) = h2 exactly, so the founders' expected Var(A)
is h2 given the allele frequencies (GCTA_Architecture leaves it random, with
relative SD sqrt(2/M)). Mating is LinearAssortativeMatingRegime on y with
correlation r (r = 0 is random pairing by the same code path), one mate per
individual and two offspring per couple, so N is constant.

Usage: run_xft.py SEED N M G R
Generation g in the output is the offspring of g rounds of mating (g = 0 are
the founders).
"""

import json
import os
import resource
import sys
import time
import warnings

warnings.filterwarnings("ignore")  # pandas 2.x FutureWarnings from xftsim/index.py

import numpy as np  # noqa: E402
import xftsim as xft  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numba as nb  # noqa: E402
from stats import generation_stats  # noqa: E402

seed, n, m, g = (int(a) for a in sys.argv[1:5])
r = float(sys.argv[5])
h2 = 0.5
np.random.seed(seed)


# xftsim's meiosis draws inside numba (xftsim/reproduce.py, _meiosis_i under a
# parallel prange), whose RNG state is separate from numpy's and per thread.
# Seed it on the main thread and in the parallel region. Reproducible only with
# NUMBA_NUM_THREADS=1, which am/pilot.sh sets.
@nb.njit
def _seed_main(s):
    np.random.seed(s)


@nb.njit(parallel=True)
def _seed_threads(s):
    for i in nb.prange(nb.get_num_threads()):
        np.random.seed(s + 7919 * (i + 1))


_seed_main(seed)
_seed_threads(seed)
xft.config.print_durations_threshold = np.inf

t0 = time.perf_counter()
haps = xft.founders.founder_haplotypes_uniform_AFs(n=n, m=m)
vind = haps.xft.get_variant_indexer()
beta = np.random.normal(size=(m, 1))
beta *= np.sqrt(h2 / np.mean(beta**2))
cind = xft.index.ComponentIndex.from_product(["y"], ["additiveGenetic"])
effects = xft.effect.AdditiveEffects(
    beta=beta, variant_indexer=vind, component_indexer=cind, scaled=True, standardized=True
)
iind, _ = xft.arch.SumTransformation.construct_cindexes(["y"])
arch = xft.arch.Architecture(
    components=[
        xft.arch.AdditiveGeneticComponent(effects),
        xft.arch.AdditiveNoiseComponent([1 - h2], phenotype_name=["y"]),
        xft.arch.SumAllTransformation(iind),
    ]
)
rmap = xft.reproduce.RecombinationMap.constant_map_from_haplotypes(haps, p=0.5)
mating = xft.mate.LinearAssortativeMatingRegime(
    component_index=xft.index.ComponentIndex.from_product(["y"], ["phenotype"]),
    r=r,
    mates_per_female=1,
    offspring_per_pair=2,
)
sim = xft.sim.Simulation(
    founder_haplotypes=haps,
    mating_regime=mating,
    recombination_map=rmap,
    architecture=arch,
    statistics=[],
    post_processors=[xft.proc.LimitMemory(n_haplotype_generations=1)],
)
sim.run(g + 1)  # generations 0..g


def frame(gen):
    """Iids, A and y of generation gen, in xftsim's sample order."""
    ph = sim.phenotype_store[gen]
    comp = list(ph.component_name.values)
    a = np.asarray(ph.data[:, comp.index("additiveGenetic")])
    y = np.asarray(ph.data[:, comp.index("phenotype")])
    return np.asarray(ph.iid.values), a, y


rows = []
iid0, a0, y0 = frame(0)
rows.append({"generation": 0, "var_A": float(np.var(a0, ddof=1)), "var_L": float(np.var(y0, ddof=1))})
for gen in range(1, g + 1):
    piid, pa, py = frame(gen - 1)
    iid, a, y = frame(gen)
    rep = sim.mating_store[gen - 1].get_reproduction_frame()
    pos = {k: i for i, k in enumerate(piid)}
    by_iid = rep.set_index("iid")
    mother = np.array([pos[k] for k in by_iid.loc[iid, "maternal_iid"]])
    father = np.array([pos[k] for k in by_iid.loc[iid, "paternal_iid"]])
    rows.append({"generation": gen, **generation_stats(a, y, mother, father, pa, py)})

print(
    json.dumps(
        {
            "tool": "xftsim",
            "version": xft.__version__,
            "seed": seed,
            "n": n,
            "m": m,
            "g": g,
            "r": r,
            "h2": h2,
            "wall_s": round(time.perf_counter() - t0, 2),
            "maxrss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
            "rows": rows,
        }
    )
)
