"""One xftsim run for the assortative-mating pilot (issue #4); prints one JSON line.

Runs in the `xftsim` pixi env (Python 3.9). Trait: y = A + E with Var(A) = h2
and Var(E) = 1 - h2 in the founders; M unlinked biallelic loci (xftsim's
constant map, p = 0.5 between adjacent loci). Effects are drawn N(0, 1) and
rescaled so that mean(beta^2) = h2 exactly, so the founders' expected Var(A)
is h2 given the allele frequencies (GCTA_Architecture leaves it random, with
relative SD sqrt(2/M)). Mating is LinearAssortativeMatingRegime on y with
correlation r (r = 0 is random pairing by the same code path), one mate per
individual and two offspring per couple, so N is constant.

With TRAITS = 2, two independent traits (independent effects and noise) each
have h2 = 0.5, and the regime mates on both phenotypes, which pairs on their
standardised sum: every mate correlation (y1-y1, y2-y2, y1-y2, y2-y1) is r.

With TRAITS = 2, R may be R1/R2 (both > 0): mate correlations R1 on trait 1,
R2 on trait 2 and sqrt(R1 R2) across traits, the rank-1 law xftsim's LAM
cannot reach (its general regime needs the commercial LocalSolver).
Rank1AssortativeMatingRegime below pins it (ADR 0003, decision 2); the
record's r is then the string "R1/R2".

MATE picks what the regime mates on (power mutations, ADR 0003 decision 8):
phenotype (default), additive (the additive genetic component instead of the
phenotype), or y1 (two traits, but mating on trait 1 alone).

C > 0 (one trait only) adds a simACE-style household C, built by the
xftsim-bench session: a private noise trait hh (variance 1, never summed into
y) and a LinearVerticalComponent that passes the mother's hh, normalised and
scaled by sqrt(C), into the child's y. Children of one mother share that
draw, and it is independent of the parents' own y, as simACE's household C
is (household = mother). E is then 1 - h2 - C.

Usage: python -m am.run_xft SEED N M G R [TRAITS [MATE [C]]]
Generation g in the output is the offspring of g rounds of mating (g = 0 are
the founders).
"""

import json
import resource
import sys
import time
import warnings

warnings.filterwarnings("ignore")  # pandas 2.x FutureWarnings from xftsim/index.py

import numba as nb  # noqa: E402
import numpy as np  # noqa: E402
import xftsim as xft  # noqa: E402

from am.stats import cross_trait_stats, generation_stats  # noqa: E402

seed, n, m, g = (int(a) for a in sys.argv[1:5])
traits = int(sys.argv[6]) if len(sys.argv) > 6 else 1
r1, r2 = (float(x) for x in sys.argv[5].split("/")) if "/" in sys.argv[5] else (float(sys.argv[5]),) * 2
rank1 = r1 != r2
if rank1 and (traits != 2 or min(r1, r2) <= 0):
    raise SystemExit("R1/R2 needs TRAITS 2 and R1, R2 > 0")
r = sys.argv[5] if rank1 else r1
names = ["y"] if traits == 1 else ["y1", "y2"]
mate = sys.argv[7] if len(sys.argv) > 7 else "phenotype"
mate_index = {
    "phenotype": (names, "phenotype"),
    "additive": (names, "additiveGenetic"),
    "y1": (["y1"], "phenotype"),
}[mate]
if mate == "y1" and (traits != 2 or rank1):
    raise SystemExit("MATE y1 needs TRAITS 2 and one R")
c = float(sys.argv[8]) if len(sys.argv) > 8 else 0.0
if c and traits != 1:
    raise SystemExit("C > 0 needs TRAITS 1")
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


class Rank1AssortativeMatingRegime(xft.mate.MatingRegime):
    """LinearAssortativeMatingRegime with a weighted score, pinning the mate correlation v v'.

    v = (sqrt(R1), sqrt(R2)). Each side's score is s = Z w with w = S^-1 v (Z the
    standardised traits, S their covariance on that side), so Cov(Z, s) = v.
    Pairing the noisy scores sqrt(rho) s + sqrt(1 - rho) sd(s) e by rank gives
    corr(s_f, s_m) = rho, hence corr(Z_f, Z_m) = rho v v' / sqrt(q_f q_m) with
    q = v' S^-1 v; rho = sqrt(q_f q_m) pins it at v v'. At R1 = R2 it is LAM's law
    (w proportional to (1, 1), rho LAM's R) and reproduces LAM's output exactly.
    Mirrors LAM.mate except the score, and raises where LAM would clamp.
    """

    def __init__(self, component_index, v, mates_per_female, offspring_per_pair):
        """Mate on the two components of component_index towards the target v v'."""
        super().__init__(
            self,
            offspring_per_pair=xft.utils.ConstantCount(offspring_per_pair),
            mates_per_female=xft.utils.ConstantCount(mates_per_female),
            female_offspring_per_pair="balanced",
            sex_aware=False,
            exhaustive=True,
            component_index=component_index,
            haplotypes=False,
        )
        assert component_index.k_total == 2, "the rank-1 regime mates on two traits"
        self.component_index = component_index
        self.v = np.asarray(v)

    def mate(self, haplotypes=None, phenotypes=None, control=None):
        """Pair females and males on their weighted scores."""
        female_indices, male_indices = self.get_potential_mates(haplotypes, phenotypes)
        sample_indexer = phenotypes.xft.get_sample_indexer()
        female_indexer = sample_indexer[female_indices]
        male_indexer = sample_indexer[male_indices]
        n_new = np.min([female_indexer.n, male_indexer.n])
        female_indexer = female_indexer.at_most(n_new)
        male_indexer = male_indexer.at_most(n_new)
        zf = xft.utils.standardize_array(phenotypes.xft[female_indexer, self.component_index].data)
        zm = xft.utils.standardize_array(phenotypes.xft[male_indexer, self.component_index].data)
        n_pairs = zm.shape[0]
        wf = np.linalg.solve(np.cov(zf, rowvar=False), self.v)
        wm = np.linalg.solve(np.cov(zm, rowvar=False), self.v)
        rho = np.sqrt((self.v @ wf) * (self.v @ wm))
        if not 0 < rho < 1:
            raise ValueError(f"rank-1 target infeasible at this generation: rho {rho}")
        sf, sm = zf @ wf, zm @ wm
        score_f = sf * np.sqrt(rho) + np.sqrt(1 - rho) * np.random.normal(0, np.std(sf), n_pairs)
        score_m = sm * np.sqrt(rho) + np.sqrt(1 - rho) * np.random.normal(0, np.std(sm), n_pairs)
        return self.enumerate_assignment(
            female_indices=female_indices[np.argsort(score_f)[:n_new]],
            male_indices=male_indices[np.argsort(score_m)[:n_new]],
            haplotypes=haplotypes,
            phenotypes=phenotypes,
        )


t0 = time.perf_counter()
haps = xft.founders.founder_haplotypes_uniform_AFs(n=n, m=m)
vind = haps.xft.get_variant_indexer()
beta = np.random.normal(size=(m, traits))
beta *= np.sqrt(h2 / np.mean(beta**2, axis=0))
cind = xft.index.ComponentIndex.from_product(names, ["additiveGenetic"])
effects = xft.effect.AdditiveEffects(
    beta=beta, variant_indexer=vind, component_indexer=cind, scaled=True, standardized=True
)
if c == 0:
    iind, _ = xft.arch.SumTransformation.construct_cindexes(names)
    arch = xft.arch.Architecture(
        components=[
            xft.arch.AdditiveGeneticComponent(effects),
            xft.arch.AdditiveNoiseComponent([1 - h2] * traits, phenotype_name=names),
            xft.arch.SumAllTransformation(iind),
        ]
    )
else:
    v_in = xft.index.ComponentIndex.from_product(["hh"], ["phenotype"], [0])  # vorigin 0: the mother
    v_in.comp_type = "output"
    arch = xft.arch.Architecture(
        components=[
            xft.arch.AdditiveGeneticComponent(effects),
            xft.arch.AdditiveNoiseComponent([1 - h2 - c, 1.0], phenotype_name=["y", "hh"]),
            xft.arch.LinearVerticalComponent(
                input_cindex=v_in,
                output_cindex=xft.index.ComponentIndex.from_product(["y"], ["vertical"], [-1]),
                founder_variances=[1.0],
                coefficient_matrix=np.array([[np.sqrt(c)]]),
                normalize=True,
            ),
            xft.arch.SumAllTransformation(
                xft.index.ComponentIndex(
                    ["y", "y", "y", "hh"], ["additiveGenetic", "additiveNoise", "vertical", "additiveNoise"]
                )
            ),
        ]
    )
rmap = xft.reproduce.RecombinationMap.constant_map_from_haplotypes(haps, p=0.5)
mate_cindex = xft.index.ComponentIndex.from_product(mate_index[0], [mate_index[1]])
if rank1:
    mating = Rank1AssortativeMatingRegime(mate_cindex, np.sqrt([r1, r2]), mates_per_female=1, offspring_per_pair=2)
else:
    mating = xft.mate.LinearAssortativeMatingRegime(
        component_index=mate_cindex, r=r1, mates_per_female=1, offspring_per_pair=2
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
    """Iids, then A and y of each trait, of generation gen, in xftsim's sample order."""
    ph = sim.phenotype_store[gen]
    cname = list(ph.component_name.values)
    pname = list(ph.phenotype_name.values)

    def col(p, c):
        return np.asarray(ph.data[:, next(i for i in range(len(cname)) if cname[i] == c and pname[i] == p)])

    return np.asarray(ph.iid.values), [(col(p, "additiveGenetic"), col(p, "phenotype")) for p in names]


def suffixed(d, t):
    """Trait 1 keys unchanged; trait 2 keys get `_t2`."""
    return d if t == 0 else {f"{k}_t{t + 1}": v for k, v in d.items()}


rows = []
_, v0 = frame(0)
row0 = {"generation": 0}
for t, (a0, y0) in enumerate(v0):
    row0.update(suffixed({"var_A": float(np.var(a0, ddof=1)), "var_L": float(np.var(y0, ddof=1))}, t))
rows.append(row0)
for gen in range(1, g + 1):
    piid, pv = frame(gen - 1)
    iid, v = frame(gen)
    rep = sim.mating_store[gen - 1].get_reproduction_frame()
    pos = {k: i for i, k in enumerate(piid)}
    by_iid = rep.set_index("iid")
    mother = np.array([pos[k] for k in by_iid.loc[iid, "maternal_iid"]])
    father = np.array([pos[k] for k in by_iid.loc[iid, "paternal_iid"]])
    row = {"generation": gen}
    for t in range(traits):
        row.update(suffixed(generation_stats(v[t][0], v[t][1], mother, father, pv[t][0], pv[t][1]), t))
    if traits == 2:
        row.update(cross_trait_stats(v, pv, mother, father))
    rows.append(row)

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
            "traits": traits,
            "mate": mate,
            "C": c,
            "h2": h2,
            "wall_s": round(time.perf_counter() - t0, 2),
            "maxrss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
            "rows": rows,
        }
    )
)
