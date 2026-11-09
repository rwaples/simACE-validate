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

BM, BF > 0 (C > 0, one trait) add simACE's vertical C transmission (ADR 0025 in
simACE; rwaples/simACE-validate#5): a phenotype c = BM c_mother + BF c_father +
sqrt(C (1 - BM^2 - BF^2)) hh_mother, one LinearVerticalComponent writing that value into
both c and y (normalize=False). The household term has the fixed random-mating residual
variance, where simACE tops up from the realised covariance (differences.md D5). Founders'
parental c inputs are drawn N(0, C), so founder Var(c) is C.

PAIRING halfsib (r = 0) replaces LinearAssortativeMatingRegime by
RandomHalfSibMatingRegime, which gives some individuals two or more partners, so the
pedigree has maternal and paternal half-sibs.

Usage: python -m am.run_xft SEED N M G R [TRAITS [MATE [C [BM [BF [PAIRING]]]]]]
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

from am.stats import cross_trait_stats, ctm_stats, generation_stats  # noqa: E402

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
bm = float(sys.argv[9]) if len(sys.argv) > 9 else 0.0
bf = float(sys.argv[10]) if len(sys.argv) > 10 else 0.0
pairing = sys.argv[11] if len(sys.argv) > 11 else "lam"
ctm = bm > 0 or bf > 0
if ctm and not (c > 0 and bm**2 + bf**2 < 1):
    raise SystemExit("BM, BF need C > 0 and BM^2 + BF^2 < 1")
if pairing not in ("lam", "halfsib") or (pairing == "halfsib" and (r1 != 0 or traits != 1)):
    raise SystemExit("PAIRING is lam or halfsib; halfsib needs R 0 and TRAITS 1")
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
    """Assortative mating on two components towards the rank-1 mate correlation v v' (ADR 0003, decision 2).

    Built on xftsim's MatingRegime API: get_potential_mates gives the two sides, enumerate_assignment the
    couples. Each generation, on each side: Z is the two components standardised, S its covariance, and the
    score s = Z S^-1 v has Cov(Z, s) = v and Var(s) = q = v' S^-1 v. With rho = sqrt(q_f q_m), both sides'
    noisy scores sqrt(rho) s + sqrt(1 - rho) SD(s) e, e ~ N(0, 1), are ranked and the k-th of one side
    pairs with the k-th of the other, so corr(s_f, s_m) = rho and corr(Z_f, Z_m) = v v'. A rho outside
    (0, 1) raises ValueError rather than being clamped. A larger side gives up a random excess.
    """

    def __init__(self, component_index, v, mates_per_female, offspring_per_pair):
        """Mate on the two components of component_index; v = (sqrt(R1), sqrt(R2))."""
        if component_index.k_total != 2:
            raise ValueError("the rank-1 regime mates on exactly two components")
        super().__init__(
            offspring_per_pair=xft.utils.ConstantCount(offspring_per_pair),
            mates_per_female=xft.utils.ConstantCount(mates_per_female),
            female_offspring_per_pair="balanced",
            sex_aware=False,
            exhaustive=True,
            component_index=component_index,
            haplotypes=False,
        )
        self.component_index = component_index
        self.v = np.asarray(v, dtype=np.float64)

    def _scores(self, x):
        """Score Z S^-1 v of an n x 2 side, and its variance q."""
        z = (x - x.mean(axis=0)) / x.std(axis=0)
        s = z @ np.linalg.solve(np.cov(z, rowvar=False, ddof=0), self.v)
        return s, float(np.var(s))

    def mate(self, haplotypes=None, phenotypes=None, control=None):
        """Pair the two sides by rank of their noisy scores and return xftsim's MateAssignment."""
        females, males = self.get_potential_mates(haplotypes, phenotypes)
        k = min(len(females), len(males))
        if len(females) > k:
            females = np.random.permutation(females)[:k]
        if len(males) > k:
            males = np.random.permutation(males)[:k]
        y = np.asarray(phenotypes.xft[None, self.component_index].data, dtype=np.float64)
        row = {label: i for i, label in enumerate(phenotypes.sample.values)}
        s_f, q_f = self._scores(y[[row[label] for label in females]])
        s_m, q_m = self._scores(y[[row[label] for label in males]])
        rho = np.sqrt(q_f * q_m)
        if not 0 < rho < 1:
            raise ValueError(f"rank-1 target infeasible this generation: rho = sqrt(q_f q_m) = {rho:.4f}")

        def noisy(s, q):
            return np.sqrt(rho) * s + np.sqrt((1 - rho) * q) * np.random.normal(size=len(s))

        pairs_f = females[np.argsort(noisy(s_f, q_f))]
        pairs_m = males[np.argsort(noisy(s_m, q_m))]
        return self.enumerate_assignment(pairs_f, pairs_m, haplotypes, phenotypes)


class RandomHalfSibMatingRegime(xft.mate.MatingRegime):
    """Random mating that makes half-sibs: n couples drawn with replacement, each with two offspring.

    Built on xftsim's MatingRegime API (get_potential_mates, enumerate_assignment), written for
    rwaples/simACE-validate#5 without reading any xftsim regime's source. Each of the n couples
    (n = the smaller side) takes a uniform random member of each side, independently and with
    replacement, so partner counts are about Poisson(1) and a member with two or more partners
    has half-sib offspring. A couple drawn twice is redrawn, so couples are distinct; constant N
    follows from two offspring per couple.
    """

    def __init__(self):
        """Two offspring per couple; enumerate_assignment receives each couple once."""
        super().__init__(
            offspring_per_pair=xft.utils.ConstantCount(2),
            mates_per_female=xft.utils.ConstantCount(1),
            female_offspring_per_pair="balanced",
            sex_aware=False,
            exhaustive=True,
            component_index=None,
            haplotypes=False,
        )

    def mate(self, haplotypes=None, phenotypes=None, control=None):
        """Draw the couples and return xftsim's MateAssignment."""
        females, males = self.get_potential_mates(haplotypes, phenotypes)
        n_f, n_m = len(females), len(males)
        k = min(n_f, n_m)
        f = np.random.randint(n_f, size=k)
        m = np.random.randint(n_m, size=k)
        while True:
            _, first = np.unique(f * n_m + m, return_index=True)
            dup = np.setdiff1d(np.arange(k), first)
            if not len(dup):
                break
            f[dup] = np.random.randint(n_f, size=len(dup))
            m[dup] = np.random.randint(n_m, size=len(dup))
        return self.enumerate_assignment(females[f], males[m], haplotypes, phenotypes)


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
elif ctm:
    v_in = xft.index.ComponentIndex(
        ["c", "c", "hh"], ["phenotype"] * 3, [0, 1, 0]
    )  # mother's c, father's c, mother's hh
    v_in.comp_type = "output"
    coef = [bm, bf, np.sqrt(c * (1 - bm**2 - bf**2))]
    arch = xft.arch.Architecture(
        components=[
            xft.arch.AdditiveGeneticComponent(effects),
            xft.arch.AdditiveNoiseComponent([1 - h2 - c, 1.0], phenotype_name=["y", "hh"]),
            xft.arch.LinearVerticalComponent(
                input_cindex=v_in,
                output_cindex=xft.index.ComponentIndex(["c", "y"], ["vertical", "vertical"], [-1, -1]),
                founder_variances=[c, c, 1.0],
                coefficient_matrix=np.array([coef, coef]),
                normalize=False,
            ),
            xft.arch.SumAllTransformation(
                xft.index.ComponentIndex(
                    ["y", "y", "y", "c", "hh"],
                    ["additiveGenetic", "additiveNoise", "vertical", "vertical", "additiveNoise"],
                )
            ),
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
if pairing == "halfsib":
    mating = RandomHalfSibMatingRegime()
elif rank1:
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


def c_of(gen):
    """The c phenotype of generation gen, in xftsim's sample order (C transmission only)."""
    ph = sim.phenotype_store[gen]
    keys = zip(ph.phenotype_name.values, ph.component_name.values, ph.vorigin_relative.values)  # noqa: B905  (Python 3.9)
    return np.asarray(ph.data[:, list(keys).index(("c", "phenotype", -1))])


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
    if ctm:
        row.update(ctm_stats(c_of(gen), v[0][1], mother, father, c_of(gen - 1), pv[0][1]))
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
            "bm": bm,
            "bf": bf,
            "pairing": pairing,
            "h2": h2,
            "wall_s": round(time.perf_counter() - t0, 2),
            "maxrss_mb": round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, 1),
            "rows": rows,
        }
    )
)
