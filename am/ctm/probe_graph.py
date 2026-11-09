"""Probe for rwaples/simACE-validate#5, open question 1: can xftsim build simACE's C transmission?

simACE (ADR 0025): c_child = b_m c_mother + b_f c_father + u, u one draw per mother with
variance C - S (S the transmitted part's realised variance). The fixed-residual analogue:
u = sqrt(C (1 - b_m^2 - b_f^2)) hh_mother, hh a private N(0, 1) trait never summed into y.

Two constructions, both random mating, y = A + c + E:
  issue: c = Sum(c.vertical, c.household); y reads c's proband phenotype through a
         LinearTransformationComponent (the issue's candidate), then Sum(y).
  dual:  one LinearVerticalComponent reading (c mother, c father, hh mother) writes the same
         value into c.vertical and y.vertical; one SumAllTransformation.
Each prints Var(c) by generation and the c and y correlations of FS, MO and FO, against the
expected values under random mating: Var(c) = C; c: MO b_m, FO b_f, FS 1;
y: MO A/2 + b_m C, FO A/2 + b_f C, FS A/2 + C.

Usage (xftsim env, repo root): python -m am.ctm.probe_graph {issue|dual} SEED N M G BM BF C
"""

import json
import sys
import warnings

warnings.filterwarnings("ignore")

import numpy as np  # noqa: E402
import xftsim as xft  # noqa: E402

from am.stats import _corr, _sib_pairs, _sym_corr  # noqa: E402

variant = sys.argv[1]
seed, n, m, g = (int(a) for a in sys.argv[2:6])
bm, bf, c = (float(a) for a in sys.argv[6:9])
h2, u = 0.5, c * (1 - bm**2 - bf**2)
np.random.seed(seed)

haps = xft.founders.founder_haplotypes_uniform_AFs(n=n, m=m)
vind = haps.xft.get_variant_indexer()
beta = np.random.normal(size=(m, 1))
beta *= np.sqrt(h2 / np.mean(beta**2))
effects = xft.effect.AdditiveEffects(
    beta=beta,
    variant_indexer=vind,
    component_indexer=xft.index.ComponentIndex.from_product(["y"], ["additiveGenetic"]),
    scaled=True,
    standardized=True,
)
v_in = xft.index.ComponentIndex(["c", "c", "hh"], ["phenotype"] * 3, [0, 1, 0])
v_in.comp_type = "output"
noise = xft.arch.AdditiveNoiseComponent([1 - h2 - c, 1.0], phenotype_name=["y", "hh"])
genetic = xft.arch.AdditiveGeneticComponent(effects)
coef = [bm, bf, np.sqrt(u)]
if variant == "dual":
    components = [
        genetic,
        noise,
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
else:
    c_in = xft.index.ComponentIndex(["c"], ["phenotype"], [-1])
    c_in.comp_type = "outcome"  # as SumTransformation labels it; "output" duplicates the column
    components = [
        genetic,
        noise,
        xft.arch.LinearVerticalComponent(
            input_cindex=v_in,
            output_cindex=xft.index.ComponentIndex(["c"], ["vertical"], [-1]),
            founder_variances=[c, c, 1.0],
            coefficient_matrix=np.array([coef]),
            normalize=False,
        ),
        xft.arch.SumTransformation(
            xft.index.ComponentIndex(["c", "hh"], ["vertical", "additiveNoise"]),
            xft.index.ComponentIndex(["c", "hh"], ["phenotype", "phenotype"]),
        ),
        xft.arch.LinearTransformationComponent(
            input_cindex=c_in,
            output_cindex=xft.index.ComponentIndex(["y"], ["cfrom"], [-1]),
            coefficient_matrix=np.array([[1.0]]),
            normalize=False,
        ),
        xft.arch.SumTransformation(
            xft.index.ComponentIndex(["y", "y", "y"], ["additiveGenetic", "additiveNoise", "cfrom"]),
            xft.index.ComponentIndex(["y"], ["phenotype"]),
        ),
    ]
arch = xft.arch.Architecture(components=components)
sim = xft.sim.Simulation(
    founder_haplotypes=haps,
    mating_regime=xft.mate.RandomMatingRegime(mates_per_female=1, offspring_per_pair=2),
    recombination_map=xft.reproduce.RecombinationMap.constant_map_from_haplotypes(haps, p=0.5),
    architecture=arch,
    statistics=[],
    post_processors=[xft.proc.LimitMemory(n_haplotype_generations=1)],
)
sim.run(g + 1)


def frame(gen):
    """Iids and {name: column} of the proband components of generation gen."""
    ph = sim.phenotype_store[gen]
    keys = zip(ph.phenotype_name.values, ph.component_name.values, ph.vorigin_relative.values)  # noqa: B905  (Python 3.9)
    cols = {f"{p}.{k}.{v}": np.asarray(ph.data[:, i]) for i, (p, k, v) in enumerate(keys)}
    return np.asarray(ph.iid.values), cols


_, f0 = frame(0)
out = [{"generation": 0, "var_c": float(np.var(f0["c.phenotype.-1"], ddof=1))}]
for gen in range(1, g + 1):
    piid, pc = frame(gen - 1)
    iid, cc = frame(gen)
    if gen == 1:
        y_resid = cc["y.phenotype.-1"] - cc["y.additiveGenetic.-1"] - cc["y.additiveNoise.-1"]
        out[0]["max_abs_y_minus_A_E_c_g1"] = float(np.max(np.abs(y_resid - cc["c.phenotype.-1"])))
    rep = sim.mating_store[gen - 1].get_reproduction_frame().set_index("iid")
    pos = {k: i for i, k in enumerate(piid)}
    mo = np.array([pos[k] for k in rep.loc[iid, "maternal_iid"]])
    fa = np.array([pos[k] for k in rep.loc[iid, "paternal_iid"]])
    _, couple = np.unique(mo * (len(piid) + 1) + fa, return_inverse=True)
    fi, fj = _sib_pairs(couple)
    row = {"generation": gen, "var_c": float(np.var(cc["c.phenotype.-1"], ddof=1))}
    for name in ("c", "y"):
        x, px = cc[f"{name}.phenotype.-1"], pc[f"{name}.phenotype.-1"]
        row[f"r_FS_{name}"] = _sym_corr(x, fi, fj)
        row[f"r_MO_{name}"] = _corr(px[mo], x)
        row[f"r_FO_{name}"] = _corr(px[fa], x)
    # the transmitted inputs must be the parents' own c
    row["max_abs_input_vs_mother_c"] = float(np.max(np.abs(cc["c.phenotype.0"] - pc["c.phenotype.-1"][mo])))
    out.append(row)
expected = {
    "var_c": c,
    "r_FS_c": 1.0,
    "r_MO_c": bm,
    "r_FO_c": bf,
    "r_FS_y": h2 / 2 + c,
    "r_MO_y": h2 / 2 + bm * c,
    "r_FO_y": h2 / 2 + bf * c,
}
print(json.dumps({"variant": variant, "seed": seed, "n": n, "bm": bm, "bf": bf, "C": c, "expected": expected}))
for row in out:
    print(json.dumps({k: round(v, 4) if isinstance(v, float) else v for k, v in row.items()}))
