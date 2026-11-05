"""One xftsim run of two-trait AM with unequal mate correlations (rank 1); prints one JSON line.

Pilot for the unequal-r two-trait case (issue #4). xftsim's free regime,
LinearAssortativeMatingRegime, pairs on the standardised sum of the traits, so
every mate correlation equals one r. Its GeneralAssortativeMatingRegime takes
any K x K target but needs the commercial LocalSolver. This runs am/run_xft.py
unchanged (its sha256 pins the committed references) with LAM replaced by
Rank1AssortativeMatingRegime, which pins the rank-1 target R = v v' with
v = (sqrt(r1), sqrt(r2)): r1, r2 on the diagonal and sqrt(r1 r2) across traits.

Each side's score is s = Z w with w = S^-1 v (Z the standardised traits, S
their covariance on that side), so Cov(Z, s) = v. Pairing the noisy scores
sqrt(rho) s + sqrt(1 - rho) sd(s) e by rank gives corr(s_f, s_m) = rho, hence
corr(Z_f, Z_m) = rho v v' / sqrt(q_f q_m) with q = v' S^-1 v, and
rho = sqrt(q_f q_m) pins it at R. For r1 = r2 this is LAM's law exactly: w is
proportional to (1, 1) and rho is LAM's R.

Meiosis, effects, the architecture and the pedigree stay xftsim's; the pairing
law is this file's, written against xftsim's MatingRegime interface, and
mirrors LAM.mate except the score and an error, not LAM's clamp, when rho leaves (0, 1).

Usage: python -m am.xft_rank1 SEED N M G R1 R2
"""

import contextlib
import io
import json
import runpy
import sys
import warnings

warnings.filterwarnings("ignore")  # pandas 2.x FutureWarnings from xftsim/index.py

import numpy as np  # noqa: E402
import xftsim as xft  # noqa: E402

seed, n, m, g = sys.argv[1:5]
r1, r2 = float(sys.argv[5]), float(sys.argv[6])
label = f"{sys.argv[5]}/{sys.argv[6]}"  # the record's r, as am/run_simace.py writes it
if min(r1, r2) <= 0:
    raise SystemExit("R1, R2 > 0")
V = np.sqrt([r1, r2])
built = []


class Rank1AssortativeMatingRegime(xft.mate.MatingRegime):
    """LAM with the score weighted to pin the rank-1 mate correlation V V'."""

    def __init__(self, component_index, r=0, mates_per_female=1, offspring_per_pair=2):
        """Take run_xft.py's LAM arguments; r is unused (the target is V V')."""
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
        self.component_index = component_index
        assert component_index.k_total == 2, "rank-1 regime mates on two traits"
        built.append(self)

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
        wf = np.linalg.solve(np.cov(zf, rowvar=False), V)
        wm = np.linalg.solve(np.cov(zm, rowvar=False), V)
        rho = np.sqrt((V @ wf) * (V @ wm))
        if not 0 < rho < 1:  # LAM clamps instead; here that would silently miss the target
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


xft.mate.LinearAssortativeMatingRegime = Rank1AssortativeMatingRegime
sys.argv = ["am.run_xft", seed, n, m, g, str(r1), "2", "phenotype"]
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    runpy.run_module("am.run_xft", run_name="__main__")
assert built, "am/run_xft.py no longer builds its regime through xft.mate.LinearAssortativeMatingRegime"
(line,) = buf.getvalue().splitlines()
print(json.dumps(json.loads(line) | {"r": label, "mate": "rank1"}))
