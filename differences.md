# Known differences

Each entry explains a gap between simACE's law and an outside model. Rows
carrying the id are INFO: they are reported with the measured gap and are
never scored.

## D1. Segregation variance under inbreeding

simACE draws an offspring's A as the midparent plus N(0, A/2) whatever the
parents' inbreeding (simACE `simulate.py`, `generate_mendelian_noise` with
`sd_A/√2`). The infinitesimal model, for an additive trait at unlinked
loci, uses N(0, (A/2)(1 − (F_m + F_f)/2)), where F_m and F_f are the
parents' inbreeding coefficients. The two laws agree until parents are
inbred.

Case `deep_infinitesimal` (N = 100, 12 recorded generations) scores
simACE's pedigrees against the infinitesimal law. Its companion
`deep_constant` scores the same pedigrees against simACE's own law, so a
gap between the two cases isolates this model choice. The size of the gap
is the mean z of the `var_A*_g*` rows in `deep_infinitesimal`, by
generation. In the 2026-10-21 run, simACE's Var(A) sat above the
infinitesimal law in every generation: mean z +0.03 to +0.52, mostly
+0.2 to +0.5 from generation 5 on.

The F used is the recorded pedigree's. Recorded generation 0 hides one
generation of sibships (ADR 0001 §3), so that F is a lower bound, and the
gap is understated.

Layer 2 gives independent evidence. Case `genedrop_base` reports the
innovation variance of children of inbred recorded parents, simACE minus
gene-drop, as INFO rows `d1_innov_var_g*`; its scored `innov_var_g*` rows
use only children of non-inbred parents (rwaples/simACE#77).

simACE answered the open question in ADR 0033 (rwaples/simACE#76): A/2
stays the default (`pedigree.inbreeding.mode: "off"`), and the `approx`
and `exact` modes scale the residual by √(1 − F̄). Scoring exact mode
waits for the simace release that ships it. Tracked in status.md.

## D2. The C-transmission residual: realised or analytic

Under vertical C transmission (simACE ADR 0025), an offspring's C is
`b_m C_mother + b_f C_father + u`, with one u per maternal household. simACE
draws u with covariance `Σ − S`, where Σ is the configured C covariance and S
is the realised covariance of the transmitted part over that generation's
offspring. Each generation's C covariance is then Σ in expectation, whatever
the parents' realised C variance was. The textbook law instead uses the
analytic residual `Σ − B_m Σ B_m − B_f Σ B_f`. That is also Σ in expectation
across generations, but it carries part of the previous generation's
sampling error forward. ADR 0025 considers that rule and rejects it.

Case `ctm_base` scores simACE against its own (realised) law, and
`ctm_analytic` scores the same pedigrees against the analytic law. The
gap would show as a dispersion mismatch in `var_C*` and `var_L*` (per-seed
z spread away from 1), not as a mean shift. The 2026-10-21 run showed none
at R = 50 and M = 199. For every `var_C*_g*` and `var_L*_g*` row in
`ctm_analytic`:

- mean z −0.20 to +0.12
- SD of per-seed z 0.82 to 1.19 (ctm_base: 0.88 to 1.19)
- smallest p-value 0.0328 (`var_C1_g2` KS), against ctm_base's
  0.00922 for the same row

So these diagnostics cannot tell the two laws apart with three recorded
generations at N = 2000. ctm_base still checks simACE against the law it
documents.

## D3. Gaussian liability in simACE's phenotype and analysis stages

Some of simACE's stages assume liability is Gaussian, by design:

- the threshold models mark a case where standardised liability exceeds
  `ndtri(1 − K)` (simACE `phenotype/models/_prevalence.py:130` at v0.2.0),
  so realised prevalence is K only for normal liability;
- ADuLT and LTM onset map liability through a normal survival function
  (`phenotype/models/adult.py:224`);
- tetrachoric correlations, and the observed-scale h² built on affection
  status, read a bivariate-normal liability.

Layer 4 (ADR 0002) hands simACE a finite-locus A from a gene drop, whose
distribution is not exactly normal, so these quantities can depart from
their Gaussian targets without any defect in simACE. Case `simace_ts_p2`
reports them as INFO rows (`info/tetrachoric_*`, `info/affected_corr_*`,
`info/h2_observed_*`) and scores only quantities whose expectation needs
second moments alone. Its phenotype model is Weibull frailty, whose hazard
scales with exp(β L) and does not use `ndtri`; its prevalence
(`info/prevalence`) is INFO because no target is derived for it, not
because of this gap.

Measured on the 50 scored seeds (9000 to 9049, run of 2026-11-05) by
`pixi run python simace_ts/sweep_info.py`. Frailty onset and censoring move
these quantities far from the liability-scale targets whatever the
distribution of A (tetrachoric FS 0.077 against a liability correlation
near 0.3), so no Gaussian target is available. The measure is instead each
architecture's paired per-seed difference from c1m_a05, the architecture
with the most causal sites; c1m_a05 is itself finite-locus, so this cannot
show a departure that all architectures share.

Over the 15 computable statistics and 5 architectures (75 comparisons;
`h2_observed_falconer` is empty with no MZ twins), no single comparison
rejects: the largest |t| is 2.87 (`h2_observed_po`, c1k_a05, mean
−0.0028), below the two-sided Bonferroni cutoff of 3.63. The largest
absolute mean difference is 0.0045 (`tetrachoric_FO`, c10k_a05);
single seeds differ by far more. But c1m_a05 sits above the other
architecture in 65 of the 75 mean differences, 12 to 14 of 15 in each
architecture; only the MHS rows and `h2_observed_cousins` break the
pattern. The comparisons share c1m_a05 and the seeds' phenotypes, so they
are correlated, and the sign count is not a test. It does not separate an
unusually high c1m_a05 at these seeds from a real effect of the
architecture, and with c100_a05 no further from c1m_a05 than c10k_a05 or
c100k_a05 there is no sign of a departure that grows as the number of
causal sites falls (the series is c100, c1k, c10k, c100k and c1m at alpha
−0.5; c1k_a0 changes alpha instead). The cause is unresolved; a joint test
across statistics, with c1m_a05 against the other architectures, would
settle whether it is noise.

## D4. The xftsim matched configuration (layer 5)

Layer 5 (ADR 0003) compares simACE with xftsim at a configuration where
the two population laws agree to first order. They agree only approximately:

- **C as a household draw only.** xftsim has no household component.
  `am/run_xft.py` builds one: the mother's private noise draw, passed to
  her children by a vertical component (ADR 0003, Context). It is checked
  for one trait (`amc_r05`, C 0.2, r 0.5).
  - xftsim's built-in regimes make no half-sibs, so the maternal/paternal
    half-sib split of C is not compared.
  - Two-trait AM with C, and C transmission (rwaples/simACE-validate#5),
    are not checked.
- **Finite loci.** xftsim has M = 1000 unlinked loci, while simACE is
  infinitesimal. Case `am_finite_m` (INFO) compares xftsim at M = 250 with
  M = 4000, at r = 0.5 and 0.8. Its rows give the mean Var(A) gap and Welch
  t by generation. In the pilot, the M = 1000 gap from the Bulmer recursion
  was bounded at 0.005 to 0.009 in Var(A) by 95% CI (ADR 0003, decision 5).
- **Family structure.** simACE's `standard` mating (λ 0.5) gives variable
  sibship sizes and half-sibs. xftsim's configuration has one mate each and
  two children per couple. Only population moments are compared, never
  pedigrees. The spreads over seeds differ by design, so rows are scored by
  a Welch t on the mean, not a distribution test.
- **The pairing law.** Both simulators pair mates by a Gaussian copula of
  correlation r on liability ranks, which is the same law only for Gaussian
  liability and in the large-N limit. Neither side's realised correlation
  is exactly r:
  - xftsim's targets R = r(n−1)/n with one component;
  - simACE v0.3.0 pins its single-trait slot-level correlation to within
    5·10⁻⁴ of r by a Metropolis step after the copula (simACE ADR 0032);
    v0.2.0's copula alone fell slightly short (simACE #74). Layer 5's
    `r_mate_L` is over distinct couples with recorded children
    (`am/stats.py:74-75`), not slots, so it still scatters: one v0.3.0
    seed at r 0.8 (990001) gave 0.7946 to 0.8016 over g 2 to 6;
  - simACE's two-trait path pins all four correlations by Metropolis.
- **Two-trait AM** is checked only for rank-1 mate correlation matrices:
  every mate correlation equal to r (xftsim's linear regime), and
  r1 = 0.4, r2 = 0.1 with cross terms √(r1·r2) = 0.2. The unequal case
  pairs by `Rank1AssortativeMatingRegime` in `am/run_xft.py`, our code on
  xftsim's interface, so that case's pairing law is not independent of
  this repo (ADR 0003, decision 2). Other cross terms are unchecked,
  including simACE's derived one (`_cross_am_matrix`), which is 0 for
  independent traits.
- **Labels.** xftsim pairs a random split of the population, so its
  "mother" and "father" are labels. The mother–offspring and
  father–offspring rows are INFO.
