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
  - simACE v0.2.0's single-trait copula falls slightly short (simACE #74);
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
