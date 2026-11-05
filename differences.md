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

Measured on the 50 scored seeds (9000 to 9049, run of 2026-11-05). Frailty
onset and censoring move these quantities far from the liability-scale
targets whatever the distribution of A (tetrachoric FS 0.077 against a
liability correlation near 0.3), so the measure is each architecture's
paired per-seed difference from c1m_a05, the closest to Gaussian. Over the
15 statistics and 5 architectures (75 comparisons), the largest |t| is 2.87
(`h2_observed_po`, c1k_a05, −0.0028), below the 3.63 that Bonferroni over 75
needs at 0.05, and the paired differences are at most 0.0045 in size
(`tetrachoric_FO`, c10k_a05). They do not grow as the number of causal sites
falls: c10k_a05 and c100k_a05 differ from c1m_a05 as much as c100_a05 does.
The FO and PO rows sit below c1m_a05 in every architecture, which points
to c1m_a05 being high at these seeds rather than to the others being low.
No departure is detected at this sample size; a smaller one remains possible.
