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

Open question for simACE: whether A/2 is the intended model or should
scale with parental F. Tracked in status.md.
