# ADR 0001: Checks are conditional calibrations, scored in one Bonferroni family

**Status:** accepted; §2, §4 and §6 amended for layer 4 by ADR 0002; §1 and §2 amended for layer 5 by ADR 0003; §2 amended on 2026-11-08 (status.md G4): the chi-square's exact reference
**Date:** 2026-10-21

## Context

pedsum-validate and pg-phenotype-validate compare deterministic numbers with
a reference at a stated tolerance. simACE's output is random, so a reference
value and a tolerance do not exist. `simace validate` (simACE ADR 0028)
already tests simACE's residuals against its configured parameters, with
simACE's own code on both sides. This repo needs evidence from code that
simACE does not share, and a verdict whose false-positive rate is known.

## Decision

1. **Condition on the pedigree simACE produced.** An oracle draws M
   realisations of the quantity under test on that same pedigree, and
   computes every statistic on simACE's values and on each draw with one
   function. The truth is therefore the realised pedigree's, never a nominal
   kinship.
2. **Score ranks, not distances.** If simACE samples from the oracle's law,
   the randomised PIT of its value among the draws is exactly U(0, 1)
   (calibration.py). Over R seeds, a KS test and a chi-square test on equal
   bins score each statistic.

   *Amendment (2026-11-08, G4).* The chi-square's p-value comes from the
   exact null of Pearson's statistic under multinomial(R, 1/bins), not from
   the asymptotic chi2(bins − 1). The statistic is (bins/R)·Σc² − R, so its
   law is that of the sum of squared counts. A DP over the bins, in
   probabilities, computes that law exactly (`calibration.chi2_null_tail`).
   It holds while R·log10(bins) ≤ 300, which is R ≤ 300 at 10 bins, and it
   refuses a larger R. The p-value is P(S ≥ s_obs), with no randomisation.
   With 10 bins and α/K = 3.67e-05, the exact null pmf gives these sizes,
   P(p < α/K):
   - The asymptotic reference: 7.17e-05 at R = 50 and 5.70e-05 at R = 100,
     up to twice the nominal rate.
   - The exact test: 3.24e-05 at R = 50 and 3.46e-05 at R = 100. It is
     conservative only because the statistic is discrete.

   `tools/check_chi2_null.py` prints these sizes, and also the DP against
   enumeration and a simulation of both. Rejected alternatives:
   - A randomised exact p-value. Its size would be exactly α, but it adds a
     jitter stream per row, and F1 showed how much borderline rows move
     with jitter.
   - A seeded Monte Carlo reference. Its p-value at α/K has a relative
     standard error of 0.165 with 1e6 draws and 0.052 with 1e7, where the
     DP has none.
3. **Recorded generation 0 is given, not drawn.** simACE records generation
   0 after one round of reproduction (simACE `simulate.py`, the
   `reproduce()` call before `_fill_pedigree_slice`), so its rows hide
   sibships. Oracles condition on generation-0 values and draw generations
   1 and later.
4. **One Bonferroni family per run.** α = 0.05 / K, where K is the number
   of p-values over every scored row of every selected case: two for a
   calibration row (KS and chi-square), one for a t-test row. Layers are a
   grouping for reporting only. A run with `--only` or `--replicate` is
   partial: K covers its cases only, and it writes to `.work/results/` with
   the selection named in its results.md, so `results/` always holds a full
   run.
5. **Statistics fixed by construction are not rows.** Examples are MZ A
   correlation and the C correlation of full sibs. The maternal half-sib C
   correlation is fixed only while the father transmits no C (simACE ADR
   0025); with b_f > 0 it is a row.
6. **Power is measured, not assumed.** `power.py` runs the harness with
   stated oracle errors (mutations). Its `null` mutation measures the
   family-wise false-positive rate. R and M are chosen from that table.
7. **Oracles are called, never copied** (as in pg-phenotype-validate ADR
   0001). GPL tools are run as binaries (ped-sim) or
   loaded as libraries (nadiv), never vendored.

## Consequences

- A FAIL means simACE's statistic is not exchangeable with draws from the
  stated law on its own pedigree. It is not a measured bias. The mean z
  column shows the size and sign of the gap.
- A known gap between simACE's law and an outside model (for example, the
  infinitesimal segregation variance under inbreeding, D1) is an INFO case
  with an entry in differences.md, reported with its measured size. It is
  not scored: the D1 gap is too small to reject at α/K with these R and M,
  so a scored "expected rejection" would fail on its own. A gap that is
  detectable would need its affected statistics named, so that the rest of
  the case is still scored.
- Rows are correlated within a seed, because they are statistics of the
  same pedigree. Bonferroni stays valid under that dependence, though it is
  conservative.
