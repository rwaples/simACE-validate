# ADR 0001: Checks are conditional calibrations, scored in one Bonferroni family

**Status:** accepted
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
