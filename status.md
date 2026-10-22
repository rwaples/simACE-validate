# Status of non-PASS rows

Run of 2026-10-21 against simace 0.2.0 (v0.2.0).

**Totals:** INFO 248, PASS 361. alpha/K = 7.03e-05.

Status vocabulary:
- **Open:** needs a decision or a fix.
- **Blocked:** waiting on another repo.
- **Accepted:** a known, documented gap.
- **Info:** reported, never scored.

| Id | Rows | Result status | Status | Next step | Updated |
|---|---|---|---|---|---|
| F1 | rm_base `r_E12_g3` | PASS since the C-transmission run (KS p 1.11e-4 at alpha/K 7.03e-5; mean z +0.35). Was FAIL (KS p 1.06e-4 at alpha/K 1.08e-4) | Accepted | A false positive at the pinned seeds, confirmed by an independent Codex review (below). It now passes only because K grew and the PIT jitter changed (below), not because anything was fixed | 2026-10-21 |
| D1 | deep_infinitesimal, every row | INFO (`var_A*_g*` mean z +0.03 to +0.52, largest at g ≥ 7) | Info | Ask simACE whether A/2 should scale with parental F (differences.md D1) | 2026-10-21 |
| G1 | layer 3b, simACE `standard` mating | not checked | Open | Write `slim/standard.slim` from simACE's documented mating model (a second implementation, not an outside model) | 2026-10-21 |
| D2 | ctm_analytic, every row | INFO (`var_C*_g*`, `var_L*_g*`: mean z −0.20 to +0.12, SD of z 0.82 to 1.19) | Info | none: no gap detectable between the realised and analytic C residual at R = 50 (differences.md D2) | 2026-10-21 |
| G2a | C transmission (simACE ADR 0025) | checked by ctm_base: 80 rows, all PASS | Closed | none | 2026-10-21 |
| P1 | results/power.md, every row | Generated with the old power seeds `seed0 + 1_000_000 + 1000·rep`, so rm_base reps 5 to 10 reused ctm_base reps 1 to 6's simACE seeds; the two cases' rows are not independent | Open | Rerun `pixi run power`, which now seeds each case from its id (about 3 h 10 min) | 2026-10-22 |
| G2b | assortative mating | not checked | Open | The layer-1 oracle needs that law. Under AM the liability covariance is not 2φA | 2026-10-21 |
| G3 | genedrop_base `innov_var_g*` | PASS, but the null (simACE − gene-drop = 0) is exact only for non-inbred parents; expected shift 0.00025·A at g3, about 0.05 SE | Open | Design the comparison under inbreeding (rwaples/simACE#77) | 2026-10-21 |
| I1 | rows with fewer than 20 seeds (MZ rows of the deep cases), and genedrop `spread_*` | INFO | Info | none | 2026-10-21 |

## F1: a false positive at the pinned seeds (accepted)

The evidence that simACE is not at fault:

- **Fresh seeds pass.** `pixi run harness --only rm_base --replicate 1` (seeds 101000 to 101049; written to `.work/results/`) passes every row. Its smallest KS p is 0.046, its smallest p over both tests 0.016 (`r_C12_g2` chi-square), and its `r_E12_g3` KS p 0.72.
- **A direct check finds no gap.** `simace simulate` over 300 fresh rm_base seeds (500000 to 500299) gives this cross-trait r(E1, E2) by generation, with z against rE = 0.1:

  | Generation | r(E1, E2) | z |
  |---|---|---|
  | 0 | 0.1002 ± 0.0014 | +0.18 |
  | 1 | 0.1014 ± 0.0013 | +1.07 |
  | 2 | 0.1002 ± 0.0013 | +0.18 |
  | 3 | 0.1004 ± 0.0014 | +0.30 |

- **The pinned seeds are borderline.** Across the 50 pinned seeds (1000 to 1049), generation 3 gives 0.1075 ± 0.0029.
- **One false alarm is expected now and then.** With `rm_base` scored as its own Bonferroni family, the measured false-positive rate was 1/10 (`results/power.md`, null row, before the per-statistic PIT streams; 0/10 since, exact 95% interval 0 to 0.308). The full run's threshold is stricter, so that rate is an upper bound for rm_base there. The row sits just under its bound.
- **An independent review agrees.** Codex (2026-10-21) read the E draw in simACE v0.2.0 and in `oracles/liability.R`, the statistic, and the randomised PIT and exact KS test. It found no defect that explains F1, and gave the verdict "likely false positive". Its own checks:
  - another 300 fresh seeds (900000 to 900299): generation-3 r(E1, E2) = 0.0995 ± 0.0011
  - an exact replay of F1: KS p 1.064e-4 against alpha/K 1.080e-4 (K = 463)
  - re-jittering the same fixed ranks: KS p 0.92e-4 to 1.26e-4, either side of the cutoff
  - float32 storage of E: changes the correlation by at most 2.3e-9

  Its caveats: the null rate of 1/10 has a 95% interval of about 0.25% to 44.5%, and the chi-square test is approximate with sparse bins.

Decision: keep the pinned seeds and carry F1 as an accepted false positive. Changing seed0 after a failure would be selection after the fact. Until the C-transmission run (below), the harness exited 1 on this row.

**Update, C-transmission run (2026-10-21).** F1 passes now, at the same seeds. Two changes cause this. Both were fixed in the plan before any result was seen:

- **K grew from 463 to 711.** The causes are case ctm_base, the split of PO into MO and FO, and the new `var_C` rows. alpha/K fell from 1.08e-4 to 7.03e-5. This alone makes F1 pass. A replay with the old jitter gives KS p 1.064e-4, which is above 7.03e-5.
- **The PIT jitter now has one stream per statistic** (`run_harness.py`), so adding a statistic no longer re-jitters the others. This moved F1's KS p from 1.064e-4 to 1.11e-4 (the replay was Codex's; 1.11e-4 is the run's).

simACE and F1's ranks are unchanged. F1 stays recorded as an accepted false positive at these seeds.
