# Status of non-PASS rows

Run of 2026-11-05 against simace 0.2.0 (v0.2.0).

**Totals:** INFO 658, PASS 819. alpha/K = 4.06e-05 (K = 1231).

Status vocabulary:
- **Open:** needs a decision or a fix.
- **Blocked:** waiting on another repo.
- **Accepted:** a known, documented gap.
- **Info:** reported, never scored.
- **Closed:** the gap is now checked; kept for the record.

| Id | Rows | Result status | Status | Next step | Updated |
|---|---|---|---|---|---|
| F1 | rm_base `r_E12_g3` | PASS in the run of 2026-11-05 (KS p 1.11e-04 at alpha/K 4.06e-05; mean z +0.35). Was FAIL (KS p 1.06e-4 at alpha/K 1.08e-4) | Accepted | A false positive at the pinned seeds, confirmed by an independent Codex review (below). It now passes only because K grew and the PIT jitter changed (below), not because anything was fixed | 2026-11-05 |
| D1 | deep_infinitesimal, every row; genedrop_base `d1_innov_var_g*` | INFO (`var_A*_g*` mean z +0.03 to +0.52, largest at g ≥ 7) | Blocked | simACE kept A/2 as the default and added `approx`/`exact` modes (ADR 0033, rwaples/simACE#76). Add an exact-mode case (no burn-in, so recorded F is simACE's F) that scores every child in layer 2 and Var(A1) in layer 1 (ADR 0033 §10), once a simace release ships the modes | 2026-11-04 |
| G1 | layer 3b, simACE `standard` mating: std_mating (7 rows), std_structure (55 rows) | PASS (smallest p: std_mating 0.178, std_structure 0.017) | Closed | Checked two ways (README, layer 3b). std_structure is a second implementation of simACE's documented rules, so it cannot catch a flaw in the model itself. Doc gaps settled from simACE's code are reported in rwaples/simACE#78. The power of both cases is in results/power_std_mating.md and power_std_structure.md | 2026-10-22 |
| D2 | ctm_analytic, every row | INFO (`var_C*_g*`, `var_L*_g*`: mean z −0.20 to +0.12, SD of z 0.82 to 1.19) | Info | none: no gap detectable between the realised and analytic C residual at R = 50 (differences.md D2) | 2026-10-21 |
| G2a | C transmission (simACE ADR 0025) | checked by ctm_base: 80 rows, all PASS | Closed | none | 2026-10-21 |
| P1 | results/power.md, every row | Was generated with the old power seeds `seed0 + 1_000_000 + 1000·rep`, so rm_base reps 5 to 10 reused ctm_base reps 1 to 6's simACE seeds. Rerun with per-case seeds at 28847ef: null 0/10 on both cases, every mutation 10/10 | Closed | none. The rerun was 10 GenomeDK jobs, one per case and mutation (6 CPUs each, 9 to 15 min), merged with `--merge` | 2026-10-22 |
| G2b | assortative mating: layer 5 (ADR 0003), am_r0, am_r05, am_r08, am2_r02, am2_r04, amc_r05 (396 rows); am_finite_m (INFO) | PASS in the run of 2026-11-05 (smallest p 6.4e-4, am2_r02 `r_w_A_g6`; am_r08 `r_mate_L` p ≥ 0.247 at every g, so no D entry for simACE #74). Power (results/power_am.md, 10 reps): every mutation 10/10; null 0/10 (am_r05) and 1/10 (am2_r04, one row) | Closed | none for the scope tested: single-trait AM (C 0, and household C 0.2 at r 0.5) and two-trait AM with every mate correlation equal, independent traits (differences.md D4). Unequal r (rank 1) is G5 | 2026-11-05 |
| G5 | layer 5, every case; am2_r0401 (two traits, r1 0.4, r2 0.1, rank 1; 114 rows) | ERROR until the references are regenerated: am2_r0401 moved its pairing regime into am/run_xft.py and its Bulmer matrix into am/stats.py, so every committed reference's manifest is stale, and am2_r0401 has none. Pilot (40 seeds per side, am/pilot/rank1_summary.md): smallest Welch p 0.0147 over 139 rows | Open | Run `pixi run -e xftsim xft-reference` (every reference, about 3 h at 6 jobs; power mutations with `--case am_r05,am2_r04 --mutation r_090,r_0,mate_A,mate_y1`), then the full harness (K 1345, alpha/K 3.72e-05) | 2026-11-05 |
| G3 | genedrop_base `innov_var_g*` | Partial run of 2026-11-04: PASS (g3 mean −0.0078, t −1.39); `d1_innov_var_g3` INFO (21 seeds, mean +0.106, t +1.62) | Closed | none. The rows now score only children whose recorded parents both have F = 0, where the null is exact (rwaples/simACE#77, `genedrop/layer2.py`) | 2026-11-04 |
| G4 | every calibrated row (layers 1 and 3b): the chi-square test | PASS, but its Pearson reference is asymptotic. Codex simulated the exact multinomial null at R = 100, 10 bins, alpha/K = 0.05/587: rejection 0.000115 (95% CI 0.000108 to 0.000122) against a nominal 0.0000852 | Open | Give the chi-square an exact multinomial or Monte Carlo reference. That rescores every case, so it is a separate change | 2026-10-22 |
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
