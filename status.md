# Status of non-PASS rows

Run of 2026-11-06 against simace 0.3.0 (v0.3.0).

**Totals:** ERROR 1, INFO 748, PASS 937. alpha/K = 3.7e-05 (K = 1353). The ERROR is simace_ts_p2 (G5): its committed Parts were built under simace 0.2.0.

Status vocabulary:
- **Open:** needs a decision or a fix.
- **Blocked:** waiting on another repo.
- **Accepted:** a known, documented gap.
- **Info:** reported, never scored.
- **Closed:** the gap is now checked; kept for the record.

| Id | Rows | Result status | Status | Next step | Updated |
|---|---|---|---|---|---|
| F1 | rm_base `r_E12_g3` | PASS in the run of 2026-11-06 against simace v0.3.0 (KS p 0.388, chi2 0.911; mean z −0.11): v0.3.0's standard mating (per-birth MZ twins and family order, simACE #78) changes the pedigrees at the pinned seeds. PASS in the run of 2026-11-05 (KS p 1.11e-04 at alpha/K 4.06e-05; mean z +0.35). Was FAIL (KS p 1.06e-4 at alpha/K 1.08e-4) | Accepted | A false positive at the pinned seeds, confirmed by an independent Codex review (below). It now passes only because K grew and the PIT jitter changed (below), not because anything was fixed | 2026-11-05 |
| D1 | deep_infinitesimal, every row; genedrop_base `d1_innov_var_g*` | INFO (`var_A*_g*` mean z +0.03 to +0.52, largest at g ≥ 7) | Blocked | simACE kept A/2 as the default and added `approx`/`exact` modes (ADR 0033, rwaples/simACE#76). Add an exact-mode case (no burn-in, so recorded F is simACE's F) that scores every child in layer 2 and Var(A1) in layer 1 (ADR 0033 §10), once a simace release ships the modes | 2026-11-04 |
| G1 | layer 3b, simACE `standard` mating: std_mating (7 rows), std_structure (55 rows) | PASS (smallest p: std_mating 0.178, std_structure 0.017). Under simace 0.3.0, before the nulls followed its per-birth MZ model, both cases were ERROR ("a couple has two MZ pairs"); with the v0.3.0 nulls, `--only std_mating,std_structure` on 2026-11-06 is PASS (smallest p: std_mating 0.0279, std_structure 0.00872, at alpha/K 0.000403) | Closed | Checked two ways (README, layer 3b). std_structure is a second implementation of simACE's documented rules, so it cannot catch a flaw in the model itself. Doc gaps settled from simACE's code are reported in rwaples/simACE#78; v0.3.0 also fixed the MZ model and id order there (rwaples/simACE-validate#1). The power of both cases at v0.3.0 is in results/power_std_mating.md and power_std_structure.md: every mutation rejects 40/40 (20/20); null 3/40, null_slim 2/40, std_structure null 1/20 | 2026-11-06 |
| D2 | ctm_analytic, every row | INFO (`var_C*_g*`, `var_L*_g*`: mean z −0.20 to +0.12, SD of z 0.82 to 1.19) | Info | none: no gap detectable between the realised and analytic C residual at R = 50 (differences.md D2) | 2026-10-21 |
| G2a | C transmission (simACE ADR 0025) | checked by ctm_base: 80 rows, all PASS | Closed | none | 2026-10-21 |
| P1 | results/power.md, every row | Was generated with the old power seeds `seed0 + 1_000_000 + 1000·rep`, so rm_base reps 5 to 10 reused ctm_base reps 1 to 6's simACE seeds. Rerun with per-case seeds at 28847ef: null 0/10 on both cases, every mutation 10/10 | Closed | none. The rerun was 10 GenomeDK jobs, one per case and mutation (6 CPUs each, 9 to 15 min), merged with `--merge` | 2026-10-22 |
| G2b | assortative mating: layer 5 (ADR 0003), am_r0, am_r05, am_r08, am2_r02, am2_r04, am2_r0401, amc_r05 (510 rows); am_finite_m (INFO) | PASS in the run of 2026-11-06 against simace v0.3.0 (smallest p 0.00177, am_r08 `r_FS_L_g6`; am2_r0401 smallest p 0.00178, `r_mate_L_t2_g2`, against its clean-room reference (G7); am_r08 `r_mate_L` p ≥ 0.0503 at every g). Was PASS against v0.2.0 in the run of 2026-11-05. Power at v0.3.0 (results/power_am.md, 10 reps): every mutation 10/10; null 0/10 (am_r05) and 1/10 (am2_r04, one row, `r_FS_A_t2_g3`) | Closed | none for the scope tested: single-trait AM (C 0, and household C 0.2 at r 0.5) and two-trait AM with a rank-1 mate correlation matrix (every correlation equal, or r1 0.4, r2 0.1), independent traits (differences.md D4) | 2026-11-06 |
| G7 | layer 5, every case; am2_r0401 (two traits, r1 0.4, r2 0.1, rank 1; 114 rows) | Every reference regenerated on 2026-11-06 after am2_r0401 moved its pairing regime into am/run_xft.py and its Bulmer matrix into am/stats.py; the full run against simace v0.3.0 passes every layer-5 row (G2b). Regenerated again the same day: the rank-1 regime was rewritten clean-room (its first version followed xftsim's GPL LAM code, ADR 0003 decision 2), the manifest now hashes am/reference.py and no longer am/bulmer.py, and ruff left the xftsim env. Only am2_r0401's run data changed (the other 17 references are identical bar their manifests); the full rerun passes every layer-5 row | Closed | none | 2026-11-06 |
| G3 | genedrop_base `innov_var_g*` | Partial run of 2026-11-04: PASS (g3 mean −0.0078, t −1.39); `d1_innov_var_g3` INFO (21 seeds, mean +0.106, t +1.62) | Closed | none. The rows now score only children whose recorded parents both have F = 0, where the null is exact (rwaples/simACE#77, `genedrop/layer2.py`) | 2026-11-04 |
| G4 | every calibrated row (layers 1 and 3b): the chi-square test | PASS, but its Pearson reference is asymptotic. Codex simulated the exact multinomial null at R = 100, 10 bins, alpha/K = 0.05/587: rejection 0.000115 (95% CI 0.000108 to 0.000122) against a nominal 0.0000852 | Open | Give the chi-square an exact multinomial or Monte Carlo reference. That rescores every case, so it is a separate change | 2026-10-22 |
| I1 | genedrop `spread_*` | INFO | Info | none. Under simace v0.3.0's per-birth MZ model every deep-case seed has MZ pairs (n 50, was 4), so deep_constant's four MZ rows are now scored and PASS; deep_infinitesimal's stay INFO under D1 | 2026-11-06 |
| G5 | layer 4, simace_ts_p2 (simACE-ts, rwaples/simACE-validate#2) | Partial run of 2026-11-07 under simace 0.3.0 (`--only simace_ts_p2`, 50 seeds 9000 to 9049; output in `.work/results/`, not committed): the 11 scored rows PASS (largest absolute t 2.06, `lc_FS`, p 0.0447; alpha/K 0.00455 at K = 11) and the 22 exact rows PASS with 0 mismatches (largest difference 7.2e-14); the sweep's 22 exact rows have 0 mismatches over 5 architectures (largest 7.4e-14). At a full run's alpha/K of 5.91e-5 (K = 846, ADR 0002 item 10) the minimum detectable gap is 0.0016 (`dz_L`) to 0.0067 (`po_g9_slope`) (`simace_ts/sweep_info.py --mde-k 846`). The 0.2.0 Parts (1f8bd5d) also passed | Open | Score it in a full run, where K covers every case. This layer's oracle consumes simACE's pedigree and feeds simACE's own stages (README, layer 4) | 2026-11-07 |
| D3 | simace_ts_p2 `info/tetrachoric_*`, `info/affected_corr_*`, `info/h2_observed_*` | INFO. Against c1m_a05 no single comparison rejects (simace 0.3.0: largest absolute paired t 2.76 of 75, Bonferroni cutoff 3.63). c1m_a05 sits above the other architecture in 51 of the 75 mean differences, down from 65 under 0.2.0 (differences.md D3) | Open | Run a joint test across statistics of c1m_a05 against the other architectures, to settle whether the pattern is noise | 2026-11-07 |
| G6 | layer-2 cross-check on the layer-4 pedigree (check 4 of #2) | not checked | Open | Run genedrop/drop.py on the same N = 35,951 pedigree and compare its IBD and A correlations with the tskit drop's, as INFO. Its cost at this size is unmeasured | 2026-11-05 |

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
