# Status of non-PASS rows

Run of 2026-10-21 against simace 0.2.0 (v0.2.0).

**Totals:** FAIL 1, INFO 138, PASS 236. alpha/K = 0.000108.

Status vocabulary:
- **Open:** needs a decision or a fix.
- **Blocked:** waiting on another repo.
- **Accepted:** a known, documented gap.
- **Info:** reported, never scored.

| Id | Rows | Result status | Status | Next step | Updated |
|---|---|---|---|---|---|
| F1 | rm_base `r_E12_g3` | FAIL (KS p 1.06e-4 at alpha/K 1.08e-4; mean z +0.35) | Open | Decide below | 2026-10-21 |
| D1 | deep_infinitesimal, every row | INFO (`var_A*_g*` mean z +0.03 to +0.52, largest at g ≥ 7) | Info | Ask simACE whether A/2 should scale with parental F (differences.md D1) | 2026-10-21 |
| G1 | layer 3b, simACE `standard` mating | not checked | Open | Write `slim/standard.slim` from simACE's documented mating model (a second implementation, not an outside model) | 2026-10-21 |
| G2 | C transmission (simACE ADR 0025), assortative mating | not checked | Open | The layer-1 oracle needs those laws. Under AM the liability covariance is not 2φA | 2026-10-21 |
| I1 | rows with fewer than 20 seeds (MZ rows of the deep cases), and genedrop `spread_*` | INFO | Info | none | 2026-10-21 |

## F1: likely a false positive at the pinned seeds

The evidence that simACE is not at fault:

- **Fresh seeds pass.** `pixi run harness --only rm_base --replicate 1` (seeds 101000 to 101049; written to `.work/results/`) passes every row, and its smallest p is 0.046.
- **A direct check finds no gap.** `simace simulate` over 300 fresh rm_base seeds (500000 to 500299) gives this cross-trait r(E1, E2) by generation, with z against rE = 0.1:

  | Generation | r(E1, E2) | z |
  |---|---|---|
  | 0 | 0.1002 ± 0.0014 | +0.18 |
  | 1 | 0.1014 ± 0.0013 | +1.07 |
  | 2 | 0.1002 ± 0.0013 | +0.18 |
  | 3 | 0.1004 ± 0.0014 | +0.30 |

- **The pinned seeds are borderline.** Across the 50 pinned seeds (1000 to 1049), generation 3 gives 0.1075 ± 0.0029.
- **One false alarm is expected now and then.** With `rm_base` scored as its own Bonferroni family, the measured false-positive rate is 1/10 (`results/power.md`, null row). The full run's threshold is stricter, so that rate is an upper bound for rm_base there. The row sits just under its bound.

Changing seed0 because a row failed is selection after the fact, so it has not been done. The decision for the maintainer:

1. keep the pinned seeds and carry F1 as a known false positive, with the harness exiting 1, or
2. re-pin every case's seeds once, by a rule fixed in advance (for example, seed0 from `case_seed(id)`), and accept that run's outcome.
