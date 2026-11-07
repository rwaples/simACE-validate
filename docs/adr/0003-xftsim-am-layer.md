# ADR 0003: Assortative mating is scored as a two-sample comparison against xftsim

**Status:** accepted 2026-11-05 (rwaples/simACE-validate#4).
- The user chose option (b) for the simace pin (decision 8) and added
  two-trait AM at two levels (decision 2).
- Decision 8 is built: `am/manifest.py`, `am/reference.py`,
  `run_harness.py` (`run_xftsim`), `power.py` and cases.yaml.
- The user added a scored two-trait case with unequal mate correlations
  (`am2_r0401`, rank 1; decision 2) and accepted regenerating every
  reference for it.

**Date:** 2026-11-04

Amends ADR 0001 §1 (conditioning) and §2 (the test) for layer 5 only.
Answers the open questions of rwaples/simACE-validate#4.

## Context

G2b (status.md) is open because layer 1 has no law for assortative mating
(AM). ADR 0001 conditions on the pedigree simACE produced. Under AM that
pedigree was chosen on the liabilities. So, given the pedigree, simACE's A
is not drawn by midparent plus Mendelian sampling, and the layer-1 oracle
would have to condition on the pairing too. No tractable law does that.

[xftsim](https://github.com/border-lab/xftsim) (GPLv3; v0.3.1, commit
`8a4cdf0` of 2026-07-22) is a forward-time simulator that shares no code with simACE. It
has finite loci and its own AM regime, so it can serve as an empirical
reference at the population level.

Citation: Richard Border, *xftsim: eXtensible Forward Time SIMulator*,
https://github.com/border-lab/xftsim (formerly `rborder/xftsim`, which
redirects; `pixi.toml` installs from that URL), commit `8a4cdf0`, GPL-3.0.
Author, licence and commit date are from the repository and its `setup.py`
at that commit (checked 2026-11-06). xftsim has no paper or citation file:
its README and docs ask for none, and Border et al. (2022), *Science*
378:754–761, doi:10.1126/science.abo2059, the closest paper, does not name
it (its simulations use the authors' own R code).

The AM laws at the pinned versions:

- **simACE, single trait** (`assort1`; v0.2.0
  `simace/simulation/simulate.py:699-726`, v0.3.0 `:809-858`).
  - It ranks mating slots by liability A + C + E in each sex, draws M
    bivariate normals with correlation r, and pairs the female and male
    slots whose ranks match each draw's ranks: a Gaussian copula of
    correlation r on the liability ranks.
  - It then swaps out duplicate couples (`:728-749`). A parent with k
    partners fills k slots.
  - In v0.2.0 the realised correlation falls slightly below r (simACE #74;
    `plans/issue-69-mate-correlation-test.md` reports a bias of −0.0017
    to −0.0029 at N 2k, λ 0.5). v0.3.0 pins it with a 1-D Metropolis
    step after the copula, to within `METROPOLIS_TOL` = 5·10⁻⁴ (simACE ADR
    0032; v0.3.0 `simulate.py:886`, `_refine_single_trait`).
- **simACE, both traits** (`assort1` and `assort2` both non-zero,
  optional `assort_matrix`; v0.2.0 `simulate.py:601-696`, v0.3.0
  `:725-808`).
  - Each female's ideal mate is drawn from the Gaussian conditional, given
    the target 2×2 mate-correlation matrix R_mf and the realised
    within-sex trait correlations. Males are matched to those draws, and a
    Metropolis step pins all four mate correlations to within 5·10⁻⁴.
  - Both cross terms take `R_mf[0, 1]`.
  - The run stops if the realised correlations make R_mf infeasible
    (ADR 0024 in simACE).
- **xftsim `LinearAssortativeMatingRegime`** (`xftsim/mate.py:843-887` at
  `8a4cdf0`).
  - It splits the population at random into two halves, standardises each
    listed phenotype, and builds the score √R·Σy + √(1−R)·ε with
    independent noise ε in each half. Each half is sorted by its score and
    the two halves are paired rank to rank.
  - With one component, R = r, so for Gaussian y and large n the mates'
    phenotypes are BVN(r): the same Gaussian copula as simACE's single-trait
    path.
  - With K components it pairs on the standardised sum and sets every
    cross-mate correlation (y_i of one mate with y_j of the other, all i, j)
    to r. That is feasible only for r ≤ 1/K (its own warning) when the
    traits are independent.

The two laws therefore coincide in the large-n limit only under Gaussian
liability, and only when the mating variable is the same: the full
liability, noise included.
- **Single trait:** simACE `assort1 = r` matches xftsim with one component.
- **Two traits:** a 4-variate Gaussian is fixed by its covariance, so simACE
  with `assort1 = assort2 = r` and `assort_matrix [[r, r], [r, r]]`
  matches xftsim on both phenotypes at r.
- **C:** xftsim has no household component, but one can be built
  (construction by the xftsim-bench session, 2026-11-05). Each person gets a
  private noise trait hh (variance 1), never summed into y. A
  `LinearVerticalComponent` reads the mother's hh (input vorigin 0,
  `xftsim/index.py:873`), normalised and scaled by √C, into the child's y.
  Children of one mother then share one draw, independent of the parents'
  own y. That is simACE's household C with household = mother. `am/run_xft.py`
  builds it for one trait when C > 0, and leaves the C = 0 path unchanged.

The match is approximate, not exact, in several ways:

- At fixed finite M, xftsim's A is not exactly Gaussian. Pairing on noisy
  scores and simACE's rank copula need not give identical laws.
- xftsim standardises with the population SD but takes its target from
  `np.cov` (sample variance). With one component this gives
  R = r(n−1)/n: 0.79984 at r = 0.8 and n = 5000 (pointed out by Codex).
- Finite rank matching lowers both simulators' realised correlations
  slightly.
- simACE's slots (one per partnership), its dedup swaps, and the
  different family structures can move the means at finite N.

The pilot measures the combined effect. It supports approximate agreement
within the measured MDDs, not equality.

## Decision

1. **Subject: population moments, two-sample.** Layer 5 compares per-seed,
   per-generation statistics of simACE's pedigrees with those of
   independent xftsim runs at a matched configuration. Neither side is
   conditioned on the other's pedigree. The infinitesimal Bulmer
   recursions are reported beside each row as INFO:
   - `am/bulmer.py::bulmer` for one trait;
   - `bulmer2` for two traits. It is the multivariate form: mates' A
     cross-covariance is Σ_A Σ_y⁻¹ R_y Σ_y⁻¹ Σ_A.

   Alternatives considered:
   - *Score against the Bulmer recursion alone* (one-sample). It is a
     first-order law: it ignores finite-N drift, the copula's departure from
     exact normality at finite M, and dedup. It is theory that simACE's own
     docs share, not independent code. At harness R its approximation error
     could become detectable. Kept as an INFO column.
   - *Extend the layer-1 oracle to AM.* That needs the law of A given a
     pedigree that was itself selected on liability. Rejected as
     intractable.
2. **Matched configurations.** Two traits, each with A = 0.5, C = 0 and
   E = 0.5, except where a case sets trait 1's household C (then
   E = 0.5 − C). The traits are independent: `rA = rC = rE = 0`, passed
   explicitly, because simACE's CLI defaults are `rA` 0.5 and `rC` 0.3
   (`simulate.py:1669-1671`). Mating is on the full liability. There are no
   MZ twins, N = 10⁴ per generation, and there are 6 rounds of mating from
   founders.
   - **Single trait, r ∈ {0, 0.5, 0.8}.** simACE `assort1 = r`; trait 2 is
     a passenger. xftsim mates on y alone. r = 0 is check 3 of the issue,
     the random-mating baseline, run through the same xftsim code path.
   - **Household C, single trait, C = 0.2, r = 0.5** (`amc_r05`). simACE
     `C1 = 0.2`, `E1 = 0.3`; xftsim builds C as in Context. Under AM, both
     sides mate on the liability with C included, and the full-sib
     liability covariance gains C.
   - **Two traits, r ∈ {0.2, 0.4}, every mate correlation r.** simACE uses
     `assort1 = assort2 = r` and `assort_matrix [[r, r], [r, r]]`. xftsim
     mates on y1 and y2. With independent traits, R_mf = r·11ᵀ is feasible
     only for r < (1 + ρ_w)/2 = 0.5, where ρ_w is the within-person trait
     correlation (it starts at 0 and grows under this AM). 0.4 leaves
     margin; 0.5 or more would put generation 0 at or past the bound.
   - **Two traits, unequal r, rank 1: r1 = 0.4, r2 = 0.1** (`am2_r0401`).
     The mate correlation matrix is R_mf = v·vᵀ with v = (√r1, √r2), so the
     cross terms are √(r1·r2) = 0.2. simACE uses `assort1 = 0.4`,
     `assort2 = 0.1` and `assort_matrix [[0.4, 0.2], [0.2, 0.1]]`.
     - xftsim's LAM gives every mate correlation one r (it pairs on the
       equally weighted sum). Its `GeneralAssortativeMatingRegime` takes
       any target but needs LocalSolver, a commercial solver
       (`xftsim/mate.py:1026-1053` at `8a4cdf0`).
     - So `am/run_xft.py` adds `Rank1AssortativeMatingRegime`, written
       against xftsim's `MatingRegime` interface (`get_potential_mates`,
       `enumerate_assignment`). It was written clean-room on 2026-11-06,
       from a spec of this law and of that interface alone, by an agent
       that never read an xftsim regime's source; the first version had
       followed LAM's `mate` step by step, which put GPL-derived code in
       this MIT repo. Each side pairs on the
       score s = Z·S⁻¹v (Z the standardised traits, S their covariance on
       that side) at correlation ρ = √(q_f·q_m), with q = vᵀS⁻¹v. That
       pins all four mate correlations at v·vᵀ in every generation, as LAM
       re-solves its R for the realised ρ_w.
     - The first version reproduced LAM's output exactly at r1 = r2 (103
       statistics, one seed); the clean-room one draws its noise
       differently, so the check is now statistical. Its 40 fresh seeds
       (29000 to 29039) against the 100 seeds of the first version's
       `am2_r0401` reference give 162 Welch rows with smallest p 0.0015
       (`var_A_t2_g3`; Bonferroni cutoff 3.09·10⁻⁴), and 12 rows below
       0.05 where 8.1 are expected. The low rows are trait 2's Var(A) at
       g 2 to 4, where the old reference sat below the Bulmer recursion
       (z −4.15 and −3.74 at g 2 and 3) and the new code does not (|z| ≤
       2.26 at every g). Meiosis, effects, the architecture and the
       pedigree stay xftsim's; the pairing law of this case is ours.
     - Pilot, run with the first version: 40 seeds per side (`am/pilot/rank1_summary.md`). simACE
       against xftsim gives 139 Welch rows, smallest p 0.0147. Against the
       Bulmer recursion with this R (`bulmer2`), the largest |z| is 2.45
       for simACE and 1.94 for xftsim.
     - Feasible while ρ < 1; the regime raises otherwise rather than clamp
       as LAM does.
   - **simACE:** `G_ped = G_sim = 6` and `mating_lambda` 0.5. Recorded
     generation k follows k + 1 matings (ADR 0001 §3), so it is reported as
     g = k + 1.
   - **xftsim:** M = 1000 unlinked loci (`constant_map_from_haplotypes`,
     p = 0.5), one mate each, two offspring per couple, and fixed-scale
     effects (decision 3). In the two-trait configuration the traits'
     effects and noise are independent.
   - **Compared statistics.** One function runs on both sides
     (`am/stats.py::generation_stats` per trait; `cross_trait_stats`).
     - Per trait: Var(A) and Var(L) for g = 1 to 6.
     - For g = 2 to 6: the mate correlation on L and on A (over distinct
       couples), the full-sib correlation of A and L, and the
       parent–offspring correlation of A and L (mother and father pooled).
     - Two traits also: the within-person r(A1, A2) and r(L1, L2), the
       cross-trait mate correlations on L (both orders), and on A.
   - **Not compared:**
     - MO and FO separately. xftsim's `sex_aware=False` makes "maternal" a
       random label, so those rows are INFO.
     - Half-sibs. xftsim's built-in regimes cannot make them:
       `MatingRegime.enumerate_assignment` expands the female and male lists
       with the same per-female counts in the same order, so a female's
       second mating repeats the same male, and the duplicate couple is
       merged (`exhaustive=False` raises `NotImplementedError`;
       `xftsim/mate.py` and `utils.exhaustive_enumerate` at `8a4cdf0`).
       So the maternal/paternal half-sib split of household C cannot be
       checked against xftsim without writing our own regime.
     - Generation 0 of simACE, which is not recorded.

   Alternatives considered:
   - *C = 0 throughout* (the first draft). Lifted once the household-C
     construction was found; the pilot matches with C = 0.2 (below).
   - *C transmission (b_m, b_f > 0).* Moved to its own issue,
     rwaples/simACE-validate#5, which chains the parents' c through a
     vertical component. Not built here (Open items).
   - *Match family structure* by setting xftsim's mates per female to a
     ZTP. The population moments do not depend on it to first order.
     Rejected: only moments are compared, never pedigrees.
   - *Two-trait AM through simACE's derived cross term* (no
     `assort_matrix`; `_cross_am_matrix` derives it from ρ_w). xftsim has
     no regime with that structure. A cross-term-0 law (simACE's default
     with unequal r and `rA = 0`) is rank 2. No single-score pairing
     reaches it, so it would need xftsim's `GeneralAssortativeMatingRegime`
     and LocalSolver. Not piloted; unequal r is covered at rank 1 only.
3. **xftsim effects are fixed-scale.** β ~ N(0, 1) per trait is rescaled so
   that mean(β²) = h² exactly, through `AdditiveEffects(scaled=True,
   standardized=True)`. The founders' expected Var(A) is then h² given the
   allele frequencies. `GCTA_Architecture` leaves the scale random, with a
   relative SD of √(2/M) = 4.5% at M = 1000, which would dominate every
   Var(A) row's spread. Alternative: keep GCTA and score Var(A)/Var(A₀).
   Rejected, because simACE's founders are not recorded, so the ratio has
   no simACE counterpart.
4. **The test is a Welch t on the mean, per row. KS is not used.** In the
   single-trait cells, simACE's SD over seeds is mostly larger than
   xftsim's: the median ratio is 0.94 to 1.18 by statistic, and individual
   cells range from 0.76 to 1.72 (`am/pilot/summary.md`, "Spread").
   - In the two-trait cells, simACE's pinned mate correlations have about
     half xftsim's SE.
   - Both differences are by design. λ = 0.5 gives variable sibship sizes,
     and simACE's Metropolis pinning removes most of the spread of the mate
     correlations.
   - A two-sample KS test would detect such dispersion gaps with no bug
     present. Each row is one p-value in K (ADR 0001 §4).

   Alternatives considered:
   - *KS plus Welch,* as in layer 3. Rejected for the reason above.
   - *Equivalence tests (TOST)* with a stated margin. They need a margin
     with a rationale. Deferred.
     - A pass does not rule out a gap the size of the minimum detectable
       difference (MDD): at that size, a row still passes 20% of the time.
     - The MDDs here are normal approximations, not finite-sample Welch
       power.
     - Phase 2 reports each row's 95% CI for the mean gap next to its p,
       and `power.py` measures real power (decision 8).
5. **M = 1000, with the finite-M gap bounded only as tightly as the pilot
   allows.** Over 40 pilot seeds per cell, the 95% t-CI of Var(A) minus
   the Bulmer prediction, at g = 0 to 6, has its largest endpoint
   (`am/pilot/summary.md`, "Finite M") at:

   | | xftsim M = 250 | M = 1000 | M = 4000 | simACE |
   |---|---|---|---|---|
   | r = 0.5 | 0.0067 | 0.0054 | 0.0054 | 0.0071 |
   | r = 0.8 | 0.0080 | 0.0094 | 0.0082 | 0.0071 |

   - **What the pilot shows.** It shows no trend in M. The M = 250 minus
     M = 4000 sentinel has 95% Welch CIs that cover 0 at every g, at both
     r. At r = 0.8 every xftsim mean is below the prediction, at all M
     (−0.0000 to −0.0046).
   - **What it does not show.** The pilot bounds the M = 1000 gap at about
     0.005 to 0.009 in Var(A), no tighter. It does not prove the gap is
     negligible. Codex notes that the general infinitesimal approximation
     error is O(M^-1/2), with some moment corrections O(M^-1) and no
     universal constant (UNVERIFIED; we have not read the source it cites).
     Phase 2 keeps the M sweep as an INFO sentinel at harness R and
     reports CIs, not p-values.
   - **Inbreeding.** simACE ignores it in the segregation variance (D1;
     the default `off` mode at v0.2.0 and v0.3.0). After 6 generations at N = 10⁴, F is of order
     3·10⁻⁴, a change of about 10⁻⁴ in Var(A). That figure is an
     order-of-magnitude estimate, UNVERIFIED by a run.
   - **Cost per run, 6 jobs in parallel (median wall time; max RSS).** The
     two pilots measured different times on the same host; the later one
     is slower. The cause (host load from other sessions) is UNVERIFIED,
     and `tools/host_check.py` was not run.

     | Configuration | Earlier pilot | Final pilot | Max RSS |
     |---|---|---|---|
     | M = 250 | 11.9 s | 24.7 s | 543 MB |
     | M = 1000 | 16.5 s | 34.2 s | 775 MB |
     | M = 4000 | 35.2 s | 72.8 s | 1472 MB |
     | two traits, M = 1000 | — | 36.5 s | 776 MB |
     | simACE | — | 3.3 s (4.2 s two-trait) | — |

   Alternatives considered:
   - *M = 4000.* It costs about twice as much, and the pilot shows no gain.
   - *M = 250.* It also shows no measurable gap, but its margin is thinner
     as R grows.
6. **Seeds and RNG state.**
   - **Pilot seeds** are 9000001 to 9000040 (GLOSSARY "Pilot seed"); they
     are never scored. The follow-up used 9100001 to 9100200.
   - **xftsim draws meiosis inside numba,** whose RNG state is separate
     from numpy's and kept per thread (`xftsim/reproduce.py`, `_meiosis_i`;
     found by Codex). `am/run_xft.py` therefore seeds numba's state on the
     main thread and in a parallel region, and the pilot runs with
     `NUMBA_NUM_THREADS=1`.
     - The final pilot reproduced all 280 single-trait xftsim records of
       the run before it, exactly.
     - Output at 4 threads differs from output at 1 thread, so the thread
       count belongs in the manifest.
   - **Cells must not share seeds.** xftsim draws E noise and its mating
     noise from one numpy stream in a fixed order, so the same seed at a
     different r reuses the same draws.
     - In the pilot, the per-seed correlation of xftsim `r_FS_L` at g = 3
       is +0.64 between the r = 0 and r = 0.5 cells, and +0.39 between
       r = 0 and r = 0.8. For simACE it is +0.15 and +0.19
       (`am/pilot/cross_cell.md`).
     - Each scored case therefore gets its own seed range. Scored runs may
       still use the same integers on the two sides, because simACE's
       generator is separate.
   - **Fixed before scoring.** Each scored range is fixed in cases.yaml
     before any scored run (ADR 0002 §10).
7. **R.** The approximate MDD at 80% power and α/K = 5·10⁻⁵ (a
   harness-sized K), by row group, over r and g (`am/pilot/summary.md`;
   normal approximation):

   | group, single trait | R = 50 | R = 100 |
   |---|---|---|
   | Var(A) | 0.0070 to 0.0167 | 0.0049 to 0.0118 |
   | Var(L) | 0.0133 to 0.0227 | 0.0094 to 0.0161 |
   | r_mate_L | 0.0042 to 0.0150 | 0.0030 to 0.0106 |
   | r_mate_A | 0.0094 to 0.0164 | 0.0066 to 0.0116 |
   | FS (A, L) | 0.0080 to 0.0147 | 0.0056 to 0.0104 |
   | PO (A, L) | 0.0046 to 0.0093 | 0.0033 to 0.0066 |

   Over all two-trait rows, the MDD at R = 100 is 0.0035 to 0.0200.

   A real AM bug is much larger than these. Two examples:
   - Mating on A instead of L would move `r_mate_A` from r·V/V_L ≈ 0.28
     (single trait, r = 0.5) to r.
   - Pairing that ignored r would leave Var(A) at g = 6 near 0.50, not
     0.68 (r = 0.5) or 0.89 (r = 0.8).

   An error in the two-trait cross structure would also be large. If the
   traits had stayed correlated (simACE's default `rA` 0.5), one pilot seed
   at N 2000 gave r(A1, A2) = 0.614 at g = 3, against 0.188 under the
   matched law.

   **R = 100 per side.** The references are 600 runs for the 6 scored cases
   and 400 for the sentinel's 4 sets.

## Pilot results (40 seeds per side and cell; `am/pilot/summary.md`)

Built by `am/pilot.sh` and `am/summarize.py`. The raw records are in
`am/pilot/*.jsonl`.

Two earlier runs were superseded and are not kept:
- 20 seeds, without numba seeding;
- 40 seeds, with simACE at its default `rA` 0.5 and `rC` 0.3. Trait 1's
  law is the same with or without these defaults, but the passenger was
  not independent as claimed.

**Random mating (check 3).** Every A statistic matches 2φA:

| Statistic | Range | Expected |
|---|---|---|
| Var(A) | 0.4972 to 0.5039 | 0.5 |
| Full-sib A | 0.4957 to 0.5021 | 0.5 |
| Parent–offspring A | 0.4976 to 0.5014 | 0.5 |

The smallest Welch p at r = 0 is 0.00499 (`r_FO_L`, an INFO row, g = 2).

The superseded `rA` 0.5 run had one outlier: `r_FS_L` at g = 3, simACE
0.2418 (z −4.04) against xftsim 0.2555, Welch p 1.04·10⁻⁵. That was below
the pilot's Bonferroni bound of 1.60·10⁻⁴.
- **simACE.** A rerun of simACE alone on 200 fresh seeds
  (`am/followup_fs.py`, `am/pilot/followup_fs.md`, independent traits)
  shows no deviation. It gives `r_FS_L` at g = 3 of 0.2504 (z +0.38), and
  over g = 2 to 6, |z| ≤ 1.71 for the full-sib A, E and L correlations. The
  final pilot, with different simACE draws, has simACE `r_FS_L` 0.2475 to
  0.2511.
- **xftsim.** Its g = 3 excursion (+2.65 z) is shared across its r cells
  (decision 6).

So that row was chance at those draws, and a harness-sized run could fail
on such chance. That is the false-positive rate the power nulls (decision
8) must measure.

**Single-trait AM.** Both sides follow the Bulmer recursion:

| r | Statistic at g = 6 | simACE | xftsim | Predicted | Welch p |
|---|---|---|---|---|---|
| 0.5 | Var(A) | 0.6839 | 0.6812 | 0.6821 | 0.433 |
| 0.5 | Full-sib A | 0.6348 | 0.6335 | 0.6335 | — |
| 0.8 | Var(A) | 0.8880 | 0.8840 | 0.8885 | 0.275 |

- The smallest Welch p is 0.013 at r = 0.5 and 0.0626 at r = 0.8.
- **Mate correlation at r = 0.8.** simACE's realised mate correlation is
  0.7979 to 0.7997 (mean over g 0.7987). xftsim's is 0.7990 to 0.8005.
  That is consistent with the v0.2.0 copula shortfall of simACE #74. By
  the recursion it lowers g = 6 Var(A) by 0.0011.
- Over all 156 single-trait rows the smallest Welch p is 0.00499, against
  a pilot Bonferroni bound of 1.60·10⁻⁴.

**Household C = 0.2, single trait, r = 0 and 0.5** (added 2026-11-05;
`am/pilot/summary.md`, "Household C"). Both sides match the Bulmer
expectation, with C added to the full-sib liability covariance only:

| Statistic | r | simACE | xftsim | Predicted |
|---|---|---|---|---|
| `r_FS_L`, g = 2 to 6 | 0 | 0.4482 to 0.4517 | 0.4477 to 0.4529 | 0.45 |
| `r_PO_L`, g = 2 to 6 | 0 | 0.2472 to 0.2528 | 0.2493 to 0.2518 | 0.25 |
| Var(A), g = 6 | 0.5 | 0.6837 | 0.6794 | 0.6821 |
| `r_FS_L`, g = 6 | 0.5 | 0.5377 | 0.5344 | 0.5347 |

- Over the 104 rows, the smallest Welch p is 1.75·10⁻⁴: `r_FO_L`, an INFO
  row, at r = 0 and g = 5. simACE gives 0.2440 (z −3.67). That is below
  that table's Bonferroni bound of 4.81·10⁻⁴.
- A rerun of simACE alone on 200 fresh seeds (`am/followup_c.py`,
  `am/pilot/followup_c.md`) gives `r_FO_L` z −1.34 to +0.31 over
  g = 2 to 6, so that row was chance.
- The same rerun shows simACE's full-sib `r_FS_L` slightly low at every g:
  0.4478 to 0.4497 against 0.45, z −2.23 to −0.26. That is a gap of about
  0.001, far below the FS MDD at R = 100 (decision 7). It is recorded here,
  not chased. Layer 1 scores this law conditionally on the pedigree.

**Two-trait AM, r = 0.2 and 0.4.** simACE and xftsim agree with each other
and with `bulmer2`. At r = 0.4, g = 6:

| Statistic | simACE | xftsim | Predicted |
|---|---|---|---|
| Var(A) | 0.6638 | 0.6609 | 0.6636 |
| r(A1, A2) | 0.2489 | 0.2463 | 0.2465 |
| r(L1, L2) | 0.1420 | 0.1398 | 0.1406 |

- All four mate correlations on L are within 0.0036 of r on both sides.
- Over the 166 two-trait rows, the smallest Welch p is 0.028 (`r_FS_L`,
  r = 0.2, g = 4), against a Bonferroni bound of 3.01·10⁻⁴.
- The largest |z| against `bulmer2` is 2.37 for simACE (`var_A_t2`,
  r = 0.2, g = 6) and 2.28 for xftsim (`var_A`, r = 0.4, g = 6).

## Decision 8: wiring into the harness

- **A committed xftsim reference sample.** xftsim's output does not depend
  on simACE. Each case's R xftsim runs are made once, by
  `pixi run -e xftsim xft-reference`, and committed as Parts (ADR 0002 §9).
  Each Part has a manifest, and the case is ERROR if any field differs from
  the current value:
  - the xftsim commit;
  - the `xftsim` env's lock entries;
  - the config;
  - the hashes of `am/stats.py` (the statistics, and `flatten`),
    `am/run_xft.py` (the driver that builds effects, seeds the RNGs and
    pairs) and `am/reference.py` (each run's argv and thread settings);
    the Bulmer predictions live in `am/bulmer.py`, outside the hashes,
    since no sample is built from them;
  - `NUMBA_NUM_THREADS`;
  - the seed range.

  `pixi run harness` reruns only simACE (a few seconds per seed) and scores
  it against the committed sample. This keeps `harness` in the default env
  and keeps xftsim off its critical path.
  - The reference's own sampling error stays in the Welch SE, so the test
    is calibrated over draws of the reference.
  - It is not calibrated conditional on the one committed sample. Every
    later run reuses that sample, so a sample that happens to be extreme
    biases all of them the same way.
  - The reference is regenerated only when a pin changes, never after a
    result is seen.
- **Cases.** Each case has `oracle: xftsim`, R = 100, its own seed range,
  and `rA = rC = rE = 0` stated explicitly.
  - `am_r0` (check 3), `am_r05` and `am_r08`, single trait: 42 Welch rows
    each. That is Var(A) and Var(L) at g = 1 to 6, and 6 statistics at
    g = 2 to 6.
  - `am2_r02` and `am2_r04`, two traits: 114 Welch rows each:
    - both traits' 42-row sets (84);
    - r(A1, A2) and r(L1, L2) at g = 2 to 6 (10);
    - the cross-trait mate correlations on L, in both orders, at g = 2 to 6
      (10);
    - the cross-trait mate correlations on A, in both orders, at g = 2 to 6
      (10).
  - `am2_r0401`, two traits with r1 0.4 and r2 0.1 (rank 1): the same 114
    Welch rows.
  - `amc_r05`, single trait with household C 0.2: 42 Welch rows.
  - Single-trait rows add 168 p-values to K, and the first two two-trait
    cases add 228: 396 in all. With the 835 of the last full run, α/K fell
    from 5.99·10⁻⁵ to 0.05/1231 = 4.06·10⁻⁵ (run of 2026-11-05).
    `am2_r0401` adds 114 more. With main merged (layer 4, layer 3b
    rescored), the run of 2026-11-06 against simace v0.3.0 has K = 1353,
    α/K 3.7·10⁻⁵.
  - INFO rows: the MO and FO rows, and one `bulmer:<row>` row per scored
    row, a one-sample t of simACE against the Bulmer expectation. As
    built, no KS is computed at all.
  - An INFO sentinel case compares xftsim at M = 250 and M = 4000, at
    r = 0.5 and 0.8, and reports the mean gap and its Welch t by generation
    (the 95% CI follows from them). Agreement does not
    prove closeness to the large-M limit; the CI bounds the gap only as
    tightly as R allows.
- **simACE pin: option (b), decided.** Layer 5 scored simace v0.2.0 at
  first (run of 2026-11-05). The pin moved to v0.3.0, which ships ADR 0032,
  with the merge of main on 2026-11-06, and the run of 2026-11-06 scores it.
  - At r = 0.8 its mean `r_mate_L` over g is 0.7987, against xftsim's
    0.7997. That gap is below the `r_mate_L` MDD at R = 100 (0.0030 at its
    smallest), so a FAIL is possible but not likely.
  - If a row fails and the cause is simACE #74:
    - differences.md gains a D entry for it;
    - the rows it affects become INFO (ADR 0001, Consequences);
    - status.md records the row.
  - Dropping r = 0.8 after seeing its shortfall would be selection after
    the fact, so it is not an option.
  - When the simace pin moves to a release that ships ADR 0032, the
    reference stays valid; only simACE's side changes. The D entry is then
    re-tested.
  - The two-trait path already pins its four mate correlations by
    Metropolis at v0.2.0, so #74 does not apply there.
- **Power.** ADR 0001 §6 puts mutations on the oracle side. Here the oracle
  side is the reference sample, so `power.py` builds mutated xftsim
  references and scores the pinned simACE against them:
  - xftsim at 0.9·r;
  - xftsim mating on the additive component instead of the phenotype;
  - xftsim with r = 0;
  - for two traits, xftsim mating on y1 alone.

  There are two nulls:
  - simACE on fresh seeds against the committed reference (`null`), the
    false-positive rate given that sample;
  - xftsim against xftsim (`null_xft`), a second xftsim sample in simACE's
    place. It needs R xftsim runs per rep and was not run.

  Run on am_r05 and am2_r04 at 10 reps against simace v0.2.0: every
  mutation rejected 10/10; `null` rejected 0/10 (am_r05) and 1/10
  (am2_r04, one row, `r_FS_L_g4`). Rerun against v0.3.0 on 2026-11-06
  (results/power_am.md): every mutation 10/10; `null` 0/10 (am_r05) and
  1/10 (am2_r04, one row, `r_FS_A_t2_g3`).
- **status.md G2b** becomes Closed with the cases' results, for the scope
  tested: single-trait AM (C = 0, and household C 0.2 at r = 0.5), and
  two-trait AM with every mate correlation equal (C = 0), with independent
  traits. differences.md gains a D entry for the matched configurations,
  covering:
  - C only as a household draw by mother, single trait, no half-sibs;
  - the family-structure difference behind decision 4;
  - the finite-M CI of decision 5;
  - the restriction of two-trait AM to rank-1 laws: R_mf = r·11ᵀ, and
    v·vᵀ with unequal r through our pairing regime.

## Consequences

- A FAIL in layer 5 means that simACE's population moments under AM differ
  from xftsim's at the matched configuration. It does not say which side is
  wrong. The Bulmer columns are a diagnostic, not a referee, since they are
  approximations themselves.
- These stay unchecked; each needs its own matched law:
  - two-trait AM with C > 0, and the half-sib split of C (xftsim's regimes
    make no half-sibs);
  - two-trait AM with any cross term other than √(r1·r2), including
    simACE's derived one (`_cross_am_matrix`), which is 0 for independent
    traits;
  - correlated traits;
  - AM combined with C transmission (rwaples/simACE-validate#5).
- xftsim's environment is Python 3.9 with numba 0.56.4 and numpy 1.23 (pixi
  env `xftsim`, no default feature). It is installed, never vendored (ADR
  0001 §7), and no xftsim code is copied: `Rank1AssortativeMatingRegime` is
  our own, written clean-room against xftsim's API (decision 2). The env
  holds only what the xftsim runs need (no ruff), since its lock entries are
  hashed into every reference.

## Open items

- **rwaples/simACE-validate#5: C transmission against xftsim** (ADR 0025;
  check 2 of #4, promoted). It chains the parents' c through a vertical
  component, plus the household term above. Not built, pending approval.
  Constraints found here that apply to it:
  - xftsim's built-in regimes make no half-sibs (decision 2), so the MHS
    (1 − b_f²) and PHS (b_f²) rows need a regime of our own. That would be
    our code on xftsim's side, which weakens independence.
  - Its reference samples would follow decision 8's manifest design. The
    household construction already lives in `am/run_xft.py`.
