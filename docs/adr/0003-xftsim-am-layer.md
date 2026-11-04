# ADR 0003: Assortative mating is scored as a two-sample comparison against xftsim

**Status:** proposed (draft, phase 1 of rwaples/simACE-validate#4). Decisions
1 to 7 are locked for the pilot. Decision 8 (the harness wiring) awaits
approval.
**Date:** 2026-11-04

Amends ADR 0001 §1 (conditioning) and §2 (the test) for layer 5 only.
Answers the open questions of rwaples/simACE-validate#4.

## Context

G2b (status.md) is open because layer 1 has no law for assortative mating
(AM). ADR 0001 conditions on the pedigree simACE produced. Under AM that
pedigree was chosen on the liabilities. So, given the pedigree, simACE's A
is not drawn by midparent plus Mendelian sampling, and the layer-1 oracle
would have to condition on the pairing too. No tractable law does that.

[xftsim](https://github.com/rborder/xftsim) (GPLv3; v0.3.1, commit
`8a4cdf0`) is a forward-time simulator that shares no code with simACE. It
has finite loci and its own AM regime, so it can serve as an empirical
reference at the population level.

The two AM laws at the pinned versions:

- **simACE v0.2.0** (single-trait `assort1`; installed
  `simace/simulation/simulate.py:699-726`). It ranks mating slots by
  liability A + C + E in each sex, draws M bivariate normals with
  correlation r, and pairs the female and male slots whose ranks match each
  draw's ranks: a Gaussian copula of correlation r on the liability ranks.
  It then swaps out duplicate couples (`:728-749`). A parent with k partners
  fills k slots. In v0.2.0 the realised correlation falls slightly below r
  (simACE #74; `plans/issue-69-mate-correlation-test.md` reports a bias of
  −0.0017 to −0.0029 at N 2k, λ 0.5). simACE dev now pins it with a 1-D
  Metropolis step (simACE ADR 0032, `_metropolis_single`), which no release
  ships yet.
- **xftsim `LinearAssortativeMatingRegime`** with one component
  (`xftsim/mate.py:843-887` at `8a4cdf0`). It splits the population at
  random into two halves, standardises the phenotype y, and builds the
  score √r·y + √(1−r)·ε with independent noise ε in each half. Each half is
  sorted by its score and the two are paired rank to rank. For Gaussian y
  and large n the two scores are paired almost exactly, so
  (y_f, y_m) ≈ BVN(r): the same Gaussian copula of correlation r on the
  phenotype.

So the two laws coincide in the large-n limit only under Gaussian
liability, and only when the mating variable is the same: the full
liability, noise included. xftsim's phenotype is A + E with no C, so the
matched configuration sets C = 0.

The match is approximate, not exact, in several ways:

- At fixed finite M, xftsim's A is not exactly Gaussian. Pairing on noisy
  scores and simACE's rank copula then need not give identical laws.
- xftsim standardises with the population SD but takes its target from
  `np.cov` (sample variance). With one component, this gives
  R = r(n−1)/n: 0.79984 at r = 0.8 and n = 5000 (pointed out by Codex).
- Finite rank matching lowers both simulators' realised correlation
  slightly.
- simACE's slots (one per partnership) and dedup swaps, and the two
  simulators' different family structures, can move the means at finite
  N.

The pilot measures the combined effect. It supports approximate agreement
within the measured MDDs, not equality.

## Decision

1. **Subject: population moments, two-sample.** Layer 5 compares
   per-seed, per-generation statistics of simACE's pedigrees with those of
   independent xftsim runs at a matched configuration. Neither side is
   conditioned on the other's pedigree. The infinitesimal Bulmer recursion
   (`am/stats.py::bulmer`) is reported beside each row as INFO.
   Alternatives considered:
   - *Score against the Bulmer recursion alone* (one-sample). It is a
     first-order law: it ignores finite-N drift, the copula's departure
     from exact bivariate normality at finite M, and dedup. It is theory
     that simACE's own docs share, not independent code. At harness R its
     approximation error could become detectable. Kept as an INFO column.
   - *Extend the layer-1 oracle to AM.* That needs the law of A given a
     pedigree that was itself selected on liability. Rejected as
     intractable.
2. **Matched configuration.** Single trait, h² = 0.5, C = 0, E = 0.5,
   mating on the full liability, no MZ twins, N = 10⁴ per generation, 6
   rounds of mating from founders, and r ∈ {0, 0.5, 0.8}. r = 0 is check 3
   of the issue, the random-mating baseline, run through the same xftsim
   code path. The details:
   - **simACE:** `G_ped = G_sim = 6` and `mating_lambda` 0.5 (the default).
     Trait 2 is an independent passenger. Recorded generation k follows
     k + 1 matings (ADR 0001 §3), so it is reported as g = k + 1.
   - **xftsim:** M = 1000 unlinked loci (`constant_map_from_haplotypes`,
     p = 0.5), one mate each, two offspring per couple, and fixed-scale
     effects (decision 3).
   - **Compared statistics** (`am/stats.py::generation_stats`, one function
     run on both sides): Var(A) and Var(L) for g = 1 to 6; for g = 2 to 6,
     the mate correlation on L and on A (over distinct couples), the
     full-sib correlation of A and L, and the parent–offspring correlation
     of A and L (mother and father pooled).
   - **Not compared:**
     - MO and FO separately. xftsim's `sex_aware=False` makes "maternal"
       a random label, so those rows are INFO.
     - Half-sibs. The xftsim configuration makes none.
     - Generation 0 of simACE, which is not recorded.
   Alternatives considered:
   - *Map C onto an xftsim vertical or horizontal component.* That would
     add a second law to match, and C transmission is already Closed
     (G2a). Deferred to a later INFO check (check 2 of the issue).
   - *Match family structure* by setting xftsim's mates per female to a
     ZTP. The population moments do not depend on it to first order, and it
     would make xftsim's pairing law depart from simACE's slot model in a
     different way. Rejected: only moments are compared, never pedigrees.
3. **xftsim effects are fixed-scale.** β ~ N(0, 1) is rescaled so that
   mean(β²) = h² exactly, through `AdditiveEffects(scaled=True,
   standardized=True)`. The founders' expected Var(A) is then h² given the
   allele frequencies. `GCTA_Architecture` leaves the scale random, with a
   relative SD of √(2/M) = 4.5% at M = 1000, which would dominate every
   Var(A) row's spread. Alternative: keep GCTA and score Var(A)/Var(A₀).
   Rejected, because simACE's founders are not recorded, so the ratio
   has no simACE counterpart.
4. **The test is a Welch t on the mean, per row. KS is not used.** The
   pilot's SD over seeds is mostly larger for simACE than for xftsim
   (median ratio 0.99 to 1.46 by statistic; 1.46 for full-sib A,
   `am/pilot/summary.md`, "Spread"). That is by design: λ = 0.5 gives variable sibship sizes and
   fewer full-sib pairs than xftsim's fixed pairs of two. A two-sample KS
   test would detect that dispersion gap with no bug present. Each row is
   one p-value in K (ADR 0001 §4). Alternatives considered:
   - *KS plus Welch,* as in layer 3. Rejected for the reason above.
   - *Equivalence tests (TOST)* with a stated margin. They need a margin
     with a rationale. Deferred. A pass does not rule out a gap the size
     of the minimum detectable difference (MDD): at that size, a row still
     passes 20% of the time. The MDDs here are normal approximations, not
     finite-sample Welch power. Phase 2 reports each row's 95% CI for the
     mean gap next to its p, and `power.py` measures real power (decision
     8).
5. **M = 1000, with the finite-M gap bounded only as tightly as the pilot
   allows.** Over 40 pilot seeds per cell, the 95% t-CI of xftsim Var(A)
   minus the Bulmer prediction, at g = 0 to 6, has its largest endpoint
   (`am/pilot/summary.md`, "Finite M") at:

   | | M = 250 | M = 1000 | M = 4000 | simACE |
   |---|---|---|---|---|
   | r = 0.5 | 0.0067 | 0.0054 | 0.0054 | 0.0061 |
   | r = 0.8 | 0.0080 | 0.0094 | 0.0082 | 0.0144 |

   - **What the pilot shows.** It shows no trend in M. The M = 250 minus
     M = 4000 sentinel has 95% Welch CIs that cover 0 at every g, at both
     r. At r = 0.8 every xftsim mean is below the prediction, at all M
     (−0.0000 to −0.0046). simACE is too from g = 4 on. The Bulmer
     recursion ignores finite N, so this may be a shared departure from it.
     It is not a difference between the simulators.
   - **What it does not show.** The pilot bounds the M = 1000 gap at about
     0.005 to 0.009 in Var(A), no tighter. It does not prove the gap is
     negligible. Codex notes that the general infinitesimal approximation
     error is O(M^-1/2), with some moment corrections O(M^-1) and no
     universal constant (UNVERIFIED; we have not read the source it
     cites). Phase 2 keeps the M sweep as an INFO sentinel at harness R
     and reports CIs, not p-values.
   - **Inbreeding.** simACE ignores it in the segregation variance (D1;
     off mode at v0.2.0). After 6 generations at N = 10⁴, F is of order
     3·10⁻⁴, a change of about 10⁻⁴ in Var(A). That figure is an
     order-of-magnitude estimate, UNVERIFIED by a run.

   Cost per run, at 6 jobs in parallel (median wall time and max RSS):
   - M = 250: 11.9 s and 543 MB (80 runs);
   - M = 1000: 16.5 s and 775 MB (120 runs);
   - M = 4000: 35.2 s and 1472 MB (80 runs).

   Alternatives considered:
   - *M = 4000.* It costs about twice as much, and the pilot shows no gain.
   - *M = 250.* It also shows no measurable gap, but its margin is thinner
     as R grows.
6. **Seeds and RNG state.**
   - **Pilot seeds** are 9000001 to 9000040 (GLOSSARY "Pilot seed"); they
     are never scored. The follow-up below used 9100001 to 9100200.
   - **xftsim draws meiosis inside numba,** whose RNG state is separate from
     numpy's and kept per thread (`xftsim/reproduce.py`, `_meiosis_i`; found
     by Codex). `am/run_xft.py` therefore seeds numba's state on the main
     thread and in a parallel region, and the pilot runs with
     `NUMBA_NUM_THREADS=1`. Two runs of the same seed then gave identical
     output, at 1 thread and also at 4 (one check each). Without this
     seeding, a reference sample cannot be regenerated from its seeds.
   - **Cells must not share seeds.** xftsim draws E noise and its mating
     noise from one numpy stream in a fixed order, so the same seed at a
     different r reuses the same draws. In the pilot, the per-seed
     correlation of xftsim `r_FS_L` between the r = 0 cell and the
     r = 0.5 cell is +0.64 at g = 3; between r = 0 and r = 0.8 it is
     +0.39. For simACE it is +0.05 and +0.04 (`am/pilot/cross_cell.md`).
     Each scored case therefore gets its own seed range. Scored runs may
     still use the same integers on the two sides, because simACE's
     generator is separate.
   - **Fixed before scoring.** Each scored range is fixed in cases.yaml
     before any scored run (ADR 0002 §10).
7. **R.** The approximate MDD at 80% power and α/K = 5·10⁻⁵ (a
   harness-sized K), by row group, over r and g (`am/pilot/summary.md`, "MDD
   range by row group"; normal approximation):

   | group | R = 50 | R = 100 |
   |---|---|---|
   | Var(A) | 0.0074 to 0.0200 | 0.0053 to 0.0142 |
   | Var(L) | 0.0129 to 0.0243 | 0.0091 to 0.0172 |
   | r_mate_L | 0.0042 to 0.0145 | 0.0030 to 0.0102 |
   | r_mate_A | 0.0096 to 0.0156 | 0.0068 to 0.0111 |
   | FS (A, L) | 0.0087 to 0.0145 | 0.0062 to 0.0102 |
   | PO (A, L) | 0.0055 to 0.0101 | 0.0039 to 0.0071 |

   A real AM bug is much larger than these. Two examples:
   - Mating on A instead of L would move `r_mate_A` from r·V/V_L ≈ 0.28 to
     r.
   - Pairing that ignored r would leave Var(A) at g = 6 near 0.50, not
     0.68 (r = 0.5) or 0.89 (r = 0.8).

   **Proposed R = 100 per side.** One xftsim reference set at 3 values of
   r is 300 runs. At 16.5 s per run that is about 83 minutes of job time,
   or about 14 minutes at 6 jobs in parallel. Both figures are UNVERIFIED
   extrapolations from the pilot medians, which were measured with 6 jobs
   running.

## Pilot results (40 seeds per side and cell; `am/pilot/summary.md`)

Built by `am/pilot.sh` and `am/summarize.py`. The raw records are in
`am/pilot/*.jsonl`. A first run of 20 seeds, without numba seeding, was
superseded by this one and is not kept.

- **Random mating (check 3).** Every A statistic matches 2φA:
  - Var(A) 0.4967 to 0.5023, against 0.5;
  - full-sib A 0.4957 to 0.5021;
  - parent–offspring A 0.4972 to 0.5008.

  One row stands out: `r_FS_L` at g = 3. simACE gives 0.2418 (z −4.04
  against 0.25) and xftsim 0.2555 (z +2.65). Welch p is 1.04·10⁻⁵, below
  the pilot's Bonferroni bound of 1.60·10⁻⁴ over 312 p-values. The two
  sides were checked separately:
  - **simACE.** A rerun of simACE alone on 200 fresh seeds
    (`am/followup_fs.py`, `am/pilot/followup_fs.md`) gives g = 3 `r_FS_L`
    0.2496 (z −0.39). Over g = 2 to 6, full-sib E correlations are
    −0.0012 to +0.0004 (|z| ≤ 1.61), and `r_FS_L` z is −2.18 to −0.12.
    The pilot seeds were a chance excursion, as in F1.
  - **xftsim.** Its g = 3 excursion is shared across its r cells
    (decision 6).

  At the pinned pilot seeds, then, a harness-sized run would have failed
  on chance. That is the false-positive rate a power null (decision 8)
  must measure.
- **AM, r = 0.5.** Both sides follow the Bulmer recursion:
  - Var(A) at g = 6: simACE 0.6825 and xftsim 0.6812, against 0.6821
    predicted (Welch p 0.731);
  - full-sib A: 0.6322 and 0.6335, against 0.6335.

  The smallest Welch p at r = 0.5 is 0.0112.
- **AM, r = 0.8.**
  - **Mate correlation.** simACE's realised mate correlation is 0.7978 to
    0.8004 (mean over g 0.7987; z −2.31 to +0.50 against 0.8). xftsim's is
    0.7990 to 0.8005. This is consistent with the v0.2.0 copula shortfall
    of simACE #74. By the recursion, 0.7987 lowers g = 6 Var(A) by 0.0011
    (0.8874 against 0.8885).
  - **Moments.** Var(A) at g = 6 is 0.8820 for simACE and 0.8840 for
    xftsim, against 0.8885 (Welch p 0.665). The 20-seed run's `r_FS_A` at
    g = 6 signal (Welch p 0.0014) did not hold at 40 seeds: 0.7134 against
    0.7164, Welch p 0.141. The smallest Welch p at r = 0.8 is 0.00388
    (`var_L`, g = 1).
- **Spread.** The SD over seeds, simACE over xftsim, has a median of 0.99
  to 1.46 by statistic, and 1.46 for full-sib A ("Spread").
- **Cost.** simACE takes 1.6 s per run. xftsim's costs are under
  decision 5.

## Decision 8, proposed and not built: wiring into the harness

- **A committed xftsim reference sample.** xftsim's output does not
  depend on simACE. Each matched configuration's R xftsim runs are made
  once, by `pixi run -e xftsim xft-reference`, and committed as Parts
  (ADR 0002 §9). Each Part has a manifest, and the case is ERROR if any
  field differs from the current value:
  - the xftsim commit;
  - the `xftsim` env's lock entries;
  - the config;
  - the hashes of `am/stats.py` and `am/run_xft.py` (the driver that builds
    effects and seeds the RNGs);
  - `NUMBA_NUM_THREADS`;
  - the seed range.
  `pixi run harness` reruns only simACE (about 2 s per seed) and scores it
  against the committed sample. This keeps `harness` in the default env
  and keeps xftsim off its critical path.
  - The reference's own sampling error stays in the Welch SE, so the test
    is calibrated over draws of the reference.
  - It is not calibrated conditional on the one committed sample. Every
    later run reuses that sample, so a sample that happens to be extreme
    biases all of them the same way.
  - The reference is regenerated only when a pin changes, never after a
    result is seen.
- **Cases.** `am_r0` (check 3), `am_r05` and `am_r08` in cases.yaml.
  Each has `oracle: xftsim`, R = 100, and the scored statistics of
  decision 2.
  - That is 42 Welch rows per case (Var(A) and Var(L) at g = 1 to 6, and
    6 statistics at g = 2 to 6), so 126 p-values in K.
  - The Bulmer z, the MO and FO rows, and the KS p are INFO.
  - An INFO sentinel case compares xftsim at M = 250 and M = 4000, at
    r = 0.5 and 0.8, and reports the 95% CI of the gap. Agreement does not
    prove closeness to the large-M limit. The CI bounds the gap only as
    tightly as R allows.
- **simACE pin.** At r = 0.8 against v0.2.0, simACE's mean `r_mate_L`
  over g is 0.7987, against xftsim's 0.7997. That gap is below the
  `r_mate_L` MDD at R = 100 (0.0030 at its smallest). It would be the
  known simACE #74. A FAIL there is possible, not likely. Options:
  - (a) wait for a simace release that ships ADR 0032 (recommended; D1
    also waits for one);
  - (b) score now. If a row fails, record #74 as a D entry and make the
    affected rows INFO (ADR 0001, Consequences). Dropping r = 0.8 after
    seeing its shortfall would be selection after the fact, so that option
    is not offered.
  With (a), Metropolis pinning makes simACE's `r_mate_L` nearly
  deterministic. Its spread then differs from xftsim's by design, which a
  Welch t on the mean already tolerates.
- **Power.** ADR 0001 §6 puts mutations on the oracle side. Here the
  oracle side is the reference sample, so `power.py` builds mutated xftsim
  references and scores the pinned simACE against them:
  - xftsim at 0.9·r;
  - xftsim mating on the additive component instead of the phenotype;
  - xftsim with r = 0.
  There are two nulls:
  - simACE against a second reference sample from a disjoint seed range,
    which includes simACE's own spread;
  - xftsim against xftsim.
- **status.md G2b** becomes Closed with the case's results, for the scope
  tested: single-trait AM with C = 0. differences.md gains a D entry for
  the matched configuration, covering C = 0, the family-structure
  difference behind decision 4, and the finite-M CI of decision 5.

## Consequences

- A FAIL in layer 5 means that simACE's population moments under AM
  differ from xftsim's at the matched configuration. It does not say
  which side is wrong. The Bulmer column is a diagnostic, not a referee,
  since it is an approximation itself.
- AM with C > 0, two-trait AM (`assort2`, `assort_matrix`), and AM
  combined with C transmission stay unchecked. Each needs its own matched
  law.
- xftsim's environment is Python 3.9 with numba 0.56.4 and numpy 1.23
  (pixi env `xftsim`, no default feature). It is installed, never
  vendored (ADR 0001 §7).
