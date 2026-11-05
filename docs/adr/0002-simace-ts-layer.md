# ADR 0002: simACE-ts scores simACE's downstream stages given a gene-drop A

**Status:** accepted
**Date:** 2026-11-03

Amends ADR 0001 §2 (the test), §4 (what K counts) and §6 (power) for
layer 4 only. Resolves the open questions of rwaples/simACE-validate#2.
Revised the same day after an independent Codex review: the first draft
scored against 2φ behind a gate, which mixed the pipeline's genetic law
with simACE's correctness.

## Context

`tskit_pipeline/` drops SimHumanity p2 founder genomes through a simACE
pedigree with msprime and a HapMapII GRCh38 map, assigns a finite causal
architecture with tstrait, and replaces A1 with the rescaled genetic value.
simACE's `cohort` and `analyze` then run on that pedigree.

Every value of A1 in this layer comes from msprime, the graft, tstrait and
the rescale, not from simACE. Comparing gene-drop A with 2φ therefore tests
the pipeline and the SimHumanity panel, not simACE. Its expectation is not
even exactly 2φ: founder relatedness and structure in the panel, LD among
causal sites, and inbreeding (2φ/√((1+F_i)(1+F_j))) all move it, by amounts
nobody has measured. Layer 2 already tests that simACE's links form the
pedigree they claim. What only this layer can test is that simACE's
downstream stages give the right answer when A is a finite-locus value
with LD rather than a Gaussian draw.

One full-scale seed (N = 35,951 founders, 10 generations) costs CPU-hours
(`tskit_pipeline/architecture-sweep.md`: about 2 CPU-hours for the drop
alone at `G_pheno` = 3). ADR 0001's rank test needs M oracle draws per
pedigree, which is out of reach.

## Decision

1. **Subject.** Layer 4 (simACE-ts) scores simACE's `cohort` and `analyze`
   output given the gene-drop A it was handed. The pipeline supplies A. The
   harness never scores A itself against 2φ.
2. **Realised truth.** For each scored `report.yaml` field, the harness
   computes the expectation on the same rows, in the same scope, from the
   realised gene-drop A and the configured C and E law:
   - a liability pair correlation over a class's pairs (i, j):
     (Cov(A_i, A_j) + C1·s) / √((Var A_i + C1 + E1)(Var A_j + C1 + E1)),
     with the realised A's population moments over exactly those pairs and
     s the realised fraction of pairs sharing a `household_id`;
   - a midparent liability slope: (Cov(A_child, A_mid) + C1·s_cm) /
     (Var A_mid + (C1·(1 + s_mf) + E1)/2), with s_mf the fraction of parent
     pairs from one household and s_cm that of child–parent pairs;
   - in both, C1 and E1 are the configured variances, not the realised
     ones: a realised V_P would absorb a wrong C or E variance in simACE;
   - an A correlation or A regression slope and R²: none is computed. These
     are statistics of the supplied A alone, so the harness recomputes them
     exactly (item 3);
   - Falconer and the twin estimators: absent, since the case has no MZ
     twins (item 6).
   Each field's target is derived on its own and listed with its
   `report.yaml` path in the case's statistic schema. A field without a
   derived target is not scored.
3. **Exact recomputation.** Fields that are deterministic functions of the
   supplied pedigree and cohort (A correlations, A slopes and R², realised
   variances, pair counts) are recomputed by harness code on the same rows
   and must match to float tolerance. A mismatch is FAIL. These rows carry
   no p-value and do not count in K.
4. **Scored rows.** For each field with a realised truth (item 2), the
   per-seed gap (simACE's value − realised truth) gets a t-test over seeds
   that its mean is 0 (one p-value, ADR 0001 §4). Liability correlations
   come from the analysis sample. The quality estimators come from the
   recorded pedigree (`simace/analysis/analyze.py:111-150` at v0.2.0). Each
   row names its scope.
5. **INFO rows.** These are reported, never scored:
   - quantities that assume Gaussian liability by simACE's design: threshold
     prevalence (`ndtri(1 − K)`, `simace/phenotype/models/_prevalence.py:130`
     at v0.2.0), the ADuLT and LTM onset mapping through a normal survival
     function (`simace/phenotype/models/adult.py:220-225`), tetrachoric
     correlations, and observed-scale h², interpreted as Gaussian-liability
     h². Explained in differences.md D3.
   - the pipeline diagnostic: gene-drop A pair correlation − 2φ by class
     and generation (check 1), and the gene-drop variance by generation
     before the rescale (LD decay and drift).
   - the architecture sweep (check 3) and the layer-2 cross-check on the
     same pedigree (check 4).
6. **Case.** `p_mztwin` = 0. The pipeline drops co-twins as separate
   offspring, so it refuses a pedigree with twins. The case also uses
   pass-through ascertainment, no C transmission (`Ctm-*` = 0) and no
   assortative mating, so that the analysis sample is the phenotyped
   population and γ is household membership alone. Each restriction can be
   lifted later with a derived target.
7. **Pipeline changes.**
   - Every generation is a sample (`--g-pheno` = G_ped in
     `build_pedigree_tables.py`), so no row mixes gene-drop and Gaussian A.
   - The rescale is per generation (Var = A1 in each), so simACE's
     downstream stages see the configured scale in every generation.
8. **Seeds.** Each seed gets its own pedigree, drop and causal
   architecture. No `drop_from` reuse, since seeds that share a drop are
   not independent. Seeds share the p2 founder panel, so every row holds
   given that panel. Each stage and chromosome takes its seed from its own
   namespace, `SeedSequence([s, stage, chrom]).generate_state(1)[0]`, not
   from integer offsets of s (which already overlap: chromosome 1's
   s + 100 equals trait 2's phenotype seed, `simace/phenotype/runner.py:148`).
9. **Parts.** Seeds run on GenomeDK. Each writes one Part: simACE's value
   and the realised truth for every scored row, the exact-recomputation
   results, the INFO statistics, and a manifest. The manifest holds the
   case config, the statistic schema, simace's version and git ref, the
   simACE-validate commit, the SimHumanity data pins (below), the seed and
   the chromosomes run. Parts are committed. `pixi run harness` scores a
   Case's Parts like any Case. It makes the Case ERROR when:
   - a seed in 0..R−1 has no Part, or two Parts share a seed;
   - any manifest field differs from the current case config, schema,
     simace ref or data pins;
   - a scored value or truth is not finite, or a scored row is missing in
     any seed.
   K counts the scored rows (item 4) fixed by the schema, not the rows that
   happen to be present. Exact-recomputation and INFO rows are not in K.
10. **Sizing and power.** Three pilot seeds (from a separate range, never
    scored) measure CPU-hours and peak memory. R is then fixed from the
    budget, at most 99, and written into cases.yaml before any scored seed
    runs. It is never revised after scored results are seen. A measured
    power study (ADR 0001 §6) would need repeated full sets of R seeds,
    which is out of reach. Instead, results.md reports each row's minimum
    detectable mean gap at α/K, computed from its observed SE across seeds.

    **Fixed 2026-11-05: R = 50** (seeds 9000 to 9049), from the pilot seeds
    990000 to 990002. Each took 32 to 33 min on 8 CPUs, about 4.0 CPU-hours
    summed over stages (drops about 2.0, genetic values about 2.0, of which
    the 1M-site sweep architecture took 1.66), with a job memory peak of
    36.0 GB. At α/K = 5.91e-5 (K = 846), R = 50 detects a mean gap of
    0.74 per-seed SDs at 80% power, against 1.01 at R = 30 and 0.51 at
    R = 99. The pilot's gaps were not used to choose R beyond that spread.

## Data

| Input | Pin |
|---|---|
| `simhumanity_trees_RO.zip` | http://benhaller.com/slim/simhumanity_trees_RO.zip, 9,573,536,479 bytes, sha256 `003fa571a76ea646d50f21cfc3be383922a78aa661baafc5d86cdd00f42193e8` |
| Recombination maps | github.com/bhaller/SimHumanity at `aa2242e`, `stdpopsim extraction/extracted/` |
| Preprocessed p2 | rebuilt on GenomeDK with `--pop-name p2`, chromosomes 1 to 22; `all_chroms.trees` fingerprint must equal the laptop's `ddea736699a31c66ade5a140051c3385` |

## Consequences

- A FAIL in layer 4 means that simACE's downstream stages are wrong given a
  finite-locus A of the stated variance on this panel. It says nothing
  about whether gene-drop A follows 2φ. That is the INFO diagnostic.
- The issue's "done when" asked for check (1) as a calibrated row. Check
  (1) is INFO here, because its expectation has no derivation. Check (2) is
  the calibrated row instead.
- The layer depends on simACE more than layers 1 to 3. Its Oracle consumes
  simACE's pedigree and feeds simACE's own stages. README.md and the
  status.md rows say so.
- A simace bump makes every Part stale. A full run then needs R seeds on
  GenomeDK before `results/` can be rewritten.
- MZ twins, ascertainment, C transmission and assortative mating are
  unchecked in this layer until each has a derived target.
- Scored rows are means only. As in layer 2, spread is INFO.
