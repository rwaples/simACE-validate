# simACE against independent simulators

Runs [simACE](https://github.com/rwaples/simACE) and checks what it
simulates against code simACE does not share:

- an R liability oracle
- [msprime](https://tskit.dev/msprime/) gene-drop
- [ped-sim](https://github.com/williamslab/ped-sim)
- [SLiM](https://messerlab.org/slim/)
- [xftsim](https://github.com/border-lab/xftsim), for assortative mating

The latest results are in [`results/results.md`](results/results.md), and
the harness's measured power is in [`results/power.md`](results/power.md)
(layer 1), [`results/power_std_mating.md`](results/power_std_mating.md) and
[`results/power_std_structure.md`](results/power_std_structure.md) (layer
3b, simace v0.3.0; each merged from one `--parts` run per mutation) and
[`results/power_am.md`](results/power_am.md) (layer 5).
Every known model gap is explained in [`differences.md`](differences.md),
and where each non-PASS row stands is in [`status.md`](status.md).

simACE's output is random, so a check is a calibrated statistical test,
not a tolerance ([ADR 0001](docs/adr/0001-conditional-calibration.md)).
That is the main difference from
[pedsum-validate](https://github.com/rwaples/pedsum-validate) and
[pg-phenotype-validate](https://github.com/rwaples/pg-phenotype-validate).

## Run it

With [pixi](https://pixi.sh). `pixi.toml` pins simace to the `v0.3.0`
family tag.

```bash
pixi run harness                          # every case, to results/ (builds nadiv and ped-sim first)
pixi run harness --only rm_base           # a subset; K and so alpha/K shrink; to .work/results/
pixi run harness --replicate 1            # the same cases on independent seeds; to .work/results/
pixi run power                            # rejection rates at stated errors (rm_base, ctm_base)
pixi run power --case std_mating          # any case whose oracle has mutations; to .work/results/power.md
pixi run lint
pixi run -e xftsim xft-reference          # layer 5's xftsim reference samples (skips current ones)
```

The script exits 1 if any row is FAIL or ERROR. There is no CI.

To check a local simACE checkout instead, without committing the change,
replace the `simace` entry in `pixi.toml` with
`simace = { path = "../..", editable = true }` and run `pixi install`.
`results/results.md` records the simace version and git ref of each run.

## What it checks

Each case simulates R pedigrees with `simace simulate` (seeds fixed in
[`cases.yaml`](cases.yaml)) and scores them in one of four layers. Layer 4
runs its seeds on GenomeDK and is scored from committed results.

**Layer 1: the liability law on simACE's own pedigree** (`oracles/liability.R`).
The cases cover random mating, with and without vertical C transmission
(simACE ADR 0025). For each pedigree, the oracle draws M = 199 realisations
of A, C and E on the same pedigree, conditional on recorded generation 0
(whose hidden parents simACE erases):

- A is the midparent plus Mendelian sampling, with MZ co-twins sharing it.
- C is b_m times the mother's C plus b_f times the father's, plus one draw
  per maternal household. That draw's covariance is the configured C minus
  the realised covariance of the transmitted part. Each generation's C
  variance is therefore the configured one in expectation, whatever the
  parents' was. With b_m = b_f = 0, C is the household draw alone.
- E is independent.

Every statistic is computed by one R function on simACE's values and on
each draw:

- correlations of A, C, E and liability for MZ, FS, MHS, PHS,
  mother–offspring (MO), father–offspring (FO) and 1C pairs
- Var(A), Var(C) and Var(liability) per generation
- the cross-trait correlations

The `ctm_analytic` case scores the C-transmission pedigrees against the
analytic residual covariance instead ([D2](differences.md)).

The randomised rank of simACE's value among the draws is uniform when
simACE samples from the stated law. Two tests score the ranks over the
seeds: an exact KS test, and a chi-square test on 10 equal bins referred to
its exact multinomial null (ADR 0001 §2). The asymptotic chi-square
reference rejects up to twice the nominal rate at the run's alpha/K (G4).
The deep-pedigree cases also run the infinitesimal-model oracle, whose
segregation variance shrinks with parental inbreeding ([D1](differences.md)).
The inbreeding coefficients come
from nadiv (`makeAinv`), installed by `pixi run build-nadiv`.

**Layer 2: gene-drop through simACE's pedigree** (`genedrop/`). Each
pedigree is dropped once through msprime's fixed-pedigree model, over 22
chromosomes of 1.5 M, with MZ co-twins collapsed to one genome. The case
checks two things:

- Realised IBD minus 2φ from kinship2 (`oracles/pairs.R`) averages 0 in
  every class. This tests that simACE's mother, father and twin links form
  the pedigree they claim.
- The Mendelian-sampling term of simACE's A matches that of a finite-locus
  genetic value from the dropped genomes, in variance per generation and
  sib correlations, in mean only. The variance rows score only children
  whose recorded parents both have F = 0, where the two laws agree;
  children of inbred parents are INFO rows `d1_innov_var_g*`
  ([D1](differences.md)).

The `pedsim_ibd` case checks the engine itself: msprime against ped-sim
IBD fractions on 2000 copies each of FS, HS and 1C families.

**Layer 3: pedigree structure** (`structure/`, `slim/wf.slim`). simACE's
`wright_fisher` pedigrees are compared with SLiM nonWF pedigrees under the
same law: mothers and fathers drawn uniformly with replacement, and
Bernoulli sex for founders and offspring. SLiM's WF model fixes the sex
ratio, so it is not used. One value per pedigree and generation is
compared with a two-sample KS test and a Welch t-test:

- sex ratio
- distinct mothers and fathers
- offspring-number variance per parent
- FS, MHS and PHS pair counts
- offspring of half- or full-sib matings

**Layer 3b: simACE's `standard` mating** (`structure/mating.py`,
`slim/standard.slim`). Two cases check the default mating model, both at
λ = 0.5 and p_mz = 0.02, an MZ rate per birth since simace v0.3.0
(rwaples/simACE#78).

`std_mating` tests simACE's own pedigrees against the documented rules,
with no second simulator. Each rule gets a statistic whose null conditions
on as much as that rule allows:

- sex ratio: Bin(N − MZ pairs, 1/2), exact
- MZ pairs: Bin(N, p_mz / (1 + p_mz)) capped at N / 2, exact
- children per couple and their split over a mother's couples, given the
  couples and the MZ pairs: the N − (MZ pairs) births by a multinomial
  conditioned on every couple having a birth, which does not depend on λ,
  then the twin births uniform among the births
- where the MZ pairs fall, given each couple's births (its children less
  its MZ pairs, all seen): couples with ≥ 1 MZ pair, and Σ b_k x_k over
  couples (births times MZ pairs), each against the twin births drawn
  uniformly among the births, multivariate hypergeometric
- couples, and mothers and fathers with ≥ 2 partners, given the parents'
  sexes and the MZ pairs: ZTP(λ) slots, trimmed to the smaller sex's
  total, then the matings the births occupy

A generation's children are the next one's parents, so each statistic
gets one randomised PIT per generation given the generations before it.
These are independent U(0, 1), and their Stouffer combination is one
exactly uniform PIT per pedigree, scored over the seeds like layer 1. The
derivations are in `structure/mating.py`'s docstring. They assume no
duplicate couple survives simACE's 5 swap passes:
`tools/residual_duplicates.py` found none in 200000 generations at
N = 2000 or at N = 500, a 95% bound of 1.5e-5 per generation.
`tools/check_mating_nulls.py` checks the samplers against exact
enumeration, every null against simACE's own mating, and the invariants.

Some rules are deterministic, not random:

- MZ co-twins share parents and sex, and have consecutive ids.
- Each couple's children have one block of consecutive ids, and mothers
  are in id order within a generation.

A couple can have several MZ pairs.

A violation raises an error, so the case becomes an ERROR row.

`std_structure` compares simACE's pedigrees with `slim/standard.slim`, the
layer-3 statistics plus couples, multi-partner parents, MZ pairs and the
sum of squared family sizes, one row per generation. The script is a
**second implementation of simACE's own rules, written from simACE's
documentation; it is not an outside model**. It can catch simACE's code
departing from its documented model. It cannot catch a flaw in the model
itself. A misreading of the docs shared by both implementations would
hide a bug.

Where the docs leave a point open, the script follows simACE's code and
says so. Those points went to simACE as
[rwaples/simACE#78](https://github.com/rwaples/simACE/issues/78):

- the larger sex loses random slots, not random persons;
- a duplicate couple has its father swapped with a uniformly random other
  mating for at most 5 passes, where the docs say "a nearby partner";
- children's sex is Bernoulli(1/2), stated in methods.md for founders
  only.

The same issue reported two id-order artifacts. simace v0.3.0 documents
the order (methods.md, "Offspring order"), and both cases now check it as
invariants above.

**Layer 4: simACE's downstream stages on a gene-drop A** (`simace_ts/`,
[ADR 0002](docs/adr/0002-simace-ts-layer.md)). Case `simace_ts_p2` drops
the 35,951 SimHumanity p2 genomes through simACE's 10-generation pedigree
with [`tskit_pipeline/`](tskit_pipeline/README.md) (real LD and allele
frequencies, the HapMapII GRCh38 map), assigns a finite causal architecture
with tstrait, rescales A1 to the configured variance in each generation, and
runs `simace cohort` and `simace analyze` on it. It checks that those
stages are right when A is a finite-locus value rather than a Gaussian
draw. The pipeline supplies A, so this layer depends on simACE more than
layers 1 to 3: its oracle consumes simACE's pedigree and feeds simACE's own
stages.

Each seed gives three kinds of row:

- **Scored:** a `report.yaml` field against its realised truth, computed
  from the gene-drop A simACE was handed and the configured C and E law:
  liability correlations of FS, MHS, PHS, MO and FO pairs, the full-sib
  correlation, and midparent slopes and correlations. A t-test over the
  R = 50 seeds says the mean gap is 0. `results.md` adds each row's minimum
  detectable gap.
- **Exact:** fields that depend only on the supplied values (pair counts,
  A correlations and slopes, the liability correlations themselves),
  recomputed with simACE's arithmetic. A mismatch is FAIL; these rows have
  no p-value and are not in K.
- **INFO:** quantities that assume Gaussian liability by simACE's design
  ([D3](differences.md)), the gene-drop A correlation minus 2φ by class and
  generation, the architecture sweep (100 to 1M causal sites, α = 0 or
  −0.5), and the layer-2 cross-check (G6, `g6/` rows; ADR 0002 item 11): the
  IBD of layer 2's drop and of the tskit layer's drop, each against the
  recorded-pedigree 2φ of a sample of pairs, and layer 2's A correlations
  minus the tskit layer's.

The case has no MZ twins (the pipeline drops co-twins as separate
offspring, so it refuses them), no ascertainment, no C transmission and no
assortative mating. A seed costs about 4 CPU-hours and 36 GB. Its data are
pinned in `cases.yaml`: the SimHumanity trees archive by sha256, the
SimHumanity checkout (recombination maps) by commit, and the preprocessed
founders by fingerprint.

```bash
# on GenomeDK, one job per seed (8 CPUs, 48 GB, 1 h); the ~10 GB of grafted .trees go to
# node-local disk, since home has a quota that 50 seeds at once exceed
pixi run simace-ts --seed 9000 --data ~/data/humanity_sim --work ~/simace_ts_work/scored \
  --trees-dir "${TMPDIR:-/tmp}/sts_$SLURM_JOB_ID"
# on the laptop: fetch each seed's part.json into simace_ts/parts/, then score
mkdir -p simace_ts/parts/simace_ts_p2
for d in $(ssh gdk 'ls -d ~/simace_ts_work/scored/simace_ts_p2/seed*'); do
  s=$(basename "$d"); rsync -a "gdk:$d/part.json" "simace_ts/parts/simace_ts_p2/$s.json"; done
# G6, after the seeds: one job per seed (4 CPUs, 40 GB, 6 h), written to the seed directory as g6.json
pixi run python simace_ts/g6.py ~/simace_ts_work/scored/simace_ts_p2/seed9000 --data ~/data/humanity_sim --threads 4
mkdir -p simace_ts/g6/simace_ts_p2
for d in $(ssh gdk 'ls -d ~/simace_ts_work/scored/simace_ts_p2/seed*'); do
  s=$(basename "$d"); rsync -a "gdk:$d/g6.json" "simace_ts/g6/simace_ts_p2/$s.json"; done
pixi run harness --only simace_ts_p2
pixi run test-pipeline                    # tests of tskit_pipeline/ and simace_ts/
```

A G6 file is refused when a seed has none, when it comes from another
`G6_VERSION` or from uncommitted code, when its `part_sha256` is not the
committed Part's or its tskit A correlations differ from the Part's by more
than 1e-12, or when seeds differ in pairs per class or genome. A rebuilt Part
therefore needs G6 rerun. A refusal is one `g6` ERROR row; the case's scored
rows are kept.

The harness never simulates this case. A seed without a Part, or a Part
whose case, simace, data pins or chromosomes differ from this run's, makes
the case ERROR, so a simace bump needs all R seeds rerun on GenomeDK. So does
a Part computed by another `PART_VERSION` of `simace_ts/part.py`, or from
uncommitted code; that needs only the Part recomputed from the seed's
outputs on GenomeDK, a few seconds each:
`pixi run python simace_ts/part.py ~/simace_ts_work/scored/simace_ts_p2/seed9000`.

**Layer 5: assortative mating against xftsim** (`am/`, [ADR
0003](docs/adr/0003-xftsim-am-layer.md)). Under assortative mating (AM),
simACE's pedigree is chosen on the liabilities, so no layer-1 law exists.
Instead, simACE's population moments are compared with those of
[xftsim](https://github.com/border-lab/xftsim), a finite-locus forward
simulator by Richard Border that shares no code with simACE. xftsim is
pinned to commit `8a4cdf0` in its own pixi environment (Python 3.9) and is
called, never vendored or copied (GPLv3); ADR 0003 gives the citation.

The two simulators run at a matched configuration:
- A 0.5, C 0, E 0.5 per trait (C 0.2, E 0.3 in `amc_r05`), with
  independent traits;
- N 10⁴ and 6 rounds of mating;
- xftsim with M = 1000 unlinked loci.

The cases are:
- `am_r0`: random mating, through xftsim's regime at r = 0;
- `am_r05` and `am_r08`: single-trait AM (`assort1`);
- `amc_r05`: single-trait AM at 0.5 with a household C of 0.2. xftsim
  builds C from the mother's private draw, passed to her children by a
  vertical component;
- `am2_r02` and `am2_r04`: two-trait AM, with every mate correlation r
  (`assort1 = assort2 = r` and `assort_matrix [[r, r], [r, r]]`). That is
  the law of xftsim's `LinearAssortativeMatingRegime` on two phenotypes.

Per generation, the cases compare:
- Var(A) and Var(L);
- the mate correlations on L and on A;
- the full-sib and parent–offspring correlations of A and L;
- for two traits, the within-person and cross-trait mate correlations.

**How it is scored:**
- xftsim's side is a committed **reference sample** (`am/reference/`), R
  runs per case made by `pixi run -e xftsim xft-reference`. Its manifest
  pins the xftsim lock, the scripts' hashes, the thread count and the
  seeds, and `pixi run harness` refuses a stale one.
- The harness reruns only simACE and scores each row by a Welch t over
  seeds. There is no KS test, because the two spreads differ by design.
- Each row carries an INFO companion `bulmer:<row>`, simACE against the
  infinitesimal Bulmer recursion.
- `am_finite_m` (INFO) compares xftsim at M = 250 with M = 4000, to bound
  the finite-locus gap ([D4](differences.md)).

The pilot behind the design is in `am/pilot/`.

**Multiplicity.** Every p-value of every scored row in a run forms one
Bonferroni family, α = 0.05/K. The `null` row of `results/power.md`
measures the false-positive rate of `rm_base` scored as its own family; a
full run's K is larger, so its threshold is stricter.

## Scope

[`simace validate`](https://github.com/rwaples/simACE) already tests
simACE's residuals against its configured parameters, with simACE's own
code. This repo adds independent code on the same pedigree. Some overlap is
deliberate: per-generation variance, and the variance and sib correlation
of the Mendelian residual.

Not covered yet ([`status.md`](status.md)):

- assortative mating with two traits and C > 0, with correlated traits, with two-trait
  cross terms other than √(r1·r2) (including simACE's derived one), or combined with C
  transmission; layer 5 covers single-trait AM and two-trait AM with a rank-1 mate
  correlation matrix (every correlation equal, or r1 0.4 and r2 0.1)
