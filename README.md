# simACE against independent simulators

Runs [simACE](https://github.com/rwaples/simACE) and checks what it
simulates against code simACE does not share:

- an R liability oracle
- [msprime](https://tskit.dev/msprime/) gene-drop
- [ped-sim](https://github.com/williamslab/ped-sim)
- [SLiM](https://messerlab.org/slim/)

The latest results are in [`results/results.md`](results/results.md), and
the harness's measured power is in [`results/power.md`](results/power.md).
Every known model gap is explained in [`differences.md`](differences.md),
and where each non-PASS row stands is in [`status.md`](status.md).

simACE's output is random, so a check is a calibrated statistical test,
not a tolerance ([ADR 0001](docs/adr/0001-conditional-calibration.md)).
That is the main difference from
[pedsum-validate](https://github.com/rwaples/pedsum-validate) and
[pg-phenotype-validate](https://github.com/rwaples/pg-phenotype-validate).

## Run it

With [pixi](https://pixi.sh). `pixi.toml` pins simace to the `v0.2.0`
family tag.

```bash
pixi run harness                          # every case, to results/ (builds nadiv and ped-sim first)
pixi run harness --only rm_base           # a subset; K and so alpha/K shrink; to .work/results/
pixi run harness --replicate 1            # the same cases on independent seeds; to .work/results/
pixi run power                            # rejection rates at stated errors
pixi run lint
```

The script exits 1 if any row is FAIL or ERROR. There is no CI.

To check a local simACE checkout instead, without committing the change,
replace the `simace` entry in `pixi.toml` with
`simace = { path = "../..", editable = true }` and run `pixi install`.
`results/results.md` records the simace version and git ref of each run.

## What it checks

Each case simulates R pedigrees with `simace simulate` (seeds fixed in
[`cases.yaml`](cases.yaml)) and scores them in one of three layers.

**Layer 1: the liability law on simACE's own pedigree** (`oracles/liability.R`).
The case covers random mating with no C transmission. For each pedigree,
the oracle draws M = 199 realisations of A, C and E on the same pedigree,
conditional on recorded generation 0 (whose hidden parents simACE erases):

- A is the midparent plus Mendelian sampling, with MZ co-twins sharing it.
- C is one draw per maternal household.
- E is independent.

Every statistic is computed by one R function on simACE's values and on
each draw:

- correlations of A, C, E and liability for MZ, FS, MHS, PHS, PO and 1C
  pairs
- Var(A) and Var(liability) per generation
- the cross-trait correlations

The randomised rank of simACE's value among the draws is uniform when
simACE samples from the stated law. A KS test and a chi-square test score
the ranks over the seeds. The deep-pedigree cases also run the
infinitesimal-model oracle, whose segregation variance shrinks with
parental inbreeding ([D1](differences.md)). The inbreeding coefficients
come from nadiv (`makeAinv`), installed by `pixi run build-nadiv`.

**Layer 2: gene-drop through simACE's pedigree** (`genedrop/`). Each
pedigree is dropped once through msprime's fixed-pedigree model, over 22
chromosomes of 1.5 M, with MZ co-twins collapsed to one genome. The case
checks two things:

- Realised IBD minus 2φ from kinship2 (`oracles/pairs.R`) averages 0 in
  every class. This tests that simACE's mother, father and twin links form
  the pedigree they claim.
- The Mendelian-sampling term of simACE's A matches that of a finite-locus
  genetic value from the dropped genomes, in variance per generation and
  sib correlations, in mean only.

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

- simACE's `standard` mating model in layer 3, which needs a
  re-implementation, not an outside model
- C transmission (simACE ADR 0025)
- assortative mating, whose liability law is not 2φA
