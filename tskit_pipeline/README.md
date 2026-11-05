# Gene drop and tstrait phenotyping

The standard simACE pipeline draws additive genetic values $A$ from a Gaussian under an infinitesimal model. The gene-drop pipeline replaces that parametric $A$ with a genetic value computed from genotypes. It drops an ancestral tree sequence from SimHumanity through a simACE pedigree, sums causal-site effects per individual, and rescales the result to the configured $A$ variance. The phenotype, sampling, and plotting stages then run unchanged.

This page explains the design. The worked sweep over `alpha` and `num_causal` is in [Architecture sweep](architecture-sweep.md). The parameters are in [Parameters](#parameters).

This pipeline moved here from simACE's `scripts/gene_drop/`. Unlike the oracles in `genedrop/`, `oracles/` and `structure/`, it is not independent of simACE: it consumes a simACE `pedigree.parquet` and feeds the result back to `simace cohort` and `simace analyze`. Layer 4 (`simace_ts/`, ADR 0002) runs it per seed on GenomeDK, and `pixi run harness` scores the committed results; the harness itself never runs the pipeline. simACE's config no longer has the `tstrait`, `tskit_preprocess`, `drop_from` or `use_gene_drop` keys, so `scenarios.yaml` records the parameters of past runs and is not a simACE config file. Each script takes its parameters as flags.

## Why gene drop

The Gaussian $A$ model treats trait genetics as a sum of infinitely many tiny effects. The resulting distribution is smooth and says nothing about how heritability arises from a finite set of causal variants.

Gene drop simulates that mechanism:

1. Founder haplotypes carry SLiM-derived mutations from the SimHumanity ancestral simulation.
2. Mutations pass through the simACE pedigree by Mendelian inheritance, with recombination under the HapMapII GRCh38 rate map.
3. A finite subset of those mutations is declared causal. Each effect size is drawn from $\mathcal{N}(0, \sigma_\beta^2)$ and multiplied by $[2p(1-p)]^{\alpha}$. With $\alpha=-0.5$ this is the LDAK-thin form of Speed et al.
4. Each individual's genetic value is the sum of dosage times effect across causal sites.

The shape and tails of the realised $A$ distribution then follow from the genetic architecture. That is what you want when stress-testing methods that estimate variance components from pedigree data.

After computing $A$, the pipeline rescales it to match the configured `A1` variance. The phenotype models then see the same scale as a Gaussian $A$ run, and every downstream plot compares the two directly.

## Pipeline

The gene-drop branch sits between `simace simulate` and the phenotype stage. It has one preprocessing pass, run once per SimHumanity dataset, and five per-replicate steps. Each box is a script in `scripts/`:

```
Preprocessing (once per SimHumanity dataset)
────────────────────────────────────────────
SimHumanity p2 .trees ──► canonicalize_chrom.py (×22)
                                        │
                                        ▼
                          concat_chroms.py ──► verify.py
                                        │
                                        ▼
                          tstrait_site_catalog_chrom.py (×22)
                                        │
                                        ▼
                          tstrait_site_catalog_concat.py ──► site_catalog.parquet


Per replicate
─────────────
simace simulate ──► pedigree.parquet
              │
              ▼
    build_pedigree_tables.py ──► pedigree_tables.trees
              │
              ▼
    genotype_drop_chrom.py (×22) ──► genotypes_chrom_{n}.trees
              │
              ▼
    tstrait_assign_effects.py ──► causal_effects.parquet
              │
              ▼
    tstrait_gv_chrom.py (×22) ──► gv_chrom_{n}.parquet
              │
              ├──► tstrait_augment_pedigree.py ──► pedigree.full.tstrait.parquet
              │              │
              │              ▼
              │    simace cohort → simace analyze (Step 5)
              │
              └──► tstrait_phenotype.py ──► tstrait_phenotype.parquet (optional)
```

`simace run` never ran gene drop ([ADR 0020](https://github.com/rwaples/simACE/blob/master/docs/adr/0020-standalone-simace-cli.md)). Run each script by hand in this repo's pixi env, for example `pixi run python tskit_pipeline/scripts/canonicalize_chrom.py --help`. Each script's docstring shows the command its step needs, including its seed offset.

### Step 1: preprocess the SimHumanity ancestry

`canonicalize_chrom.py` canonicalises each per-chromosome `.trees` file from SimHumanity: it keeps one mutation per site and sorts samples into a fixed order. `concat_chroms.py` concatenates the chromosomes and writes a fingerprint sidecar. `verify.py` checks the result against that sidecar and against a golden summary; the first run for a dataset writes the golden file (`SIMACE_TSKIT_WRITE_GOLDEN=1`). `tstrait_site_catalog_chrom.py` and `tstrait_site_catalog_concat.py` then build `site_catalog.parquet`, which records the founder allele frequency of every site. Every gene-drop scenario reads that one catalog.

### Step 2: drop founders through the simACE pedigree

`build_pedigree_tables.py` builds the msprime fixed pedigree from the simACE pedigree once per replicate, tagging the latest `G_pheno` generations as samples. It refuses a pedigree with MZ twins, since co-twins would be dropped as separate offspring, and one whose ids are not its row order. `genotype_drop_chrom.py` loads it for each chromosome and runs `msprime.sim_ancestry(model="fixed_pedigree")` with the per-chromosome SimHumanity recombination map, read from the SimHumanity checkout named by `--simhumanity-dir`. It then grafts the dropped tree onto the preprocessed ancestral tree at the founder generation. The result is a chromosome-scale tree sequence whose present-day samples carry alleles inherited from the SimHumanity founder haplotypes.

This is the slowest step. To reuse one set of drops across several architectures, see [Sharing drops across variants](#sharing-drops-across-variants).

### Step 3: assign causal effects and compute genetic values

`tstrait_assign_effects.py` applies the MAF filter to the site catalog and samples `num_causal` sites from what remains. If `frac_causal` is set instead, it samples that fraction of the eligible sites. For each site it draws a raw effect $\beta \sim \mathcal{N}(\mu, \sigma_\beta^2)$ and multiplies by $[2p(1-p)]^{\alpha}$. With a per-replicate seed each replicate gets its own causal sites; with the scenario-level seed (`share_architecture`) one effects file serves every replicate.

`tstrait_gv_chrom.py` computes each individual's genetic value with a single-pass kernel compiled by numba on `tskit.jit.numba.NumbaTreeSequence`. The kernel walks each chromosome's trees once, keeps `left_child` and `right_sib` arrays up to date as it goes, and for each causal site does a depth-first walk from the mutation node to add the effect to every descendant sample. Its output matches `tstrait.genetic_value` up to floating-point summation order.

Only sample-tagged individuals, the generations covered by `G_pheno`, get a genetic value. Older generations have no sample nodes in the dropped tree.

### Step 4: rescale and overwrite A in the pedigree

`tstrait_augment_pedigree.py` reads `pedigree.parquet` and sums the per-chromosome genetic values into one genome-wide value per sample individual. It then:

1. Centres the values at zero.
2. Rescales them so that $\mathrm{Var}(A_\text{new})$ equals the configured `A1`. The variance is the sample variance with `ddof=0`. With `--per-generation`, steps 1 and 2 are done within each generation, so every generation has variance `A1` (ADR 0002).
3. Overwrites the `A1` column for sample individuals.
4. Recomputes `liability1` as $A_\text{new} + C + E$.
5. Writes `pedigree.full.tstrait.parquet` next to the original. The original is unchanged.

Older ancestors keep their parametric `A1`. That asymmetry is deliberate. Only sample-generation individuals reach the downstream model fits that read the phenotype.

`A2` is untouched because gene drop is single-trait. Sharing the gene-drop $A$ between both traits would need an extension to the augment script.

### Step 5: feed the standard pipeline

A gene-drop replicate then runs the standard stages on `pedigree.full.tstrait.parquet` instead of `pedigree.parquet`. Run the stage subcommands by hand: `simace cohort --pedigree pedigree.full.tstrait.parquet` writes `cohort.parquet` and `phenotyped_population.yaml`, and `simace analyze` takes the same pedigree with those two files. The phenotype model, censoring, and ascertainment then behave as they would for a Gaussian $A$ scenario, on the gene-drop $A$ column.

`simace cohort` and `simace analyze` read only pedigrees that carry the results-layout marker `simace_layout=2` ([ADR 0021](https://github.com/rwaples/simACE/blob/master/docs/adr/0021-two-canonical-replicate-parquets.md)). `tstrait_augment_pedigree.py` writes through simace's `write_pedigree`, so its output carries the marker and its input's C-transmission clip records. It refuses an input `pedigree.parquet` without the marker, from an older replicate.

## Sharing drops across variants

The drop and graft step is expensive. To vary the architecture while holding the genotypes fixed, point the variant at the base scenario's outputs (recorded as `drop_from: <base_scenario>` in `scenarios.yaml`). The architecture is the `tstrait_assign_effects.py` flags (`--alpha`, `--num-causal`, ...). Run `tstrait_assign_effects.py` per variant, then `tstrait_gv_chrom.py` and `tstrait_augment_pedigree.py` on the base scenario's `.trees` and `pedigree.parquet`. Only the tstrait steps run per variant.

A six-scenario sweep over `num_causal` in {100, 1k, 10k, 100k, 1m} and `alpha` in {0, -0.5} therefore runs the drop once instead of six times.

## Heritability

No step chooses $h^2$ on its own. `tstrait_augment_pedigree.py --target-var-a` takes `A1`, and `tstrait_phenotype.py --h2` takes the value from the standard simACE variance components:

$$h^2 = \frac{A_1}{A_1 + C_1 + E_1}$$

Gene-drop scenarios therefore sit on the same heritability scale as Gaussian $A$ scenarios. The `h2_realized` field in `tstrait_phenotype_meta.json` reports how close each replicate came to that value.

## Comparison

| | Gaussian $A$ (standard simACE) | Gene drop |
|---|---|---|
| $A$ comes from | $\mathcal{N}(0, A_1)$ per individual | sum of dosage times effect over causal sites, rescaled to $A_1$ |
| Architecture | infinitesimal | finite causal variants under an $\alpha$-MAF model |
| Cost per replicate | seconds | minutes, dominated by the drop and graft step |
| Scale | same by construction | same after the rescale step |
| Multi-trait | yes, trait1 and trait2 | trait1 only. trait2 stays parametric |
| Parameters | `A1`, `C1`, `E1`, `assort1`, and the other pedigree keys | the same, plus the [tstrait parameters](#parameters) |

## Parameters

`tstrait_assign_effects.py` takes the architecture as flags. Except for
`--num-causal` and `--frac-causal`, exactly one of which is set, all are
required. The values below are the defaults simACE's config used to supply.

| Flag | Type | Old default | Description |
|---|---|---:|---|
| `--num-causal` | int | `1000` | Number of causal sites. Set exactly one of this flag and `--frac-causal` |
| `--frac-causal` | float | unset | Fraction of sites that pass the MAF filter to use as causal |
| `--maf-threshold` | float | `0.01` | Minimum minor-allele frequency. `0` disables the filter |
| `--alpha` | float | `-0.5` | Exponent of the effect-size dependence on allele frequency |
| `--effect-mean` | float | `0.0` | Mean of the raw effect sizes before frequency scaling |
| `--effect-var` | float | `1.0` | Variance of the raw effect sizes before frequency scaling |
| `--trait-id` | int | `0` | Trait that gets the genetic value. Only trait 1, ID `0`, is supported |

Sharing one architecture across replicates is not a flag of
`tstrait_assign_effects.py`: run it once with the scenario-level seed and
pass its output to every replicate. `tstrait_phenotype.py --share-architecture`
records that choice in its metadata.

The preprocessing scripts take these flags. Nothing defaults the population:
without `--pop-name`, `canonicalize_chrom.py` keeps the samples of every
population as founders.

| Script | Flag | Description |
|---|---|---|
| `canonicalize_chrom.py` | `--src` | One SimHumanity per-chromosome `.trees` file (autosomes only) |
| `canonicalize_chrom.py`, `concat_chroms.py` | `--pop-name` | Founder population to keep; the sweep used `p2` |
| `concat_chroms.py` | `--source-dir` | Directory of the SimHumanity `.trees` files, recorded in the summary (required) |
| `concat_chroms.py` | `--chroms` | Autosomes to include; default `1` through `22` |
| `genotype_drop_chrom.py` | `--simhumanity-dir` | SimHumanity checkout holding the recombination maps (required) |
| `verify.py` | `--golden` | Golden summary to check against, or to write on the first run (required) |

The maintainer's SimHumanity trees for the sweep were in
`/data/Documents/humanity_sim/simhumanity_trees_RO`, preprocessed into
`/data/Documents/humanity_sim/preprocessed_p2`.
