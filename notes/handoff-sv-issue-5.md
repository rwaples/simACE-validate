# Handoff: sv-issue-5 (simACE's vertical C transmission against xftsim)

rwaples/simACE-validate#5, branch `sv-issue-5` in
`/data/Documents/simACE/external/simACE-validate/.claude/worktrees/sv-issue-5`,
from `7ba2e5e`. Design B (layer 5 extended in place, our own half-sib regime,
every reference regenerated on GenomeDK) was chosen by the user on 2026-11-09.

## Open items

- **Not merged, not pushed.** Merge `sv-issue-5` into main and `/push` are the user's.
- **#5 is still open.** Close it only on the user's word. Every "Done when" item is met (below).
- **GenomeDK scratch.** `~/sv-issue-5` on gdk (628M) holds the archive of 544f816, a
  private pixi 0.84.1 (`~/sv-issue-5/pixi/bin/pixi`, copied from the laptop; the
  cluster's own `~/.pixi/bin/pixi` is still 0.82.0) and the job logs. Remove it when no
  longer needed: `ssh gdk 'rm -rf ~/sv-issue-5'`.
- **The cluster's simACE-validate checkout** (`~/Documents/simACE/external/simACE-validate`)
  is at 49a923c and pins `requires-pixi >=0.82,<0.83`. The repo now requires 0.84.1, so
  a cluster run from that checkout needs a pull and a pixi update (runbook
  `/data/Documents/simACE/notes/genomedk.md`, which still says 0.82).
- **Not done:** `null_xft` power for the new cases; a second (b_m, b_f, r) under AM; two-trait C
  transmission. None was asked for.

## What was done

Commits on `sv-issue-5` (`git log -1 --format='%h %s'` for each):

- `544f816 Add C-transmission cases to layer 5 (simACE-validate#5)`
- `13ce008 Regenerate layer-5 references and score C transmission`
- `8aa6bba Record C transmission against xftsim (simACE-validate#5)`
- this note, committed after them

### The issue's open questions (probes in `am/ctm/`)

1. **Graph.** `python -m am.ctm.probe_graph {dual|issue} 1 20000 200 4 0.5 0.3 0.2`
   (xftsim env, random mating). Both constructions run; the issue's needs the proband c input
   labelled `comp_type = "outcome"` (`"output"` duplicates the column and founder
   initialisation raises `InvalidIndexError: Reindexing only valid with uniquely valued Index
   objects`). Dual output, expected Var(c) 0.2, r_MO_c 0.5, r_FO_c 0.3, r_FS_y 0.45,
   r_MO_y 0.35, r_FO_y 0.31:

   ```
   {"generation": 0, "var_c": 0.1983, "max_abs_y_minus_A_E_c_g1": 0.0}
   {"generation": 1, "var_c": 0.2042, "r_FS_c": 1.0, "r_MO_c": 0.5086, "r_FO_c": 0.2919, "r_FS_y": 0.4592, "r_MO_y": 0.3467, "r_FO_y": 0.3152, "max_abs_input_vs_mother_c": 0.0}
   {"generation": 4, "var_c": 0.1982, "r_FS_c": 1.0, "r_MO_c": 0.4909, "r_FO_c": 0.3027, "r_FS_y": 0.4497, "r_MO_y": 0.3494, "r_FO_y": 0.3244, "max_abs_input_vs_mother_c": 0.0}
   ```

2. **Households with several partners.** simACE v0.3.0 assigns households by mother
   (`simulate.py:1136`, `:1184`), the transmitted part uses each child's own father
   (`:1332`). `python -m am.ctm.probe_halfsib SEED 20000 4 0.5 0.3 0.2 DIR`, seeds 11 to 18:

   ```
   MHS n_seeds 8 mean 0.9075 se 0.0019 expected 0.91 z -1.36 share_hh [1.0]
   PHS n_seeds 8 mean 0.0746 se 0.0127 expected 0.09 z -1.22 share_hh [0.0]
   MO  n_seeds 8 mean 0.4979 se 0.0033 expected 0.5  z -0.65
   FO  n_seeds 8 mean 0.2963 se 0.0031 expected 0.3  z -1.21
   ```

   (Seed 1 alone gave PHS 0.0474; the 8-seed mean shows that was chance.)
   xftsim's regimes make no half-sibs, so `RandomHalfSibMatingRegime` (`am/run_xft.py`)
   was written on the `MatingRegime` interface without reading any xftsim regime. A probe
   showed `enumerate_assignment` keeps a female listed with two distinct males as two couples.

3. **Founders.** simACE draws founder C iid N(0, C) (`simulate.py:1746`); xftsim's
   `founder_variances` (C, C, 1) gives founder Var(c) = C (0.1983 at g0 above).

### Pilot (40 seeds per side, 9300001 to 9300040; `am/pilot/ctm_summary.md`)

Random mating: 93 rows, smallest Welch p 0.00341 (`r_PHS_C_g4`), 6 below 0.05 (4.7 expected).
Under AM r 0.5, Var(C) − C at g 6: simACE −0.0008 (CI −0.0026, +0.0011), xftsim +0.0122
(+0.0109, +0.0135).

### References (GenomeDK jobs 7742108 to 7742111, each COMPLETED 0:0, 4:42 to 6:24)

`tools/compare_reference_records.py` (committed refs against the cluster's): all 18 existing
references hold the same draws, max |diff| 6.66e-16 to 1.33e-15, none byte-identical (CPU
rounding); "0 reference(s) differ beyond 1e-12 or are missing". On the laptop, seed 0 of
am_r05 and amc_r05 reproduced the old records exactly.

## Checks run

- `pixi run lint`: "All checks passed!", "84 files already formatted".
- `pixi run test-pipeline`: "120 passed, 3 warnings", exit 0.
- `pixi run harness` (full, exit 0): "Totals: INFO 1185, PASS 1065" against the
  baseline INFO 914, PASS 972 (K 1368). K is 1461, alpha/K 3.42e-05. Of 1886 shared rows, 0
  changed status, 0 were removed; 364 were added: ctm_rm 93 PASS and 186 INFO, ctm_am 85
  INFO. ctm_rm's smallest Welch p is 0.0145 (`r_FS_A_g3`). ctm_am `law:var_C_g*`: simACE
  −0.000758 to +0.000609 (|t| ≤ 1.63); xftsim +0.006043 to +0.01155 (t +14.06 to +26.53).
- `power.py --case ctm_rm --mutation null,ctm_swap,ctm_half --reps 10` (results/power_ctm.md):
  null 0/10, ctm_swap 10/10 (40 rows per run), ctm_half 10/10 (47.8 rows per run).
