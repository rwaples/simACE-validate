# simACE-validate

Checks what simACE simulates against code simACE does not share: an R
liability oracle, msprime gene-drop, ped-sim and SLiM. Terms are in
GLOSSARY.md, decisions in docs/adr/.

- Everything runs through pixi: `pixi run harness`, `pixi run power`,
  `pixi run lint`. `pixi.lock` changes only on purpose.
- Checks are calibrations, not tolerances (ADR 0001). Never raise α,
  shrink K or drop a statistic to make a row pass. A known model gap is an
  INFO case with a D entry in differences.md.
- Never edit an oracle or the harness while a `harness` or `power` run is
  live. Each Rscript call rereads its script.
- Never copy GPL code (ped-sim, nadiv) into the repo. tools/build_pedsim.sh
  builds the pinned ped-sim tag into .pedsim/; tools/build_nadiv.sh installs
  the pinned nadiv version into .rlib/.
- Don't commit the local-checkout override of `simace` in pixi.toml.
- Every non-PASS row has an entry in status.md. Update it when a run changes
  a row's result.
- Before reporting, quote:
  - `pixi run lint`
  - `pixi run harness` with its exit code
  - the totals line of results/results.md
  - for a change to R, M or an oracle law, `pixi run power` and its table
