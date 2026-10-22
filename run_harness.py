"""Run simACE-validate's cases and write results/results.md and results.tsv.

    pixi run harness                    # every case, written to results/
    pixi run harness --only rm_base     # a subset (K, and so alpha, shrinks)

A run with --only or --replicate is partial: it is written to .work/results/
and its results.md names the selection, so results/ always holds the full
run. Exits 1 if any row is FAIL or ERROR.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import subprocess
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
import polars as pl
import yaml

from calibration import Score, calibrate, mean_zero, randomised_pit, two_sample
from genedrop.drop import Genome
from genedrop.layer2 import seed_stats
from genedrop.pedsim import msprime_ibd, pedsim_ibd
from structure.mating import REPLACES_SIMACE, SLIM_MUTATIONS, mating_law, mating_stats
from structure.stats import run_slim, structure_stats

ROOT = Path(__file__).resolve().parent
WORK = ROOT / ".work"
RESULTS = ROOT / "results"
PARTIAL_RESULTS = WORK / "results"
PEDSIM = ROOT / ".pedsim" / "ped-sim"  # built by tools/build_pedsim.sh
ALPHA = 0.05
JOBS = 6
MIN_SEEDS = 20  # as simace validate's min_reps
SINGLE_THREAD = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1"}


def case_seed(case_id: str, salt: str) -> int:
    """Derive a stable 31-bit seed from the case id, so runs reproduce exactly."""
    return int.from_bytes(hashlib.sha256(f"{case_id}:{salt}".encode()).digest()[:4], "big") >> 1


def simulate(case: dict, seed: int) -> Path:
    """Run `simace simulate` for one seed and return the pedigree path."""
    out = WORK / case["id"] / f"seed{seed}.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    flags = [f"--{k}={v}" for k, v in case["simulate"].items()]
    cmd = ["simace", "simulate", f"--seed={seed}", *flags, f"--output-pedigree={out}"]
    subprocess.run(cmd, check=True, capture_output=True, env=os.environ | SINGLE_THREAD)
    return out


def run_oracle(case: dict, chunk: int, seeds: list[int], pedigrees: list[Path]) -> list[dict]:
    """Run the case's oracle over one chunk of pedigrees."""
    stem = WORK / case["id"] / f"oracle{chunk}"
    payload = {
        "params": {k: case["simulate"][k] for k in ("A1", "C1", "E1", "A2", "C2", "E2", "rA", "rC", "rE")},
        "M": case["M"],
        "law": case["law"],
        "mutation": case.get("mutation"),
        "oracle_seed": case_seed(case["id"], f"oracle{chunk}"),
        "seeds": seeds,
        "pedigrees": [str(p) for p in pedigrees],
    }
    Path(f"{stem}.in.json").write_text(json.dumps(payload))
    script = ROOT / "oracles" / f"{case['oracle']}.R"
    subprocess.run(
        ["Rscript", str(script), f"{stem}.in.json", f"{stem}.out.json"],
        check=True,
        capture_output=True,
        env=os.environ | SINGLE_THREAD,
    )
    return json.loads(Path(f"{stem}.out.json").read_text())


def z_score(row: dict) -> float:
    """Return simACE's value in oracle SDs; infinite when the oracle law is degenerate."""
    gap = row["simace"] - row["oracle_mean"]
    if row["oracle_sd"] > 0:
        return gap / row["oracle_sd"]
    return float(np.sign(gap)) * np.inf


def row(case: dict, stat: str, score: Score, info: bool = False) -> dict:
    """Return one results row; a statistic present in fewer than MIN_SEEDS seeds is INFO."""
    info = info or score.n < MIN_SEEDS
    return {"case": case["id"], "layer": case["layer"], "stat": stat, "score": score, "case_def": case, "info": info}


def run_liability(case: dict, pool: ThreadPoolExecutor, pedigrees: list[Path], seeds: list[int]) -> list[dict]:
    """Layer 1: run the liability oracle and calibrate every statistic over the seeds."""
    chunks = [(i, seeds[i::JOBS], pedigrees[i::JOBS]) for i in range(min(JOBS, len(seeds)))]
    per_seed = [r for out in pool.map(lambda c: run_oracle(case, *c), chunks) for r in out]
    rng = np.random.default_rng(case_seed(case["id"], "pit"))
    by_stat: dict[str, list[dict]] = {}
    for rep in sorted(per_seed, key=lambda r: r["seed"]):
        for r in rep["rows"]:
            by_stat.setdefault(r["stat"], []).append(r)
    rows = []
    for stat, reps in by_stat.items():
        n_less = np.array([r["n_less"] for r in reps])
        n_equal = np.array([r["n_equal"] for r in reps])
        z = np.array([z_score(r) for r in reps])
        rows.append(row(case, stat, calibrate(randomised_pit(n_less, n_equal, case["M"], rng), z)))
    return rows


def genedrop_seed(case: dict, seed: int, pedigree: Path) -> dict[str, float]:
    """Layer 2, one seed: classify pairs and their kinship in R, then gene-drop."""
    pairs_path = pedigree.with_suffix(".pairs.parquet")
    subprocess.run(
        ["Rscript", str(ROOT / "oracles" / "pairs.R"), str(pedigree), str(pairs_path)],
        check=True,
        capture_output=True,
        env=os.environ | SINGLE_THREAD,
    )
    genome = Genome(**case["genome"])
    A = case["simulate"]["A1"]
    return seed_stats(
        pl.read_parquet(pedigree), pl.read_parquet(pairs_path), genome, A, case_seed(case["id"], str(seed))
    )


def present(per_seed: list[dict], stat: str) -> list:
    """Return the stat's value in every seed that has it; a NaN value (too few pairs) counts as absent."""
    return [s[stat] for s in per_seed if stat in s and not np.any(np.isnan(s[stat]))]


def run_genedrop(case: dict, pool: ThreadPoolExecutor, pedigrees: list[Path], seeds: list[int]) -> list[dict]:
    """Layer 2: one gene-drop per pedigree; t-test each difference over the seeds that have it."""
    per_seed = list(pool.map(lambda sp: genedrop_seed(case, *sp), zip(seeds, pedigrees, strict=True)))
    rows = []
    for stat in dict.fromkeys(k for s in per_seed for k in s):
        values = np.array(present(per_seed, stat))
        if stat.startswith("spread_"):
            sim, gd = values.T
            ratio = float(np.std(sim, ddof=1) / np.std(gd, ddof=1))
            rows.append(row(case, stat, Score(n=len(sim), center=ratio, effect=np.nan, p_values={}), info=True))
        else:
            rows.append(row(case, stat, mean_zero(values)))
    return rows


def pedsim_check(case: dict, pool: ThreadPoolExecutor, pedigrees: list[Path], seeds: list[int]) -> list[dict]:
    """Layer 2 engine check: msprime against ped-sim IBD on fixed families (no simACE pedigrees)."""
    genome = Genome(**case["genome"], causal_per_chromosome=0)
    copies = case["copies"]
    work = WORK / case["id"]
    jobs = {
        "pedsim": pool.submit(pedsim_ibd, PEDSIM, work, genome, copies, case_seed(case["id"], "pedsim")),
        "msprime": pool.submit(msprime_ibd, genome, copies, case_seed(case["id"], "msprime")),
    }
    ped, ms = (jobs[k].result() for k in ("pedsim", "msprime"))
    return [row(case, f"pedsim_ibd_{name}", two_sample(ms[name], ped[name])) for name in ped]


def slim_pedigree(
    case: dict, script: str, label: str, seed: int, defines: dict[str, float] | None = None
) -> pl.DataFrame:
    """Run slim/<script> for one simACE seed at the case's N and G-ped; label keeps its SLiM seed apart."""
    sim = case["simulate"]
    out = WORK / case["id"] / f"{label}{seed}.tsv"
    return run_slim(
        ROOT / "slim" / script, out, sim["N"], sim["G-ped"], case_seed(case["id"], f"{label}{seed}"), defines
    )


def standard_defines(case: dict, mutation: str | None = None) -> dict[str, float]:
    """Return standard.slim's LAMBDA and PMZ: the case's, or a slim_standard mutation's."""
    sim = case["simulate"]
    base = {"LAMBDA": sim["mating-lambda"], "PMZ": sim["p-mztwin"]}
    return base | (SLIM_MUTATIONS[mutation][1] if mutation else {})


def run_structure(case: dict, pool: ThreadPoolExecutor, pedigrees: list[Path], seeds: list[int]) -> list[dict]:
    """Layer 3: structure statistics of simACE pedigrees against SLiM pedigrees, one per seed each.

    slim_wf compares wright_fisher pedigrees with slim/wf.slim. slim_standard
    compares standard pedigrees with slim/standard.slim and adds the mating
    statistics; its power mutations (structure/mating.py SLIM_MUTATIONS)
    change SLiM's LAMBDA or PMZ, and `null` replaces simACE with a second
    SLiM run.
    """
    mutation = case.get("mutation")
    standard = case["oracle"] == "slim_standard"
    if mutation and not standard:
        raise ValueError(f"{case['id']}: {case['oracle']} has no power mutations")

    def stats_of(ped: pl.DataFrame) -> dict[str, float]:
        return structure_stats(ped) | mating_stats(ped) if standard else structure_stats(ped)

    def slim_one(seed: int) -> dict[str, float]:
        if standard:
            return stats_of(slim_pedigree(case, "standard.slim", "slim", seed, standard_defines(case, mutation)))
        return stats_of(slim_pedigree(case, "wf.slim", "slim", seed))

    if mutation == "null":
        ours = list(
            pool.map(lambda s: stats_of(slim_pedigree(case, "standard.slim", "null", s, standard_defines(case))), seeds)
        )
    else:
        ours = list(pool.map(lambda p: stats_of(pl.read_parquet(p)), pedigrees))
    theirs = list(pool.map(slim_one, seeds))
    return [
        row(case, stat, two_sample(np.array([o[stat] for o in ours]), np.array([t[stat] for t in theirs])))
        for stat in ours[0]
    ]


def run_mating_law(case: dict, pool: ThreadPoolExecutor, pedigrees: list[Path], seeds: list[int]) -> list[dict]:
    """Layer 3b: score each standard pedigree against the mating rules' own nulls (structure/mating.py)."""
    sim = case["simulate"]
    mutation = case.get("mutation")

    def one(i: int) -> dict[str, tuple[float, float]]:
        if mutation == "null_slim":  # power only: a generator of the documented law in place of simACE
            ped = slim_pedigree(case, "standard.slim", "slim", seeds[i], standard_defines(case))
        else:
            ped = pl.read_parquet(pedigrees[i])
        rng = np.random.default_rng(case_seed(case["id"], f"law{seeds[i]}"))
        return mating_law(ped, sim["mating-lambda"], sim["p-mztwin"], case["M"], rng, mutation)

    per_seed = list(pool.map(one, range(len(seeds))))
    rows = []
    for stat in dict.fromkeys(k for s in per_seed for k in s):
        u, z = np.array([s[stat] for s in per_seed if stat in s]).T
        rows.append(row(case, stat, calibrate(u, z)))
    return rows


RUNNERS = {
    "liability": run_liability,
    "genedrop": run_genedrop,
    "pedsim": pedsim_check,
    "slim_wf": run_structure,
    "slim_standard": run_structure,
    "mating_law": run_mating_law,
}


def run_case(case: dict) -> list[dict]:
    """Simulate the case's pedigrees and score them with its oracle."""
    seed0 = case.get("seed0", 0)
    seeds = list(range(seed0, seed0 + case.get("R", 0)))
    replaced = case.get("mutation") in REPLACES_SIMACE.get(case["oracle"], ())
    with ThreadPoolExecutor(JOBS) as pool:
        pedigrees = [] if replaced else list(pool.map(lambda s: simulate(case, s), seeds))
        return RUNNERS[case["oracle"]](case, pool, pedigrees, seeds)


def assign_status(rows: list[dict]) -> float:
    """Set each row's status against one Bonferroni family over every p-value; return alpha/K."""

    def scored(r: dict) -> bool:
        return not r["info"] and r["case_def"].get("expect", "pass") != "info"

    K = sum(len(r["score"].p_values) for r in rows if scored(r))
    threshold = ALPHA / max(K, 1)
    for r in rows:
        if not scored(r):
            r["status"] = "INFO"
        else:
            r["status"] = "FAIL" if r["score"].p < threshold else "PASS"
    return threshold


def simace_ref() -> str:
    """Return the simace version and the git ref it was installed from."""
    dist = importlib.metadata.distribution("simace")
    ref = ""
    direct = dist.read_text("direct_url.json")
    if direct:
        info = json.loads(direct)
        vcs = info.get("vcs_info", {})
        ref = vcs.get("requested_revision") or vcs.get("commit_id", "")[:12] or info.get("url", "")
    return f"simace {dist.version}" + (f" ({ref})" if ref else "")


def write_results(rows: list[dict], threshold: float, errors: list[tuple[dict, str]], selection: str) -> None:
    """Write results.md and results.tsv: a full run to results/, a partial one (``selection``) to .work/results/."""
    out = PARTIAL_RESULTS if selection else RESULTS
    out.mkdir(parents=True, exist_ok=True)
    totals: dict[str, int] = {}
    for r in rows:
        totals[r["status"]] = totals.get(r["status"], 0) + 1
    if errors:
        totals["ERROR"] = len(errors)
    total_line = ", ".join(f"{k} {v}" for k, v in sorted(totals.items()))
    header = ["case", "layer", "stat", "n", "center", "effect", "p_values", "alpha_K", "status", "id"]
    lines = ["\t".join(header)]
    md = [
        "# Results",
        "",
        "Generated by `run_harness.py`; do not edit. A row is PASS when every p-value is at least alpha/K,",
        "with K the number of p-values over every scored row of the run (ADR 0001). Layer 1: center is the",
        "mean PIT (0.5 under the null) and effect the mean z of simACE against the oracle draws. Layer 2:",
        "center is the mean per-seed difference (0 under the null) and effect its t statistic; a spread row",
        "is the SD ratio simACE / gene-drop over seeds.",
        "",
        f"**Totals:** {total_line}",
        "",
        f"{simace_ref()}; alpha/K = {threshold:.3g}.",
        "",
        *([f"**Partial run:** {selection}. K covers these cases only.", ""] if selection else []),
        "| Case | Layer | Statistic | n | Center | Effect | p | Status | Id |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        c, case = r["score"], r["case_def"]
        did = case.get("id_diff", "")
        pv = ", ".join(f"{k} {v:.3g}" for k, v in c.p_values.items())
        vals = [r["case"], r["layer"], r["stat"], c.n, f"{c.center:.4g}", f"{c.effect:+.2f}", pv]
        lines.append("\t".join(map(str, [*vals, f"{threshold:.3g}", r["status"], did])))
        md.append("| " + " | ".join(map(str, [*vals, r["status"], did])) + " |")
    for case, message in errors:
        lines.append(
            "\t".join(map(str, [case["id"], case["layer"], message, "", "", "", "", f"{threshold:.3g}", "ERROR", ""]))
        )
        md.append(f"| {case['id']} | {case['layer']} | {message} | | | | | ERROR | |")
    (out / "results.tsv").write_text("\n".join(lines) + "\n")
    (out / "results.md").write_text("\n".join(md) + "\n")
    print(f"Totals: {total_line}; written to {out.relative_to(ROOT)}/")


def main() -> int:
    """Run the selected cases."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", help="comma-separated case ids")
    parser.add_argument(
        "--replicate", type=int, default=0, help="shift every seed (simACE and oracle) to an independent replicate run"
    )
    args = parser.parse_args()
    cases = yaml.safe_load((ROOT / "cases.yaml").read_text())["cases"]
    if args.only:
        wanted = set(args.only.split(","))
        unknown = wanted - {c["id"] for c in cases}
        if unknown:
            parser.error(f"unknown case ids: {', '.join(sorted(unknown))}")
        cases = [c for c in cases if c["id"] in wanted]
    if args.replicate:
        # The new id also moves every seed derived by case_seed (oracle, gene-drop, ped-sim, SLiM).
        shift = 100_000 * args.replicate
        cases = [
            c | {"id": f"{c['id']}~{args.replicate}"} | ({"seed0": c["seed0"] + shift} if "seed0" in c else {})
            for c in cases
        ]
    parts = [f"--only {args.only}"] if args.only else []
    if args.replicate:
        parts.append(f"--replicate {args.replicate}")
    selection = "; ".join(parts)

    rows: list[dict] = []
    errors: list[tuple[dict, str]] = []
    for case in cases:
        print(f"{case['id']}: layer {case['layer']}", flush=True)
        try:
            rows += run_case(case)
        except subprocess.CalledProcessError as exc:
            errors.append((case, f"{' '.join(map(str, exc.cmd[:2]))} exited {exc.returncode}"))
            print(exc.stderr.decode(errors="replace")[-2000:], file=sys.stderr)
        except Exception as exc:  # one broken case is an ERROR row, not a lost run
            errors.append((case, f"{type(exc).__name__}: {exc}"))
            traceback.print_exc()
    threshold = assign_status(rows)
    write_results(rows, threshold, errors, selection)
    return 1 if errors or any(r["status"] == "FAIL" for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
