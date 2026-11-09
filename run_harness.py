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
import math
import os
import subprocess
import sys
import traceback
from concurrent.futures import ThreadPoolExecutor
from functools import partial
from pathlib import Path

import numpy as np
import polars as pl
import yaml
from scipy import stats

from am.bulmer import bulmer2, predict1, predict_ctm
from am.manifest import case_spec, derived, manifest, read_reference, schema
from am.simace_stats import pedigree_rows
from am.stats import flatten
from calibration import Score, calibrate, mean_zero, randomised_pit, two_sample, welch
from genedrop.drop import Genome
from genedrop.layer2 import seed_stats
from genedrop.pedsim import msprime_ibd, pedsim_ibd
from simace_ts.score import case_rows, load_g6, load_parts
from structure.mating import REPLACES_SIMACE, SLIM_MUTATIONS, mating_law, mating_stats
from structure.stats import run_slim, structure_stats
from subproc import run

ROOT = Path(__file__).resolve().parent
WORK = ROOT / ".work"
RESULTS = ROOT / "results"
PARTIAL_RESULTS = WORK / "results"
PEDSIM = ROOT / ".pedsim" / "ped-sim"  # built by tools/build_pedsim.sh
ALPHA = 0.05
JOBS = 6
MIN_SEEDS = 20  # as simace validate's min_reps
SINGLE_THREAD = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "NUMBA_NUM_THREADS": "1"}
# The `simulate` flags the liability oracle's law accounts for. Any other flag
# (assortative mating, per-generation C or E) would change simACE's law unseen.
LIABILITY_FLAGS = {
    *("N", "G-ped", "G-sim", "mating-model", "mating-lambda", "p-mztwin"),
    *("A1", "C1", "E1", "A2", "C2", "E2", "rA", "rC", "rE"),
    *(f"Ctm-{side}{t}" for side in ("maternal", "paternal") for t in (1, 2)),
}


def case_seed(case_id: str, salt: str) -> int:
    """Derive a stable 31-bit seed from the case id, so runs reproduce exactly."""
    return int.from_bytes(hashlib.sha256(f"{case_id}:{salt}".encode()).digest()[:4], "big") >> 1


# Oracles whose cases are scored from committed results and never simulated here (layer 4, ADR 0002).
NOT_SIMULATED = {"simace_ts"}

SIMULATED: dict[tuple, Path] = {}  # one run's pedigrees by (seed, flags), shared across cases


def simulate(case: dict, seed: int) -> Path:
    """Run `simace simulate` for one seed and return the pedigree path.

    A pedigree depends only on its seed and flags, so every case of one run
    with the same seed and flags (in any order) reuses the first one's
    pedigree; the power study's mutations, which share seeds, do too.
    """
    flags = [f"--{k}={v}" for k, v in case["simulate"].items()]
    key = (seed, *sorted(flags))
    if key in SIMULATED:
        return SIMULATED[key]
    out = WORK / case["id"] / f"seed{seed}.parquet"
    out.parent.mkdir(parents=True, exist_ok=True)
    cmd = ["simace", "simulate", f"--seed={seed}", *flags, f"--output-pedigree={out}"]
    run(cmd, env=os.environ | SINGLE_THREAD)
    SIMULATED[key] = out
    return out


def run_oracle(case: dict, chunk: int, seeds: list[int], pedigrees: list[Path]) -> list[dict]:
    """Run the case's oracle over one chunk of pedigrees."""
    stem = WORK / case["id"] / f"oracle{chunk}"
    sim = case["simulate"]
    params = {k: sim[k] for k in ("A1", "C1", "E1", "A2", "C2", "E2", "rA", "rC", "rE")}
    for side, b in (("maternal", "bm"), ("paternal", "bf")):
        params[b] = [sim.get(f"Ctm-{side}{t}", 0) for t in (1, 2)]
    payload = {
        "params": params,
        "M": case["M"],
        "law": case["law"],
        "c_law": case.get("c_law", "realised"),
        "mutation": case.get("mutation"),
        "oracle_seed": case_seed(case["id"], f"oracle{chunk}"),
        "seeds": seeds,
        "pedigrees": [str(p) for p in pedigrees],
    }
    Path(f"{stem}.in.json").write_text(json.dumps(payload))
    script = ROOT / "oracles" / f"{case['oracle']}.R"
    run(["Rscript", str(script), f"{stem}.in.json", f"{stem}.out.json"], env=os.environ | SINGLE_THREAD)
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
    clipped = sum(r["c_clipped"] for r in per_seed)
    if clipped:
        print(f"{case['id']}: the oracle clipped U to PSD in {clipped} draw-generations", flush=True)
    by_stat: dict[str, list[dict]] = {}
    for rep in sorted(per_seed, key=lambda r: r["seed"]):
        for r in rep["rows"]:
            by_stat.setdefault(r["stat"], []).append(r)
    rows = []
    for stat, reps in by_stat.items():
        n_less = np.array([r["n_less"] for r in reps])
        n_equal = np.array([r["n_equal"] for r in reps])
        z = np.array([z_score(r) for r in reps])
        # One stream per statistic, so adding or removing a statistic re-jitters no other.
        rng = np.random.default_rng(case_seed(case["id"], f"pit:{stat}"))
        rows.append(row(case, stat, calibrate(randomised_pit(n_less, n_equal, case["M"], rng), z)))
    return rows


def genedrop_seed(case: dict, seed: int, pedigree: Path) -> dict[str, float]:
    """Layer 2, one seed: classify pairs and their kinship in R, then gene-drop."""
    pairs_path = pedigree.with_suffix(".pairs.parquet")
    f_path = pedigree.with_suffix(".F.parquet")
    run(
        ["Rscript", str(ROOT / "oracles" / "pairs.R"), str(pedigree), str(pairs_path), str(f_path)],
        env=os.environ | SINGLE_THREAD,
    )
    genome = Genome(**case["genome"])
    A = case["simulate"]["A1"]
    F = pl.read_parquet(f_path)["F"].to_numpy()
    return seed_stats(
        pl.read_parquet(pedigree), pl.read_parquet(pairs_path), F, genome, A, case_seed(case["id"], str(seed))
    )


def present(per_seed: list[dict], stat: str) -> list:
    """Return the stat's value in every seed that has it; a NaN value (too few pairs) counts as absent."""
    return [s[stat] for s in per_seed if stat in s and not np.any(np.isnan(s[stat]))]


def run_genedrop(case: dict, pool: ThreadPoolExecutor, pedigrees: list[Path], seeds: list[int]) -> list[dict]:
    """Layer 2: one gene-drop per pedigree; t-test each difference over the seeds that have it (d1_ rows INFO)."""
    per_seed = list(pool.map(lambda sp: genedrop_seed(case, *sp), zip(seeds, pedigrees, strict=True)))
    rows = []
    for stat in dict.fromkeys(k for s in per_seed for k in s):
        values = np.array(present(per_seed, stat))
        if not len(values):  # e.g. d1_ rows of generations with no inbred parents
            continue
        if stat.startswith("spread_"):
            sim, gd = values.T
            ratio = float(np.std(sim, ddof=1) / np.std(gd, ddof=1))
            rows.append(row(case, stat, Score(n=len(sim), center=ratio, effect=np.nan, p_values={}), info=True))
        else:
            rows.append(row(case, stat, mean_zero(values), info=stat.startswith("d1_")))
    return rows


def pedsim_check(case: dict, pool: ThreadPoolExecutor, pedigrees: list[Path], seeds: list[int]) -> list[dict]:
    """Layer 2 engine check: msprime against ped-sim IBD on fixed families (no simACE pedigrees)."""
    genome = Genome(**case["genome"], causal_per_chromosome=0)
    copies = case["copies"]
    work = WORK / case["id"]
    pedsim = pool.submit(pedsim_ibd, PEDSIM, work, genome, copies, case_seed(case["id"], "pedsim"))
    try:  # at most JOBS chromosomes in flight
        ms = msprime_ibd(genome, copies, case_seed(case["id"], "msprime"), partial(pool.map, buffersize=JOBS))
    except Exception as ms_exc:
        if (ped_exc := pedsim.exception()) is None:
            raise
        raise ExceptionGroup("msprime and ped-sim both failed", [ms_exc, ped_exc]) from None
    ped = pedsim.result()
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


def check_matched(case: dict) -> None:
    """Refuse a layer-5 case whose simACE flags do not match its xftsim spec (ADR 0003, decision 2)."""
    sim, x = case["simulate"], case_spec(case)
    r, r2, traits = x["r"], x["r2"], x["traits"]
    c = x["C"]
    want = {"N": x["N"], "G-ped": x["G"], "G-sim": x["G"], "mating-model": "standard", "mating-lambda": 0.5}
    want |= {"p-mztwin": 0, "A1": 0.5, "C1": c, "E1": 0.5 - c}
    want |= {"A2": 0.5, "C2": 0, "E2": 0.5, "rA": 0, "rC": 0, "rE": 0, "assort1": r}
    if ctm(x):
        want |= {"Ctm-maternal1": x["bm"], "Ctm-paternal1": x["bf"]}
    bad = []
    # am/run_xft.py takes transmission as bm > 0 or bf > 0, so a negative coefficient would split the sides
    if not (0 <= x["bm"] < 1 and 0 <= x["bf"] < 1):
        bad.append("xftsim bm and bf must lie in [0, 1)")
    # the half-sib regime pairs min(side sizes) couples, and xftsim halves an odd N unevenly, so N would shrink
    if x["pairing"] == "halfsib" and x["N"] % 2:
        bad.append("the halfsib pairing needs an even N")
    if traits == 2:
        cross = r if r2 == r else math.sqrt(r * r2)  # rank 1 (am/run_xft.py)
        want |= {"assort2": r2, "assort-matrix": json.dumps([[r, cross], [cross, r2]])}
    bad += [f"{k}={sim.get(k)!r} (want {v!r})" for k, v in want.items() if sim.get(k) != v]
    extra = set(sim) - set(want) - ({"assort2"} if traits == 1 else set())
    bad += [f"{k} is not modelled by the xftsim configuration" for k in sorted(extra)]
    if traits == 1 and sim.get("assort2", 0) != 0:
        bad.append("assort2 must be 0 for a single-trait case")
    if traits == 1 and r2 != r:
        bad.append("xftsim r2 needs two traits")
    if bad:
        raise ValueError(f"{case['id']}: simulate does not match the xftsim spec: {'; '.join(bad)}")


def ctm(spec: dict) -> bool:
    """Whether an xftsim spec has vertical C transmission (rwaples/simACE-validate#5)."""
    return spec["bm"] != 0 or spec["bf"] != 0


def bulmer_prediction(case: dict, stat: str) -> float | None:
    """The first-order expectation of a scored row, or None (INFO column of ADR 0003, decision 1).

    Bulmer's recursion without C transmission. With it, the random-mating law
    (`am.bulmer.predict_ctm`) at r = 0; under AM only Var(C), which simACE's realised top-up holds
    at C (simACE ADR 0025 §3).
    """
    x = case_spec(case)
    name, g = stat.rsplit("_g", 1)
    g = int(g)
    if ctm(x):
        law = predict_ctm(0.5, x["C"], 0.5 - x["C"], x["bm"], x["bf"])
        return law.get(name) if x["r"] == 0 or name == "var_C" else None
    if x["traits"] == 1:
        return predict1(0.5, 0.5, x["r"], g, c=x["C"]).get(name)
    r, r2 = x["r"], x["r2"]
    cross = r if r2 == r else math.sqrt(r * r2)
    return bulmer2(0.5, 0.5, [[r, cross], [cross, r2]], g)[g].get(name)


def run_xftsim(case: dict, pool: ThreadPoolExecutor, pedigrees: list[Path], seeds: list[int]) -> list[dict]:
    """Layer 5: simACE's per-generation moments against the committed xftsim reference, a Welch t per row.

    Power mutations (power.py) swap the reference for a mutated one; `null` reruns
    simACE on fresh seeds; `null_xft` puts a second xftsim sample in simACE's place.
    """
    check_matched(case)
    base_id = case["id"].split("~")[0]
    base = case | {"id": base_id}
    mutation = case.get("mutation")
    ref_name, spec = (base_id, case_spec(case)) if mutation in (None, "null", "null_xft") else derived(base, mutation)
    theirs = read_reference(ref_name, manifest(ref_name, spec))
    if mutation == "null_xft":
        ours_name, ours_spec = derived(base, "null_xft", int(case["id"].rsplit("~", 1)[1]))
        ours = read_reference(ours_name, manifest(ours_name, ours_spec))
    else:
        traits, with_c = spec["traits"], ctm(spec)
        ours = list(pool.map(lambda p: flatten(pedigree_rows(pl.read_parquet(p), traits, with_c)), pedigrees))
    x = case_spec(case)
    scored, info = schema(x["traits"], x["G"], ctm(x), x["pairing"] == "halfsib")
    rows = []
    for stat in scored + info:
        a = np.array([o[stat] for o in ours])
        b = np.array([t[stat] for t in theirs])
        rows.append(row(case, stat, welch(a, b), info=stat in info))
    if mutation is not None:
        return rows
    # INFO: each side against the first-order law; under C transmission the xftsim side too,
    # since its fixed household residual is the law simACE's realised top-up departs from (D5)
    prefix = "law" if ctm(x) else "bulmer"
    sides = (("", ours), ("xftsim:", theirs)) if ctm(x) else (("", ours),)
    for stat in scored:
        pred = bulmer_prediction(case, stat)
        if pred is None:
            continue
        for side, values in sides:
            a = np.array([o[stat] for o in values])
            rows.append(row(case, f"{prefix}:{side}{stat}", mean_zero(a - pred), info=True))
    return rows


def run_xftsim_finite_m(case: dict, pool: ThreadPoolExecutor, pedigrees: list[Path], seeds: list[int]) -> list[dict]:
    """Layer 5 sentinel (INFO): xftsim Var(A) at the smallest M against the largest, at each r (ADR 0003, decision 5)."""
    base_id = case["id"].split("~")[0]  # --replicate renames the case; the references keep the base id
    refs = {ref["label"]: ref["xftsim"] for ref in case["references"]}
    by_r: dict[float, dict[int, str]] = {}
    for label, x in refs.items():
        by_r.setdefault(x["r"], {})[x["M"]] = label
    rows = []
    for _r, labels in sorted(by_r.items()):
        lo, hi = labels[min(labels)], labels[max(labels)]
        a, b = (read_reference(f"{base_id}~{lab}", manifest(f"{base_id}~{lab}", refs[lab])) for lab in (lo, hi))
        for g in range(1, refs[lo]["G"] + 1):
            stat = f"var_A_g{g}"
            score = welch(np.array([x[stat] for x in a]), np.array([x[stat] for x in b]))
            rows.append(row(case, f"{lo}-{hi}:{stat}", score, info=True))
    return rows


RUNNERS = {
    "liability": run_liability,
    "genedrop": run_genedrop,
    "pedsim": pedsim_check,
    "slim_wf": run_structure,
    "slim_standard": run_structure,
    "mating_law": run_mating_law,
    "xftsim": run_xftsim,
    "xftsim_finite_m": run_xftsim_finite_m,
}
REPLACES = REPLACES_SIMACE | {"xftsim": ("null_xft",)}  # mutations whose simACE side is not simulated


def run_simace_ts(case: dict) -> list[dict]:
    """Layer 4: score the case's committed Parts (ADR 0002); nothing is simulated here."""
    parts = load_parts(case)
    values = case_rows(case, parts)
    rows = []
    for stat, gaps in values.scored.items():
        r = row(case, stat, mean_zero(gaps))
        r["gap_sd"] = float(np.std(gaps, ddof=1))
        rows.append(r)
    # An exact row's center is the largest |simACE - harness| and its effect the number of (seed, architecture)
    # entries that mismatch or are empty; only the scored architecture's can FAIL.
    for stat, (bad, worst) in values.exact.items():
        r = row(case, f"exact/{stat}", Score(n=len(parts), center=worst, effect=bad, p_values={}))
        r["exact_mismatches"] = bad
        rows.append(r)
    for stat, (bad, worst) in values.exact_sweep.items():
        rows.append(
            row(case, f"exact_sweep/{stat}", Score(n=len(parts), center=worst, effect=bad, p_values={}), info=True)
        )
    for stat, gaps in values.info_gaps.items():
        rows.append(row(case, stat, mean_zero(gaps), info=True))
    for stat, v in values.info_values.items():
        finite = v[np.isfinite(v)]
        center = float(finite.mean()) if len(finite) else np.nan
        rows.append(
            row(case, f"info/{stat}", Score(n=len(finite), center=center, effect=np.nan, p_values={}), info=True)
        )
    try:  # the layer-2 cross-check (ADR 0002 item 11); a bad G6 file must not cost the scored rows
        g6 = load_g6(case, parts)
    except (OSError, ValueError, KeyError) as exc:
        r = row(case, f"g6: {failure(exc)}", Score(n=len(parts), center=np.nan, effect=np.nan, p_values={}), info=True)
        r["error"] = True
        return [*rows, r]
    rows += [row(case, f"g6/{stat}", mean_zero(v), info=True) for stat, v in g6.items()]
    return rows


def run_case(case: dict) -> list[dict]:
    """Simulate the case's pedigrees and score them with its oracle."""
    if case["oracle"] in NOT_SIMULATED:
        return run_simace_ts(case)
    if case["oracle"] == "liability":
        unknown = set(case["simulate"]) - LIABILITY_FLAGS
        if unknown:
            raise ValueError(f"{case['id']}: the liability oracle does not model {', '.join(sorted(unknown))}")
    seed0 = case.get("seed0", 0)
    seeds = list(range(seed0, seed0 + case.get("R", 0)))
    replaced = case.get("mutation") in REPLACES.get(case["oracle"], ())
    (WORK / case["id"]).mkdir(parents=True, exist_ok=True)
    with ThreadPoolExecutor(JOBS) as pool:
        pedigrees = [] if replaced else list(pool.map(lambda s: simulate(case, s), seeds))
        return RUNNERS[case["oracle"]](case, pool, pedigrees, seeds)


def assign_status(rows: list[dict]) -> float:
    """Set each row's status against one Bonferroni family over every p-value; return alpha/K."""

    def scored(r: dict) -> bool:
        return not r["info"] and r["case_def"].get("expect", "pass") != "info"

    # An exact-recomputation row (layer 4) has no p-value: it fails on any mismatch and does not count in K.
    K = sum(len(r["score"].p_values) for r in rows if scored(r))
    threshold = ALPHA / max(K, 1)
    for r in rows:
        if "error" in r:
            r["status"] = "ERROR"
        elif "exact_mismatches" in r:
            r["status"] = "FAIL" if r["exact_mismatches"] else "PASS"
        elif not scored(r):
            r["status"] = "INFO"
        elif not all(np.isfinite(list(r["score"].p_values.values()))):
            r["status"] = "ERROR"  # a NaN p-value would compare False against the threshold and pass
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


def minimum_detectable(rows: list[dict], threshold: float) -> list[str]:
    """Return results.md's table of each layer-4 scored row's minimum detectable mean gap (ADR 0002 item 10)."""
    layer4 = [r for r in rows if "gap_sd" in r and r["status"] != "INFO" and r["score"].n >= 2]
    if not layer4:
        return []
    md = [
        "",
        "## Layer 4: minimum detectable gap",
        "",
        "The smallest mean gap (simACE minus its realised truth) each scored row rejects at alpha/K with",
        "80% power, from its per-seed SD (two-sided t). A measured power study needs repeated full sets of",
        "seeds, which is out of reach (ADR 0002 item 10).",
        "",
        "| Case | Statistic | n | Per-seed SD | Minimum detectable gap |",
        "|---|---|---|---|---|",
    ]
    for r in layer4:
        n = r["score"].n
        factor = (stats.t.ppf(1 - threshold / 2, n - 1) + stats.t.ppf(0.8, n - 1)) / np.sqrt(n)
        md.append(f"| {r['case']} | {r['stat']} | {n} | {r['gap_sd']:.4g} | {factor * r['gap_sd']:.4g} |")
    return md


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
        "is the SD ratio simACE / gene-drop over seeds. Layer 4: a scored row is like layer 2 (simACE minus its",
        "realised truth); an exact row has no p-value and is FAIL when any entry mismatches, with center the",
        "largest difference and effect the number of mismatching (seed, architecture) entries.",
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
        did = case.get("id_diff", "") or ("D1" if r["stat"].removeprefix("spread_").startswith("d1_") else "")
        pv = ", ".join(f"{k} {v:.3g}" for k, v in c.p_values.items())
        vals = [r["case"], r["layer"], r["stat"], c.n, f"{c.center:.4g}", f"{c.effect:+.2f}", pv]
        lines.append("\t".join(map(str, [*vals, f"{threshold:.3g}", r["status"], did])))
        md.append("| " + " | ".join(map(str, [*vals, r["status"], did])) + " |")
    for case, message in errors:
        lines.append(
            "\t".join(map(str, [case["id"], case["layer"], message, "", "", "", "", f"{threshold:.3g}", "ERROR", ""]))
        )
        md.append(f"| {case['id']} | {case['layer']} | {message} | | | | | ERROR | |")
    md += minimum_detectable(rows, threshold)
    (out / "results.tsv").write_text("\n".join(lines) + "\n")
    (out / "results.md").write_text("\n".join(md) + "\n")
    print(f"Totals: {total_line}; written to {out.relative_to(ROOT)}/")


def load_cases(parser: argparse.ArgumentParser, ids: str | None = None) -> list[dict]:
    """Return cases.yaml's cases, or those of the comma-separated ids; an unknown id is a usage error."""
    cases = yaml.safe_load((ROOT / "cases.yaml").read_text())["cases"]
    if ids is None:
        return cases
    wanted = set(ids.split(","))
    unknown = wanted - {c["id"] for c in cases}
    if unknown:
        parser.error(f"unknown case ids: {', '.join(sorted(unknown))}; cases: {', '.join(c['id'] for c in cases)}")
    return [c for c in cases if c["id"] in wanted]


def failure(exc: BaseException) -> str:
    """Return the ERROR row's message for exc; a group names each of its failures."""
    if isinstance(exc, BaseExceptionGroup):
        return "; ".join(failure(e) for e in exc.exceptions)
    if isinstance(exc, subprocess.CalledProcessError):
        return f"{' '.join(map(str, exc.cmd[:2]))} exited {exc.returncode}"
    return f"{type(exc).__name__}: {exc}"


def print_failure(exc: BaseException) -> None:
    """Print a failed command's stderr tail, or else the traceback, for exc and each failure in a group."""
    if isinstance(exc, BaseExceptionGroup):
        for e in exc.exceptions:
            print_failure(e)
    elif isinstance(exc, subprocess.CalledProcessError):
        print(exc.stderr.decode(errors="replace")[-2000:], file=sys.stderr)
    else:
        traceback.print_exception(exc)


def main() -> int:
    """Run the selected cases."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--only", help="comma-separated case ids")
    parser.add_argument(
        "--replicate", type=int, default=0, help="shift every seed (simACE and oracle) to an independent replicate run"
    )
    args = parser.parse_args()
    cases = load_cases(parser, args.only)
    if args.replicate:
        # A case scored from committed results has no other seeds to replicate.
        skipped = [c["id"] for c in cases if c["oracle"] in NOT_SIMULATED]
        if skipped:
            print(f"--replicate skips {', '.join(skipped)}: scored from committed Parts, not simulated", flush=True)
        cases = [c for c in cases if c["oracle"] not in NOT_SIMULATED]
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
        except Exception as exc:  # one broken case is an ERROR row, not a lost run
            errors.append((case, failure(exc)))
            print_failure(exc)
    threshold = assign_status(rows)
    write_results(rows, threshold, errors, selection)
    return 1 if errors or any(r["status"] in ("FAIL", "ERROR") for r in rows) else 0


if __name__ == "__main__":
    sys.exit(main())
