"""Make layer 5's xftsim reference samples (ADR 0003, decision 8); runs in the xftsim env.

    pixi run -e xftsim xft-reference                       # every layer-5 case's reference
    pixi run -e xftsim xft-reference --case am_r05          # one case
    pixi run -e xftsim xft-reference --case am_r05 --mutation r_090,mate_A,null_xft --reps 10

A reference is R runs of am/run_xft.py at a case's `xftsim` spec, written to
am/reference/<name>.jsonl: a manifest line, then one line of flattened
statistics per seed. A reference whose file already carries the expected
manifest is skipped, so reruns converge; a stale one is rebuilt.

Names: <case> for the case itself; <case>~<label> for a finite-M sentinel's
references; <case>~<mutation> for a power mutation; <case>~null_xft<rep> for
the power null that puts a second xftsim sample in simACE's place. Derived
references take their seeds from their name (manifest.derived_seed0), never
from a range a scored reference uses.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

import yaml

from am.manifest import (
    NUMBA_THREADS,
    ROOT,
    XFT_MUTATIONS,
    case_spec,
    derived,
    flatten,
    manifest,
    reference_path,
)

ENV = {"OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1", "MKL_NUM_THREADS": "1", "NUMBA_NUM_THREADS": NUMBA_THREADS}


def specs(case: dict, mutations: list[str], reps: int) -> list[tuple[str, dict]]:
    """Return (name, xftsim spec) for every reference the case needs."""
    if case["oracle"] == "xftsim_finite_m":
        return [(f"{case['id']}~{ref['label']}", ref["xftsim"]) for ref in case["references"]]
    out = [(case["id"], case_spec(case))]
    for m in mutations:
        if m == "null_xft":
            out += [derived(case, m, rep) for rep in range(1, reps + 1)]
        elif XFT_MUTATIONS[m][1](case_spec(case)):
            out.append(derived(case, m))
    return out


def current(name: str, expected: dict) -> bool:
    """Whether the reference file exists and carries the expected manifest."""
    path = reference_path(name)
    if not path.exists():
        return False
    with path.open() as f:
        return json.loads(f.readline())["manifest"] == expected


def run_one(spec: dict, seed: int) -> dict:
    """One xftsim run of the spec at `seed`, as {seed, stats}."""
    args = [str(a) for a in (seed, spec["N"], spec["M"], spec["G"], spec["r"], spec["traits"], spec["mate"], spec["C"])]
    out = subprocess.run(
        [sys.executable, "-m", "am.run_xft", *args],
        check=True,
        capture_output=True,
        text=True,
        cwd=ROOT,
        env={**os.environ, **ENV},
    )
    rec = json.loads(out.stdout.strip().splitlines()[-1])
    return {"seed": seed, "stats": flatten(rec["rows"])}


def build(name: str, xftsim: dict, pool: ThreadPoolExecutor) -> None:
    """Run the spec's R seeds and write the reference atomically."""
    expected = manifest(name, xftsim)
    spec = expected["spec"]
    seeds = list(range(spec["seed0"], spec["seed0"] + spec["R"]))
    records = list(pool.map(lambda s: run_one(spec, s), seeds))
    path = reference_path(name)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    lines = [json.dumps({"manifest": expected})] + [json.dumps(r) for r in records]
    tmp.write_text("\n".join(lines) + "\n")
    tmp.replace(path)
    print(f"{name}: {len(records)} runs written to {path.relative_to(ROOT)}", flush=True)


def main() -> None:
    """Build the selected references."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--case", help="comma-separated case ids (default: every layer-5 case)")
    parser.add_argument("--mutation", default="", help="comma-separated power mutations (and null_xft)")
    parser.add_argument("--reps", type=int, default=10, help="null_xft replicates")
    parser.add_argument("--jobs", type=int, default=6)
    args = parser.parse_args()
    cases = [c for c in yaml.safe_load((ROOT / "cases.yaml").read_text())["cases"] if c["layer"] == 5]
    if args.case:
        wanted = set(args.case.split(","))
        unknown = wanted - {c["id"] for c in cases}
        if unknown:
            parser.error(f"not layer-5 cases: {', '.join(sorted(unknown))}")
        cases = [c for c in cases if c["id"] in wanted]
    mutations = [m for m in args.mutation.split(",") if m]
    unknown = set(mutations) - set(XFT_MUTATIONS) - {"null_xft"}
    if unknown:
        parser.error(f"unknown mutations: {', '.join(sorted(unknown))}")
    with ThreadPoolExecutor(args.jobs) as pool:
        for case in cases:
            for name, xftsim in specs(case, mutations, args.reps):
                if current(name, manifest(name, xftsim)):
                    print(f"{name}: current, skipped", flush=True)
                    continue
                build(name, xftsim, pool)


if __name__ == "__main__":
    main()
