"""Layer 5's reference samples: specs, manifests, files and the scored-row schema (ADR 0003, decision 8).

A reference sample is R xftsim runs at one configuration, made once by
`pixi run -e xftsim xft-reference` (am/reference.py) and committed under
am/reference/. The harness reads it and checks its manifest against the one
it expects, so a stale sample is an ERROR, never a silent comparison.

Runs in both envs (Python 3.9 xftsim, Python 3.14 default): stdlib and yaml only.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
REFERENCE_DIR = ROOT / "am" / "reference"
NUMBA_THREADS = "1"  # xftsim output depends on the thread count (ADR 0003, decision 6)
FORMAT = 1

# Power mutations of the reference side (ADR 0003, decision 8): each maps the case's
# xftsim spec to the mutated one; applicable(case_spec) says where it makes sense.
XFT_MUTATIONS = {
    "r_090": (lambda x: {**x, "r": round(0.9 * x["r"], 10), "r2": round(0.9 * x["r2"], 10)}, lambda x: x["r"] != 0),
    "r_0": (lambda x: {**x, "r": 0.0, "r2": 0.0}, lambda x: x["r"] != 0),
    "mate_A": (lambda x: {**x, "mate": "additive"}, lambda x: x["r"] != 0),
    "mate_y1": (lambda x: {**x, "mate": "y1"}, lambda x: x["traits"] == 2 and x["r"] != 0 and x["r2"] == x["r"]),
}


def derived_seed0(name: str) -> int:
    """Seed range start for a derived reference (mutation or null), from its name; 31-bit, as run_harness.case_seed."""
    return int.from_bytes(hashlib.sha256(f"{name}:xftsim".encode()).digest()[:4], "big") >> 1


def defaults(xftsim: dict) -> dict:
    """An xftsim spec with defaults filled in; r2 (trait 2's mate correlation, two traits) defaults to r."""
    return {"traits": 1, "mate": "phenotype", "C": 0, "r2": xftsim["r"], **xftsim}


def case_spec(case: dict) -> dict:
    """A scored case's own xftsim spec, defaults filled in."""
    return defaults(case["xftsim"])


def derived(case: dict, mutation: str, rep: int = 0) -> tuple[str, dict]:
    """(name, spec) of a power reference: a mutation of the case's spec, or the rep-th null_xft sample."""
    base = case_spec(case)
    if mutation == "null_xft":
        name, spec = f"{case['id']}~null_xft{rep}", base
    else:
        name, spec = f"{case['id']}~{mutation}", XFT_MUTATIONS[mutation][0](base)
    return name, {**spec, "seed0": derived_seed0(name)}


def sha256_file(path: Path) -> str:
    """Hex sha256 of a file."""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def xftsim_lock() -> tuple[str, str]:
    """Return (sha256 of the xftsim env's lock entries, the xftsim git commit) from pixi.lock."""
    env = yaml.safe_load((ROOT / "pixi.lock").read_text())["environments"]["xftsim"]
    blob = json.dumps(env, sort_keys=True).encode()
    pins = [p["pypi"] for p in env["packages"]["linux-64"] if "pypi" in p and "xftsim" in p["pypi"]]
    commit = pins[0].rsplit("#", 1)[-1] if pins else ""
    return hashlib.sha256(blob).hexdigest(), commit


def full_spec(xftsim: dict) -> dict:
    """The case's xftsim block with defaults filled in, in a fixed key order."""
    keys = ("N", "M", "G", "r", "r2", "traits", "mate", "C", "seed0", "R")
    x = defaults(xftsim)
    return {k: x[k] for k in keys}


def r_arg(spec: dict) -> str:
    """The R argument of am/run_xft.py for a spec: r, or R1/R2 when the two traits' r differ."""
    return str(spec["r"]) if spec["r2"] == spec["r"] else f"{spec['r']}/{spec['r2']}"


def manifest(name: str, xftsim: dict) -> dict:
    """The manifest a reference sample named `name` at spec `xftsim` must carry."""
    lock, commit = xftsim_lock()
    return {
        "format": FORMAT,
        "name": name,
        "spec": full_spec(xftsim),
        "xftsim_commit": commit,
        "xftsim_lock_sha256": lock,
        # The code that builds a sample: the run, its statistics, and the driver (argv, threads, flatten).
        "stats_sha256": sha256_file(ROOT / "am" / "stats.py"),
        "run_xft_sha256": sha256_file(ROOT / "am" / "run_xft.py"),
        "reference_sha256": sha256_file(ROOT / "am" / "reference.py"),
        "numba_threads": NUMBA_THREADS,
    }


def reference_path(name: str) -> Path:
    """Where the reference sample `name` lives."""
    return REFERENCE_DIR / f"{name}.jsonl"


def read_reference(name: str, expected: dict) -> list[dict[str, float]]:
    """Return the reference's per-seed statistics after checking its manifest and seeds; raise if stale."""
    path = reference_path(name)
    if not path.exists():
        raise ValueError(f"reference {name} is missing; run `pixi run -e xftsim xft-reference`")
    lines = path.read_text().splitlines()
    got = json.loads(lines[0])["manifest"]
    stale = sorted(k for k in set(got) | set(expected) if got.get(k) != expected.get(k))
    if stale:
        raise ValueError(f"reference {name} is stale ({', '.join(stale)}); rerun `pixi run -e xftsim xft-reference`")
    records = [json.loads(line) for line in lines[1:]]
    spec = expected["spec"]
    want = list(range(spec["seed0"], spec["seed0"] + spec["R"]))
    if [r["seed"] for r in records] != want:
        raise ValueError(f"reference {name} does not hold seeds {want[0]}..{want[-1]} exactly once each, in order")
    return [r["stats"] for r in records]


# The scored statistics (ADR 0003, decision 8): per trait, Var(A) and Var(L) at g = 1..G,
# and six correlations at g = 2..G; two traits add the cross-trait rows.
PER_TRAIT_ALL_G = ("var_A", "var_L")
PER_TRAIT = ("r_mate_L", "r_mate_A", "r_FS_A", "r_FS_L", "r_PO_A", "r_PO_L")
PER_TRAIT_INFO = ("r_MO_L", "r_FO_L")  # xftsim's "maternal" is a random label
CROSS = ("r_w_A", "r_w_L", "r_mate_L12", "r_mate_L21", "r_mate_A12", "r_mate_A21")


def schema(traits: int, G: int) -> tuple[list[str], list[str]]:
    """Return (scored, info) row names for a case."""
    suffixes = [""] if traits == 1 else ["", "_t2"]
    scored, info = [], []
    for t in suffixes:
        scored += [f"{s}{t}_g{g}" for s in PER_TRAIT_ALL_G for g in range(1, G + 1)]
        scored += [f"{s}{t}_g{g}" for s in PER_TRAIT for g in range(2, G + 1)]
        info += [f"{s}{t}_g{g}" for s in PER_TRAIT_INFO for g in range(2, G + 1)]
    if traits == 2:
        scored += [f"{s}_g{g}" for s in CROSS for g in range(2, G + 1)]
    return scored, info
