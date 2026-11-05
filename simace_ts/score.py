"""Score a simACE-ts Case from its committed Parts (ADR 0002 items 4, 9 and 10).

Parts are ``simace_ts/parts/<case id>/seed<s>.json``, one per seed in
``seed0 .. seed0 + R - 1``, each the ``part.json`` that ``run_seed.py`` wrote
on GenomeDK. ``load_parts`` refuses (raises, so the Case is ERROR) when:

- a seed has no Part, or a Part's manifest names another seed;
- a Part's case differs from cases.yaml in a key its outputs depend on
  (``OUTPUT_KEYS``), its chromosomes are not all 22, its data pins are not
  the case's, or its simace is not the one this run has installed.

``case_rows`` then returns, for the first (scored) architecture, one
per-seed gap vector per scored row (simACE minus its realised truth), the
exact rows, and the INFO rows: the other architectures' gaps and the first
architecture's INFO statistics. A scored value missing in any seed raises.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from importlib import metadata
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
PARTS = ROOT / "simace_ts" / "parts"
OUTPUT_KEYS = ("id", "simulate", "cohort", "architectures")  # what a seed's outputs depend on (run_seed.py too)
CHROMOSOMES = list(range(1, 23))


@dataclass(frozen=True)
class CaseRows:
    """The per-seed values a simACE-ts Case is scored from."""

    scored: dict[str, np.ndarray]  # stat -> gap (simACE - truth) per seed
    exact: dict[str, tuple[int, float]]  # stat -> (seeds that mismatch, largest |difference|)
    info_gaps: dict[str, np.ndarray]  # "<architecture>/<stat>" -> gap per seed
    info_values: dict[str, np.ndarray]  # stat -> value per seed (NaN where simACE left it empty)


def installed_simace() -> dict:
    """Return the installed simace's version and install source, as run_seed.py records them."""
    dist = metadata.distribution("simace")
    direct = dist.read_text("direct_url.json")
    return {"version": dist.version, "direct_url": json.loads(direct) if direct else None}


def load_parts(case: dict, parts_dir: Path = PARTS, simace: dict | None = None) -> list[dict]:
    """Return the case's Parts in seed order, refusing any that is missing or stale."""
    simace = installed_simace() if simace is None else simace
    seeds = range(case["seed0"], case["seed0"] + case["R"])
    missing = [s for s in seeds if not (parts_dir / case["id"] / f"seed{s}.json").exists()]
    if missing:
        raise FileNotFoundError(
            f"{case['id']}: {len(missing)} of {len(seeds)} Parts missing (first: seed {missing[0]})"
        )
    parts = []
    for s in seeds:
        part = json.loads((parts_dir / case["id"] / f"seed{s}.json").read_text())
        m = part["manifest"]
        problems = [k for k in OUTPUT_KEYS if m["case"].get(k) != case.get(k)]
        if m["seed"] != s:
            problems.append(f"seed {m['seed']}")
        if m["chromosomes"] != CHROMOSOMES:
            problems.append("chromosomes")
        if m["data"] != case["data"]:
            problems.append("data pins")
        if m["simace"] != simace:
            problems.append(f"simace {m['simace']['version']} (this run has {simace['version']})")
        if problems:
            raise ValueError(f"{case['id']} seed {s}: stale Part ({', '.join(problems)}); rerun it on GenomeDK")
        parts.append(part)
    return parts


def case_rows(case: dict, parts: list[dict]) -> CaseRows:
    """Return the scored gaps, exact summaries and INFO values of the case's Parts."""
    names = [a["name"] for a in case["architectures"]]
    first = names[0]

    def gaps(arch: str, stat: str) -> np.ndarray:
        out = []
        for p in parts:
            row = p["architectures"][arch]["scored"][stat]
            if row["simace"] is None or row["truth"] is None:
                raise ValueError(f"{case['id']} seed {p['manifest']['seed']}: {arch}/{stat} is empty")
            out.append(row["simace"] - row["truth"])
        return np.array(out)

    scored_stats = list(parts[0]["architectures"][first]["scored"])
    scored = {stat: gaps(first, stat) for stat in scored_stats}
    info_gaps = {f"{arch}/{stat}": gaps(arch, stat) for arch in names[1:] for stat in scored_stats}

    exact: dict[str, tuple[int, float]] = {}
    for stat in parts[0]["architectures"][first]["exact"]:
        bad, worst = 0, 0.0
        for p in parts:
            for arch in names:
                e = p["architectures"][arch]["exact"][stat]
                if e["simace"] is None or e["harness"] is None:
                    bad += 1
                    continue
                diff = abs(e["simace"] - e["harness"])
                worst = max(worst, diff)
                bad += diff > e["tol"]
        exact[stat] = (bad, worst)

    info_values = {
        stat: np.array([np.nan if (v := p["architectures"][first]["info"][stat]) is None else v for p in parts])
        for stat in parts[0]["architectures"][first]["info"]
    }
    return CaseRows(scored=scored, exact=exact, info_gaps=info_gaps, info_values=info_values)
