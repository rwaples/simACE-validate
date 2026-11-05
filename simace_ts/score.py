"""Score a simACE-ts Case from its committed Parts (ADR 0002 items 4, 9 and 10).

Parts are ``simace_ts/parts/<case id>/seed<s>.json``, one per seed in
``seed0 .. seed0 + R - 1``, each the ``part.json`` that ``run_seed.py`` wrote
on GenomeDK. ``load_parts`` refuses (raises, so the Case is ERROR) when:

- a seed has no Part, or a Part's manifest names another seed;
- a Part's case differs from cases.yaml in a key its outputs depend on
  (``OUTPUT_KEYS``), its chromosomes are not all 22, its data pins are not
  the case's, or its simace is not the one this run has installed;
- a Part was computed by another ``PART_VERSION`` of part.py (bumped with
  any change to a truth formula or to the statistics), or by a checkout
  with uncommitted changes, or its seed was run from one;
- a Part's scored statistics are not the case's schema (``scored_schema``),
  or its exact statistics differ from the other Parts'.

``case_rows`` then returns, for the first (scored) architecture, one
per-seed gap vector per scored row (simACE minus its realised truth) and
its exact rows, and the INFO rows: the other architectures' gaps and exact
checks, and the first architecture's INFO statistics. A scored value
missing in any seed raises; a missing INFO value is NaN.
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
CHROMOSOMES = tuple(range(1, 23))
PAIR_CLASSES = ("FS", "MHS", "PHS", "MO", "FO")
PART_VERSION = 1  # bump with any change to part.py's truth formulas or statistics; older Parts are then stale


def po_generations(case: dict) -> list[int]:
    """Return the child generations whose parents are phenotyped too."""
    g_ped, g_pheno = int(case["simulate"]["G-ped"]), int(case["cohort"]["G-pheno"])
    return list(range(g_ped - g_pheno + 1, g_ped))


def scored_schema(case: dict) -> list[str]:
    """Return the scored statistics every Part of ``case`` must carry (ADR 0002 item 9: K is fixed by the schema)."""
    po = [f"po_g{g}_{s}" for g in po_generations(case) for s in ("slope", "r")]
    return [*(f"lc_{k}" for k in PAIR_CLASSES), *po, "dz_L", "po_L_slope"]


@dataclass(frozen=True)
class CaseRows:
    """The per-seed values a simACE-ts Case is scored from."""

    scored: dict[str, np.ndarray]  # stat -> gap (simACE - truth) per seed
    exact: dict[str, tuple[int, float]]  # stat -> (seeds that mismatch or are empty, largest |difference|)
    exact_sweep: dict[str, tuple[int, float]]  # the same over the sweep architectures (INFO)
    info_gaps: dict[str, np.ndarray]  # "<architecture>/<stat>" -> gap per seed
    info_values: dict[str, np.ndarray]  # stat -> value per seed (NaN where simACE left it empty)


def installed_simace() -> dict:
    """Return the installed simace's version and install source; run_seed.py records this in each manifest."""
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
        if m["chromosomes"] != list(CHROMOSOMES):
            problems.append("chromosomes")
        if m["data"] != case["data"]:
            problems.append("data pins")
        if m["simace"] != simace:
            problems.append(f"simace {m['simace']['version']} (this run has {simace['version']})")
        if part.get("part_version") != PART_VERSION:
            problems.append(f"part_version {part.get('part_version')} (this run has {PART_VERSION})")
        commits = [part.get("part_by", "-dirty"), *(r["simace_validate"] for r in m["runs"])]
        if any(c.endswith("-dirty") for c in commits):
            problems.append("computed from uncommitted code")
        if problems:
            raise ValueError(f"{case['id']} seed {s}: stale Part ({', '.join(problems)}); rerun it on GenomeDK")
        parts.append(part)
    return parts


def case_rows(case: dict, parts: list[dict]) -> CaseRows:
    """Return the scored gaps, exact summaries and INFO values of the case's Parts."""
    names = [a["name"] for a in case["architectures"]]
    first = names[0]
    schema = scored_schema(case)
    exact_names = list(parts[0]["architectures"][first]["exact"])
    for p in parts:
        for arch in names:
            got = p["architectures"][arch]
            where = f"{case['id']} seed {p['manifest']['seed']} {arch}"
            if set(got["scored"]) != set(schema):
                raise ValueError(f"{where}: scored statistics {sorted(got['scored'])} are not the schema {schema}")
            if set(got["exact"]) != set(exact_names):
                raise ValueError(f"{where}: exact statistics differ from seed {parts[0]['manifest']['seed']}'s")

    def gaps(arch: str, stat: str) -> np.ndarray:
        out = []
        for p in parts:
            row = p["architectures"][arch]["scored"][stat]
            if row["simace"] is None or row["truth"] is None:
                raise ValueError(f"{case['id']} seed {p['manifest']['seed']}: {arch}/{stat} is empty")
            out.append(row["simace"] - row["truth"])
        return np.array(out)

    def exact_summary(archs: list[str], stat: str) -> tuple[int, float]:
        bad, worst = 0, 0.0
        for p in parts:
            for arch in archs:
                e = p["architectures"][arch]["exact"][stat]
                if e["simace"] is None or e["harness"] is None:
                    bad += 1
                    continue
                diff = abs(e["simace"] - e["harness"])
                worst = max(worst, diff)
                bad += diff > e["tol"]
        return bad, worst

    def info_value(p: dict, stat: str) -> float:
        v = p["architectures"][first]["info"].get(stat)
        return np.nan if v is None else v

    info_names = dict.fromkeys(k for p in parts for k in p["architectures"][first]["info"])
    return CaseRows(
        scored={stat: gaps(first, stat) for stat in schema},
        exact={stat: exact_summary([first], stat) for stat in exact_names},
        exact_sweep={stat: exact_summary(names[1:], stat) for stat in exact_names} if names[1:] else {},
        info_gaps={f"{arch}/{stat}": gaps(arch, stat) for arch in names[1:] for stat in schema},
        info_values={stat: np.array([info_value(p, stat) for p in parts]) for stat in info_names},
    )
