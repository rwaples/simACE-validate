"""Power of simace validate's inheritance tests: verdict rejection rates at stated errors of the scored parameters.

    pixi run power-inheritance --cell rm_base --study-reps 1-10

One study replicate of a cell is one scenario of 20 reps (the tests' min_reps)
built by ``simace run``. Each mutation scores those reps as if one scenario
parameter in their ``report.yaml`` were wrong, which is the only place
``household_draw`` and ``mendelian_residual`` read it (simace
``analysis/validate/calibrated.py`` ``RepInputs.params``). The rep's own
``report.yaml`` is fingerprinted in its ``run.yaml``, so the mutation goes to a
copy: each rep runs simace's per-rep ``simace validate --pedigree ...``
primitive in-process (``analysis/validate/runner.cli``), with the argv the
scenario command builds (``cli/stages.validate_argv``) and ``--report`` pointing
at the copy. The reps' sections are then decided by simace's own
``decide(scenario_family(...))``, the scenario verdict ``simace validate
<scenario>`` writes. A study replicate rejects when that verdict is ``fail``.

The ``null`` mutation scores the true parameters; for it the real
``simace validate <scenario>`` also runs, and its ``validation_summary.yaml``
must match the in-process decision record for record, or the run stops.

To fan out, give each job --cell, --study-reps and --parts; then --merge the
part directory. Seeds depend only on the cell and the study replicate, so the
merged table equals one run's. Writes results/power_inheritance.md when the
table is every cell at the default study reps, else .work/results/.
"""

from __future__ import annotations

import argparse
import copy
import json
import shutil
import tempfile
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from power import interval
from run_harness import PARTIAL_RESULTS, RESULTS, case_seed, simace_ref
from subproc import ToolError, run

#: simace v0.3.0's ``config/_default.yaml`` (git blob de155750a2a7), which the wheel does not ship.
CONFIG = Path(__file__).parent / "config"
FOLDER = "inhpow"
#: `simace validate <scenario>`'s exit status for a ``fail`` verdict (simace cli/validate.py).
VERDICT_FAIL_EXIT = 3
REPS_PER_VERDICT = 20  # min_reps of both tests (simace analysis/validate/on_demand.py)
DEFAULT_STUDY_REPS = 100
DELTAS = (0.05, 0.1, 0.2, 0.3, 0.5)

_LTM = {"model": "adult", "params": {"method": "ltm", "cip_x0": 50, "cip_k": 0.2, "prevalence": 0.10}, "beta": 1.0}

#: rm_base of cases.yaml, as a simace scenario; LTM phenotypes, so onset_tau does not apply.
BASE: dict[str, Any] = {
    "replicates": REPS_PER_VERDICT,
    "N": 2000,
    "G_ped": 4,
    "G_sim": 4,
    "pedigree": {
        "mating_lambda": 0.5,
        "p_mztwin": 0.02,
        "rA": 0.5,
        "rC": 0.3,
        "rE": 0.1,
        "trait1": {"A": 0.5, "C": 0.2, "E": 0.3},
        "trait2": {"A": 0.4, "C": 0.1, "E": 0.5},
    },
    "phenotype": {"trait1": _LTM, "trait2": _LTM},
}

#: Each cell changes BASE in one respect; ctm_base carries cases.yaml's ctm_base transmission.
CELLS: dict[str, tuple[str, dict[str, Any]]] = {
    "rm_base": ("rm_base: N 2000, G_ped 4, mating_lambda 0.5, random mating, no C transmission", {}),
    "ctm_base": (
        "rm_base with b_m, b_f 0.5, 0.3 (trait 1) and 0.3, 0.5 (trait 2)",
        {
            "pedigree": {
                "trait1": {"C_transmission": {"maternal": 0.5, "paternal": 0.3}},
                "trait2": {"C_transmission": {"maternal": 0.3, "paternal": 0.5}},
            }
        },
    ),
    "rm_n1000": ("rm_base at N 1000", {"N": 1000}),
    "rm_n20000": ("rm_base at N 20000", {"N": 20000}),
    "rm_gped3": ("rm_base at G_ped = G_sim = 3", {"G_ped": 3, "G_sim": 3}),
    "rm_gped6": ("rm_base at G_ped = G_sim = 6", {"G_ped": 6, "G_sim": 6}),
    "rm_lambda3": ("rm_base at mating_lambda 3", {"pedigree": {"mating_lambda": 3.0}}),
    "rm_am": ("rm_base with assortative mating 0.3 on both traits", {"pedigree": {"assort1": 0.3, "assort2": 0.3}}),
}


@dataclass(frozen=True)
class Mutation:
    """One scored parameter of trait 1 off by ``delta``: added for a transmission coefficient, a share of A."""

    key: str | None
    delta: float = 0.0
    #: The test that reads ``key``; None for the null.
    test: str | None = None

    def apply(self, params: dict[str, Any]) -> dict[str, Any]:
        """Return the scenario parameters as this mutation scores them."""
        if self.key is None:
            return params
        value = float(params[self.key])
        return params | {self.key: value * (1 + self.delta) if self.key == "A1" else value + self.delta}

    @property
    def effect(self) -> str:
        """The mutation in words, for the table."""
        if self.key is None:
            return "true parameters (false-positive rate of the verdict)"
        if self.key == "A1":
            return f"A1 scored as {1 + self.delta:g}*A1: variance cell expects {1 + self.delta:g}*A1/2"
        return f"{self.key} scored {self.delta:+g} off"


def _mutations() -> dict[str, Mutation]:
    out = {"null": Mutation(None)}
    for name, key, test in (
        ("bf", "Ctm_paternal1", "household_draw"),
        ("bm", "Ctm_maternal1", "household_draw"),
        ("a", "A1", "mendelian_residual"),
    ):
        for sign in (1, -1):
            for d in DELTAS:
                out[f"{name}_{'plus' if sign > 0 else 'minus'}_{round(d * 100):03d}"] = Mutation(key, sign * d, test)
    return out


MUTATIONS = _mutations()


def merged(base: dict[str, Any], change: dict[str, Any]) -> dict[str, Any]:
    """``base`` with ``change`` merged in, nested dicts key by key."""
    out = copy.deepcopy(base)
    for key, value in change.items():
        out[key] = merged(out[key], value) if isinstance(value, dict) and isinstance(out.get(key), dict) else value
    return out


def study_seed(cell: str, study_rep: int) -> int:
    """The first rep seed of a study replicate; its 20 reps take seed..seed+19, below 2**31."""
    return case_seed(cell, f"inheritance{study_rep}") >> 1


def failed_records(decided: dict[str, Any]) -> list[str]:
    """The names of a decided family's failed records, such as ``household_draw/trait1_PHS/mean_z``."""
    from simace.analysis.validate.verdict import record_name, records

    return sorted(record_name(path) for path, rec in records(decided) if rec.get("passed") is False)


def _cli_verdict(scenario: str, roots: list[str], results: Path) -> dict[str, Any]:
    """Run the real ``simace validate <scenario>`` and return its verdict; it exits 3 on a ``fail`` verdict."""
    from simace.analysis.validate.verdict import load_verdict

    try:
        run(["simace", "validate", scenario, "--no-plots", *roots])
    except ToolError as exc:
        if exc.returncode != VERDICT_FAIL_EXIT:
            raise
    loaded = load_verdict(results / FOLDER / scenario / "validation_summary.yaml")
    if loaded is None:
        raise RuntimeError(f"simace validate {scenario} wrote no verdict")
    return loaded


def study_replicate(cell: str, study_rep: int, names: list[str], work: Path, jobs: int) -> dict[str, dict[str, Any]]:
    """Build one study replicate of a cell and return, per mutation, its verdict and failed records."""
    from simace.analysis.validate.on_demand import ON_DEMAND_CHECKS, load_validation
    from simace.analysis.validate.runner import cli as rep_cli
    from simace.analysis.validate.verdict import decide, scenario_family
    from simace.cli.layout import Layout, RepArtifact
    from simace.cli.run import plan_runs
    from simace.cli.stages import validate_argv
    from simace.core.yaml_io import to_native

    scenario = f"{cell}_s{study_rep}"
    config_dir, results, logs = work / "config", work / "results", work / "logs"
    config_dir.mkdir(parents=True)
    shutil.copy(CONFIG / "_default.yaml", config_dir)
    spec = merged(BASE, CELLS[cell][1]) | {"seed": study_seed(cell, study_rep)}
    (config_dir / f"{FOLDER}.yaml").write_text(yaml.safe_dump({scenario: spec}, sort_keys=False))
    roots = ["--config-dir", str(config_dir), "--results", str(results), "--logs", str(logs), "--jobs", str(jobs)]
    run(["simace", "run", scenario, "--no-plots", *roots])

    layout = Layout(root=results, logs=logs)
    reps = plan_runs("validate", config_dir, [scenario], None)[0].all_reps
    checks = list(ON_DEMAND_CHECKS)
    # RepInputs.params reads only report.yaml's inputs.parameters (simace analysis/validate/calibrated.py).
    params = [
        yaml.load(rep.path(layout, RepArtifact.REPORT).read_text(), Loader=yaml.CSafeLoader)["inputs"]["parameters"]
        for rep in reps
    ]
    out = {}
    for name in names:
        sections: dict[str, list[dict[str, Any]]] = {c: [] for c in checks}
        for rep, p in zip(reps, params, strict=True):
            mdir = work / "mutated" / name / f"rep{rep.rep}"
            mdir.mkdir(parents=True)
            scored = {"inputs": {"parameters": MUTATIONS[name].apply(p)}}
            (mdir / "report.yaml").write_text(yaml.safe_dump(scored, sort_keys=False))
            argv = validate_argv(rep, layout, checks, "power_inheritance")
            argv[argv.index("--report") + 1] = str(mdir / "report.yaml")
            argv[argv.index("--output") + 1] = str(mdir / "validation.yaml")
            rep_cli([*argv, "--quiet"])
            data = load_validation(mdir / "validation.yaml")
            if data is None:
                raise RuntimeError(f"{scenario} rep{rep.rep} {name}: no validation.yaml")
            for c in checks:
                sections[c].append(data["checks"][c])
        decided = to_native(decide(scenario_family(sections, ON_DEMAND_CHECKS)))
        if name == "null":
            cli = _cli_verdict(scenario, roots, results)
            if {k: cli[k] for k in decided} != decided:
                raise RuntimeError(f"{scenario}: the in-process null verdict differs from simace validate's")
        out[name] = {"verdict": decided["verdict"], "n_tests": decided["n_tests"], "failed": failed_records(decided)}
    return out


def parse_span(text: str) -> list[int]:
    """``1-50`` or ``7`` as a list of study replicates."""
    lo, _, hi = text.partition("-")
    return list(range(int(lo), int(hi or lo) + 1))


def run_cell(cell: str, names: list[str], study_reps: list[int], parts: Path | None, jobs: int) -> dict[str, Any]:
    """Run study replicates of one cell; return the per-replicate outcomes.

    With ``parts``, the part file is rewritten after each study replicate, and a rerun resumes from it.
    """
    path = parts / f"{cell}~{study_reps[0]}-{study_reps[-1]}.json" if parts is not None else None
    part = {"cell": cell, "mutations": names, "study_reps": {}, "simace": simace_ref()}
    if path is not None and path.exists():
        part = json.loads(path.read_text())
        if part["mutations"] != names or part["simace"] != simace_ref():
            raise SystemExit(f"{path} holds other mutations or another simace; move it away to restart")
    for j in study_reps:
        if str(j) in part["study_reps"]:
            continue
        with tempfile.TemporaryDirectory(prefix=f"inhpow_{cell}_{j}_") as tmp:
            outcome = study_replicate(cell, j, names, Path(tmp), jobs)
        part["study_reps"][str(j)] = outcome
        line = ", ".join(f"{n} {o['verdict']}" for n, o in outcome.items() if o["verdict"] != "pass")
        print(f"{cell} study rep {j}: {line or 'every mutation passes'}", flush=True)
        if path is not None:
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(json.dumps(part))
    return part


def summarize(cell: str, name: str, outcomes: list[dict[str, Any]]) -> str:
    """One table row: verdict rejections, the target test's, by statistic, and the other tests'."""
    mutation = MUTATIONS[name]
    n = len(outcomes)
    k = sum(o["verdict"] == "fail" for o in outcomes)
    insufficient = sum(o["verdict"] == "insufficient" for o in outcomes)

    def test_of(record: str) -> str:
        return record.split("/")[0]

    target = mutation.test
    in_target = [[r for r in o["failed"] if test_of(r) == target] for o in outcomes] if target else []
    t = sum(bool(f) for f in in_target)
    stats = {s: sum(any(r.endswith(s) for r in f) for f in in_target) for s in ("mean_z", "sd_z", "reps_over_bound")}
    others = Counter(test for o in outcomes for test in {test_of(r) for r in o["failed"]} - {target})
    target_cells = f"{t}/{n} | {stats['mean_z']} / {stats['sd_z']} / {stats['reps_over_bound']}" if target else "– | –"
    other = ", ".join(f"{test} ({c})" for test, c in others.most_common())
    note = f" ({insufficient} insufficient)" if insufficient else ""
    return f"| {cell} | {name} | {mutation.effect} | {k}/{n}{note} | {interval(k, n)} | {target_cells} | {other} |"


def write_table(parts: list[dict[str, Any]], command: str) -> None:
    """Write power_inheritance.md from the parts: to results/ when every cell has the default study reps."""
    by_cell: dict[str, dict[str, dict[str, Any]]] = {}
    for p in parts:
        by_cell.setdefault(p["cell"], {}).update(p["study_reps"])
    full = set(by_cell) == set(CELLS) and all(
        sorted(map(int, reps)) == list(range(1, DEFAULT_STUDY_REPS + 1)) for reps in by_cell.values()
    )
    refs = sorted({p["simace"] for p in parts})
    md = [
        "# Power of simace validate's inheritance tests",
        "",
        f"Generated by `{command}` against {', '.join(refs)}; do not edit.",
        f"A study replicate is one scenario of {REPS_PER_VERDICT} reps; it rejects when simace's scenario verdict",
        "(`decide(scenario_family(...))`, every registered calibrated test at familywise 0.05) is `fail` with",
        "the reps' report.yaml parameters scored as the mutation says. The interval is the exact 95% interval",
        "of the rejection rate. *Target test* is the test that reads the mutated parameter; *by statistic*",
        "counts study replicates where it failed on mean z / SD z / reps over bound.",
        "",
        "Cells:",
        "",
        *(f"- `{c}`: {CELLS[c][0]}" for c in CELLS if c in by_cell),
        "",
        "| Cell | Mutation | Effect | Verdicts failing | 95% interval | Target test failing "
        "| By statistic (mean z / SD z / reps over bound) | Other tests failing |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for cell in (c for c in CELLS if c in by_cell):
        reps = by_cell[cell]
        names = [n for n in MUTATIONS if all(n in o for o in reps.values())]
        md.extend(summarize(cell, name, [reps[j][name] for j in sorted(reps, key=int)]) for name in names)
    out = RESULTS if full else PARTIAL_RESULTS
    out.mkdir(parents=True, exist_ok=True)
    (out / "power_inheritance.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"wrote {out / 'power_inheritance.md'}")


def main() -> None:
    """Run study replicates of the selected cells, or merge part files into power_inheritance.md."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--cell", default=",".join(CELLS), help="comma-separated cells")
    parser.add_argument("--mutation", help="comma-separated mutations (default: all)")
    parser.add_argument("--study-reps", default=f"1-{DEFAULT_STUDY_REPS}", help="study replicates, e.g. 1-50")
    parser.add_argument("--jobs", type=int, default=1, help="reps simace run and simace validate compute at once")
    parser.add_argument("--parts", type=Path, help="write one JSON part per cell here, not the table")
    parser.add_argument("--merge", type=Path, help="write the table from the part files in this directory")
    args = parser.parse_args()
    if args.merge:
        files = sorted(args.merge.glob("*.json"))
        if not files:
            raise SystemExit(f"no part files in {args.merge}")
        parts = [json.loads(f.read_text()) for f in files]
        write_table(parts, f"python -m inheritance.power --merge {args.merge.name}, from {len(parts)} part runs")
        return
    cells = args.cell.split(",")
    names = args.mutation.split(",") if args.mutation else list(MUTATIONS)
    if unknown := sorted(set(cells) - set(CELLS)):
        parser.error(f"unknown cells: {', '.join(unknown)}")
    if unknown := sorted(set(names) - set(MUTATIONS)):
        parser.error(f"unknown mutations: {', '.join(unknown)}")
    parts = [run_cell(cell, names, parse_span(args.study_reps), args.parts, args.jobs) for cell in cells]
    if not args.parts:
        write_table(parts, f"python -m inheritance.power --cell {args.cell} --study-reps {args.study_reps}")


if __name__ == "__main__":
    main()
