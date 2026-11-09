"""Power study: how often the harness rejects at stated oracle errors.

    pixi run power --case rm_base,ctm_base --reps 10

To fan out (for example one cluster job per case and mutation), give each
run --case, --mutation and --parts .work/power_parts, then assemble with
--merge .work/power_parts. Seeds depend only on the case and replicate, so
the merged table equals a single run's.

Each mutation changes the oracle's law, not simACE, so a rejection is what
the harness would report if simACE had that error. ``null`` checks the
oracle against its own draws: its rejection rate is the family-wise
false-positive rate of one case scored as its own Bonferroni family
(alpha/K over that case's p-values only), which Bonferroni holds at or below
0.05. A full harness run uses a larger K, so its per-case rate is lower.
Writes results/power.md when the table is every mutation of the default
cases at the default reps (one run, or a merge of all its parts); any other
selection is partial and goes to .work/results/power.md.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from scipy import stats

from am.manifest import XFT_MUTATIONS, case_spec
from run_harness import PARTIAL_RESULTS, RESULTS, assign_status, case_seed, load_cases, run_case
from structure.mating import MUTATIONS as MATING_MUTATIONS
from structure.mating import SLIM_MUTATIONS

# results/power.md holds every mutation of these cases at these reps; any other table is partial.
DEFAULT_CASES = "rm_base,ctm_base"
DEFAULT_REPS = 10

# Mutations by oracle; a case's oracle picks its table.
MUTATIONS = {
    "liability": {
        "null": "simACE replaced by one more oracle draw (false-positive rate)",
        "seg_045": "segregation variance 0.45*A instead of 0.5*A",
        "c_by_father": "household by father instead of mother",
        "rA_plus_01": "rA off by +0.1",
        "bm_plus_01": "maternal C transmission b_m off by +0.1 (both traits)",
    },
    "mating_law": MATING_MUTATIONS,
    "slim_standard": {name: effect for name, (effect, _) in SLIM_MUTATIONS.items()},
    # Layer 5 (ADR 0003, decision 8): the reference side is mutated; build each mutated
    # reference first with `pixi run -e xftsim xft-reference --case <id> --mutation <names>`.
    "xftsim": {
        "null": "simACE on fresh seeds against the reference (false-positive rate given that sample)",
        "null_xft": "a second, independent xftsim sample in simACE's place (xftsim against xftsim)",
        "r_090": "xftsim reference at mate correlation 0.9 r",
        "r_0": "xftsim reference at random mating",
        "mate_A": "xftsim reference mating on the additive genetic value instead of the phenotype",
        "mate_y1": "two traits: xftsim reference mating on trait 1 alone",
        "ctm_swap": "C transmission: xftsim reference with b_m and b_f swapped",
        "ctm_half": "C transmission: xftsim reference with b_m and b_f halved",
    },
}


def applicable(case: dict, name: str) -> bool:
    """Whether a mutation means anything for the case (a layer-5 mutation of r = 0 does not)."""
    if case["oracle"] == "xftsim" and name in XFT_MUTATIONS:
        return XFT_MUTATIONS[name][1](case_spec(case))
    return True


def interval(k: int, n: int) -> str:
    """Return the exact (Clopper-Pearson) 95% interval of k rejections in n runs."""
    ci = stats.binomtest(k, n).proportion_ci(confidence_level=0.95, method="exact")
    return f"{ci.low:.3f} to {ci.high:.3f}"


def run_mutation(case: dict, name: str, reps: int) -> dict:
    """Run one mutation of one case over independent harness replicates and summarise the rejections."""
    rejecting = 0
    n_failed: list[int] = []
    hits: dict[str, int] = {}
    for rep in range(1, reps + 1):
        c = case | {
            "id": f"{case['id']}~{name}~{rep}",
            # From the case id, so two cases never share simACE seeds (case_seed is 31-bit; R is small).
            "seed0": case_seed(case["id"], f"power{rep}"),
            "mutation": name,
        }
        rows = run_case(c)
        threshold = assign_status(rows)
        failed = [r["stat"] for r in rows if r["status"] == "FAIL"]
        rejecting += bool(failed)
        n_failed.append(len(failed))
        for s in failed:
            hits[s] = hits.get(s, 0) + 1
        print(f"{case['id']} {name} rep {rep}: {len(failed)} rows rejected (alpha/K {threshold:.3g})", flush=True)
    return {
        "case": case["id"],
        "oracle": case["oracle"],
        "mutation": name,
        "reps": reps,
        "rejecting": rejecting,
        "n_failed": n_failed,
        "hits": hits,
    }


def effect(part: dict) -> str:
    """Return a part's mutation description."""
    return MUTATIONS[part["oracle"]][part["mutation"]]


def is_full(parts: list[dict]) -> bool:
    """Return whether parts are the committed table's: every mutation of the default cases at the default reps."""
    oracles = {p["case"]: p["oracle"] for p in parts}
    want = sorted((c, m) for c in DEFAULT_CASES.split(",") if c in oracles for m in MUTATIONS[oracles[c]])
    return (
        set(oracles) == set(DEFAULT_CASES.split(","))
        and sorted((p["case"], p["mutation"]) for p in parts) == want
        and all(p["reps"] == DEFAULT_REPS for p in parts)
    )


def write_table(parts: list[dict], command: str, note: str = "") -> None:
    """Write power.md from one summary per (case, mutation): to results/ if is_full, else to .work/results/."""
    md = [
        "# Power",
        "",
        f"Generated by `{command}`{note}; do not edit.",
        "A run rejects when any row of the case fails at alpha/K. K counts that case's p-values only, so each",
        "run is one Bonferroni family of its case alone; a full harness run's K spans every case and its",
        "threshold is stricter. The interval is the exact 95% interval of the rejection rate.",
        "",
        "| Case | Mutation | Effect | Runs rejecting | 95% interval | Rows rejected per run (mean) | Most-rejected rows |",
        "|---|---|---|---|---|---|---|",
    ]
    for p in parts:
        k, n = p["rejecting"], p["reps"]
        top = ", ".join(f"{s} ({h})" for s, h in sorted(p["hits"].items(), key=lambda kv: -kv[1])[:5])
        md.append(
            f"| {p['case']} | {p['mutation']} | {effect(p)} | {k}/{n} | {interval(k, n)} "
            f"| {sum(p['n_failed']) / n:.1f} | {top} |"
        )
    out = RESULTS if is_full(parts) else PARTIAL_RESULTS
    out.mkdir(parents=True, exist_ok=True)
    Path(out / "power.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    print(f"wrote {out / 'power.md'}")


def merge(parts_dir: Path) -> None:
    """Assemble results/power.md from the part files of separate --parts runs."""
    files = sorted(parts_dir.glob("*.json"))
    if not files:
        raise SystemExit(f"no part files in {parts_dir}")
    parts = []
    for f in files:
        part = json.loads(f.read_text())
        if part.get("oracle") not in MUTATIONS or part.get("mutation") not in MUTATIONS[part["oracle"]]:
            raise SystemExit(f"{f}: unknown oracle or mutation {part.get('oracle')}/{part.get('mutation')}; rerun it")
        parts.append(part)
    parts.sort(key=lambda p: (p["case"], list(MUTATIONS[p["oracle"]]).index(p["mutation"])))
    write_table(parts, f"power.py --merge {parts_dir.name}", f", from {len(parts)} part runs")


def main() -> None:
    """Run the selected mutations on the selected cases, or merge part files into power.md."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--case", default=DEFAULT_CASES, help="comma-separated case ids")
    parser.add_argument("--mutation", help="comma-separated mutations (default: every mutation of each case's oracle)")
    parser.add_argument("--reps", type=int, default=DEFAULT_REPS)
    parser.add_argument("--parts", type=Path, help="write one JSON part per (case, mutation) here, not power.md")
    parser.add_argument("--merge", type=Path, help="write power.md from the part files in this directory")
    args = parser.parse_args()
    if args.merge:
        merge(args.merge)
        return
    cases = load_cases(parser, args.case)
    info = [c["id"] for c in cases if c.get("expect", "pass") == "info"]
    if info:
        parser.error(f"INFO cases are never scored, so they have no power: {', '.join(info)}")
    none = [f"{c['id']} ({c['oracle']})" for c in cases if c["oracle"] not in MUTATIONS]
    if none:
        with_some = sorted(c["id"] for c in load_cases(parser) if c.get("oracle") in MUTATIONS)
        parser.error(f"no mutations for {', '.join(none)}; cases with some: {', '.join(with_some)}")
    for case in cases:
        unknown = set(args.mutation.split(",") if args.mutation else ()) - set(MUTATIONS[case["oracle"]])
        if unknown:
            parser.error(f"{case['id']}: unknown mutations for {case['oracle']}: {', '.join(sorted(unknown))}")

    parts = []
    for case in cases:
        names = args.mutation.split(",") if args.mutation else list(MUTATIONS[case["oracle"]])
        names = [n for n in names if applicable(case, n)]
        for name in names:
            part = run_mutation(case, name, args.reps)
            if args.parts:
                args.parts.mkdir(parents=True, exist_ok=True)
                (args.parts / f"{case['id']}~{name}.json").write_text(json.dumps(part))
            parts.append(part)
    if not args.parts:
        mutations = f" --mutation {args.mutation}" if args.mutation else ""
        write_table(parts, f"power.py --case {args.case}{mutations} --reps {args.reps}")


if __name__ == "__main__":
    main()
