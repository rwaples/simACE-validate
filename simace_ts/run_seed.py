"""Run one simACE-ts seed end to end (ADR 0002).

Stages, each a subprocess of this repo's pixi env (``sys.executable``, so the
``simace`` stages and ``params.yaml`` come from the same simace install):

1. ``simace simulate``
2. ``build_pedigree_tables.py`` with every generation a sample
3. ``genotype_drop_chrom.py`` per chromosome (drop and graft; the costly step)
4. per architecture: ``tstrait_assign_effects.py``, ``tstrait_gv_chrom.py``
   per chromosome, ``tstrait_augment_pedigree.py --per-generation``,
   ``simace cohort`` and ``simace analyze``
5. ``part.py``, which writes the seed's Part, ``part.json``

Every stage seed comes from its own namespace of the seed s,
``SeedSequence([s, stage, index])``, never from integer offsets of s.

A finished stage leaves a marker under ``.done/``, so a rerun (a job
resubmitted after its time limit) skips it. ``manifest.json`` records what the
seed's outputs depend on: the case config, the chromosomes, simace's version
and git ref, and the data pins. A rerun whose manifest differs is refused, so
a finished stage is never reused for a different configuration. A run on a
chromosome subset (a smoke run) writes to its own directory.

``timing.tsv`` has, per stage, its wall time, ``max_rss_mb`` (the largest
single-process peak of the stage, from ``os.wait4``; not the sum over its
processes) and ``cgroup_peak_mb`` (the peak of this run's cgroup so far, which
covers stages running at once; blank without cgroup v2). The grafted
``.trees`` are deleted after the last architecture unless ``--keep-trees``,
and so are their drop markers.

Usage::

    pixi run simace-ts --seed 990000 --data ~/data/humanity_sim --work .work/simace_ts
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import yaml
from score import CHROMOSOMES, OUTPUT_KEYS, installed_simace

if TYPE_CHECKING:
    from collections.abc import Callable

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "tskit_pipeline" / "scripts"
CASES = ROOT / "cases.yaml"
CASE_ID = "simace_ts_p2"

# Seed namespaces (ADR 0002 item 8). Append only: renumbering changes every seed.
STAGE_SIMULATE, STAGE_DROP, STAGE_EFFECTS, STAGE_COHORT, STAGE_ANALYZE = range(5)

THREAD_VARS = ("POLARS_MAX_THREADS", "NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")


def stage_seed(seed: int, stage: int, index: int = 0) -> int:
    """Return the seed of one stage (and chromosome or architecture index) of seed ``seed``, in [1, 2**31 - 1)."""
    return int(np.random.SeedSequence([seed, stage, index]).generate_state(1)[0] % (2**31 - 2)) + 1


def flags(values: dict) -> list[str]:
    """Turn ``{flag: value}`` into ``--flag value`` pairs."""
    return [part for k, v in values.items() for part in (f"--{k}", str(v))]


def seed_dir(work: Path, case: dict, seed: int, chromosomes: tuple[int, ...]) -> Path:
    """Return a seed's output directory; a chromosome subset gets its own."""
    name = f"seed{seed}" if chromosomes == CHROMOSOMES else f"seed{seed}_chr{'-'.join(map(str, chromosomes))}"
    return work / case["id"] / name


@dataclass(frozen=True)
class Data:
    """The pinned SimHumanity inputs under one directory (ADR 0002, Data)."""

    root: Path

    @property
    def canonical(self) -> Path:
        """The canonicalised per-chromosome founder trees."""
        return self.root / "preprocessed_p2" / "canonicalized"

    @property
    def catalog(self) -> Path:
        """The founder allele-frequency site catalog."""
        return self.root / "preprocessed_p2" / "site_catalog.parquet"

    @property
    def simhumanity(self) -> Path:
        """The SimHumanity checkout holding the recombination maps."""
        return self.root / "SimHumanity"

    def check(self, chromosomes: tuple[int, ...]) -> None:
        """Raise naming every input that is missing."""
        missing = [p for p in (self.catalog, self.simhumanity) if not p.exists()]
        missing += [p for c in chromosomes if not (p := self.canonical / f"chromosome_{c}.trees").exists()]
        if missing:
            raise FileNotFoundError("missing SimHumanity inputs: " + ", ".join(map(str, missing)))

    def pins(self) -> dict[str, str | None]:
        """Return the data pins this directory holds (ADR 0002, Data); None where a pin file is absent."""
        sha = self.root / "simhumanity_trees_RO.zip.sha256"
        fingerprint = self.root / "preprocessed_p2" / "all_chroms.trees.fingerprint"
        head = subprocess.run(
            ["git", "-C", str(self.simhumanity), "rev-parse", "HEAD"], capture_output=True, text=True, check=True
        )
        return {
            "trees_zip_sha256": sha.read_text().split()[0] if sha.exists() else None,
            "simhumanity_commit": head.stdout.strip(),
            "preprocessed_fingerprint": fingerprint.read_text().strip() if fingerprint.exists() else None,
        }


def validate_commit() -> str:
    """Return this checkout's commit, with ``-dirty`` when its tracked files differ."""
    return subprocess.run(
        ["git", "-C", str(ROOT), "describe", "--always", "--dirty", "--abbrev=12"],
        capture_output=True, text=True, check=True,
    ).stdout.strip()  # fmt: skip


def load_case(case_id: str) -> dict:
    """Return the case of id ``case_id`` from cases.yaml."""
    cases = {c["id"]: c for c in yaml.safe_load(CASES.read_text())["cases"]}
    if case_id not in cases or cases[case_id]["oracle"] != "simace_ts":
        raise SystemExit(f"{case_id}: not a simace_ts case in {CASES}")
    return cases[case_id]


def check_manifest(out: Path, manifest: dict) -> None:
    """Write ``manifest.json``, or refuse a rerun whose inputs differ from the first run's.

    Of the case, only the keys the outputs depend on (``OUTPUT_KEYS``) are compared, so editing R or a note
    does not invalidate a seed.
    """
    path = out / "manifest.json"
    if path.exists():
        previous = json.loads(path.read_text())
        changed = sorted(k for k in manifest if k not in ("case", "runs") and previous.get(k) != manifest[k])
        changed += [f"case.{k}" for k in OUTPUT_KEYS if previous["case"].get(k) != manifest["case"].get(k)]
        if changed:
            raise SystemExit(f"{out}: {', '.join(changed)} changed since this seed's first run; use a new --work")
        manifest["runs"] = previous["runs"] + manifest["runs"]
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(manifest, indent=2) + "\n")
    tmp.replace(path)  # atomic, so a killed job never leaves a truncated manifest


def cgroup_peak_mb() -> str:
    """Return this process's cgroup v2 ``memory.peak`` in MB, or '' when it is not available."""
    try:
        rel = Path("/proc/self/cgroup").read_text().split("::", 1)[1].strip()
        return f"{int((Path('/sys/fs/cgroup') / rel.lstrip('/') / 'memory.peak').read_text()) / 2**20:.0f}"
    except OSError, IndexError, ValueError:
        return ""


class Stages:
    """Run stages once each, recording wall time and memory."""

    def __init__(self, out: Path) -> None:
        """Record stages under ``out``."""
        self.out = out
        self.done = out / ".done"
        self.done.mkdir(parents=True, exist_ok=True)
        self.timing = out / "timing.tsv"
        self.lock = threading.Lock()
        if not self.timing.exists():
            self.timing.write_text("stage\twall_s\tmax_rss_mb\tcgroup_peak_mb\n")

    def marker(self, name: str) -> Path:
        """Return the done marker of stage ``name``."""
        return self.done / name

    def run(self, name: str, cmd: list[str], threads: int = 1) -> None:
        """Run ``cmd`` with ``threads`` threads unless ``name`` is already done; raise on a non-zero exit."""
        marker = self.marker(name)
        if marker.exists():
            return
        log = self.out / "logs" / f"{name}.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        env = {**os.environ, **dict.fromkeys(THREAD_VARS, str(threads))}
        start = time.perf_counter()
        with log.open("w") as fh:
            fh.write(shlex.join(cmd) + "\n")
            fh.flush()
            proc = subprocess.Popen(cmd, stdout=fh, stderr=subprocess.STDOUT, env=env)
            _, status, usage = os.wait4(proc.pid, 0)
            proc.returncode = os.waitstatus_to_exitcode(status)  # reaped here, so Popen must not wait on the pid
        wall = time.perf_counter() - start
        if proc.returncode != 0:
            raise RuntimeError(f"stage {name} exited {proc.returncode}; see {log}")
        with self.lock, self.timing.open("a") as fh:
            fh.write(f"{name}\t{wall:.1f}\t{usage.ru_maxrss / 1024:.0f}\t{cgroup_peak_mb()}\n")
        marker.parent.mkdir(parents=True, exist_ok=True)
        marker.touch()


def run_seed(
    case: dict, seed: int, data: Data, work: Path, jobs: int, chromosomes: tuple[int, ...], keep_trees: bool
) -> Path:
    """Run every stage of one seed and return its output directory."""
    data.check(chromosomes)
    out = seed_dir(work, case, seed, chromosomes)
    out.mkdir(parents=True, exist_ok=True)
    manifest = {
        "case": case,
        "seed": seed,
        "chromosomes": list(chromosomes),
        "simace": installed_simace(),
        "data": data.pins(),
        "runs": [{"simace_validate": validate_commit(), "started": time.strftime("%Y-%m-%dT%H:%M:%S")}],
    }
    if manifest["data"] != case["data"]:
        raise SystemExit(f"{data.root}: data pins {manifest['data']} are not the case's {case['data']}")
    check_manifest(out, manifest)
    stages = Stages(out)
    py = sys.executable
    g_ped = int(case["simulate"]["G-ped"])
    ped = out / "pedigree.parquet"
    tables = out / "pedigree_tables.trees"

    stages.run(
        "simulate",
        [py, "-m", "simace", "simulate", "--seed", str(stage_seed(seed, STAGE_SIMULATE)), *flags(case["simulate"]),
         "--output-pedigree", str(ped)],
        threads=jobs,
    )  # fmt: skip
    stages.run(
        "build_tables",
        [py, str(SCRIPTS / "build_pedigree_tables.py"), "--pedigree", str(ped), "--tables", str(tables),
         "--g-ped", str(g_ped), "--g-pheno", str(g_ped)],
    )  # fmt: skip

    def trees(c: int) -> Path:
        return out / f"genotypes_chrom_{c}.trees"

    def drop(c: int) -> None:
        stages.run(
            f"drop_chrom{c}",
            [py, str(SCRIPTS / "genotype_drop_chrom.py"), "--pedigree-tables", str(tables),
             "--trees", str(data.canonical / f"chromosome_{c}.trees"), "--out-trees", str(trees(c)),
             "--chrom", str(c), "--g-ped", str(g_ped), "--seed", str(stage_seed(seed, STAGE_DROP, c)),
             "--simhumanity-dir", str(data.simhumanity)],
        )  # fmt: skip

    # The trees are deleted once every architecture is done; drop again only if a stage still needs them.
    needs_trees = any(
        not stages.marker(f"{a['name']}/gv_chrom{c}").exists() for a in case["architectures"] for c in chromosomes
    )
    if needs_trees:
        with ThreadPoolExecutor(jobs) as pool:
            list(pool.map(drop, chromosomes))

    params = write_params(case, seed, out / "params.yaml")
    for i, arch in enumerate(case["architectures"]):
        run_architecture(stages, case, seed, i, arch, out, ped, params, trees, chromosomes, jobs, data)

    # The Part is cheap and derived only from outputs on disk, so it is recomputed on every run.
    stages.marker("part").unlink(missing_ok=True)
    stages.run("part", [py, str(Path(__file__).resolve().parent / "part.py"), str(out)])

    if not keep_trees:
        for c in chromosomes:
            trees(c).unlink(missing_ok=True)
            stages.marker(f"drop_chrom{c}").unlink(missing_ok=True)  # the drop must rerun if trees are needed again
    return out


def run_architecture(
    stages: Stages, case: dict, seed: int, i: int, arch: dict, out: Path, ped: Path, params: Path,
    trees: Callable[[int], Path], chromosomes: tuple[int, ...], jobs: int, data: Data,
) -> None:  # fmt: skip
    """Assign one architecture's effects, compute its genetic values, and run simACE's cohort and analyze on it."""
    name = arch["name"]
    d = out / name
    d.mkdir(exist_ok=True)
    py = sys.executable
    effects = d / "causal_effects.parquet"
    arch_flags = {k: v for k, v in arch.items() if k != "name"}
    stages.run(
        f"{name}/effects",
        [py, str(SCRIPTS / "tstrait_assign_effects.py"), "--catalog", str(data.catalog), "--effects", str(effects),
         "--meta", str(d / "causal_effects_meta.json"), *flags(arch_flags), "--trait-id", "0",
         "--seed", str(stage_seed(seed, STAGE_EFFECTS, i))],
        threads=jobs,
    )  # fmt: skip

    def gv(c: int) -> None:
        stages.run(
            f"{name}/gv_chrom{c}",
            [py, str(SCRIPTS / "tstrait_gv_chrom.py"), "--trees", str(trees(c)), "--effects", str(effects),
             "--gv", str(d / f"gv_chrom_{c}.parquet"), "--meta", str(d / f"gv_chrom_{c}_meta.json"),
             "--chrom", str(c), "--trait-id", "0"],
        )  # fmt: skip

    with ThreadPoolExecutor(jobs) as pool:
        list(pool.map(gv, chromosomes))

    augmented = d / "pedigree.full.tstrait.parquet"
    stages.run(
        f"{name}/augment",
        [py, str(SCRIPTS / "tstrait_augment_pedigree.py"), "--pedigree", str(ped),
         "--gv", *[str(d / f"gv_chrom_{c}.parquet") for c in chromosomes], "--out-pedigree", str(augmented),
         "--meta", str(d / "augment_meta.json"), "--target-var-a", str(case["simulate"]["A1"]), "--per-generation"],
        threads=jobs,
    )  # fmt: skip
    cohort, phenotyped = d / "cohort.parquet", d / "phenotyped_population.yaml"
    stages.run(
        f"{name}/cohort",
        [py, "-m", "simace", "cohort", "--pedigree", str(augmented), "--output-cohort", str(cohort),
         "--output-phenotyped-population", str(phenotyped), "--seed", str(stage_seed(seed, STAGE_COHORT, i)),
         *flags(case["cohort"])],
        threads=jobs,
    )  # fmt: skip
    analyze_flags = {"censor-age": case["cohort"]["censor-age"]}  # one value for both stages, as simace run does
    stages.run(
        f"{name}/analyze",
        [py, "-m", "simace", "analyze", "--pedigree", str(augmented), "--params", str(params), "--cohort", str(cohort),
         "--phenotyped-population", str(phenotyped), "--folder", "simace_ts", "--scenario", case["id"],
         "--rep", "1", "--report-output", str(d / "report.yaml"), "--plot-payload-output", str(d / "plot_payload.yaml"),
         "--samples-output", str(d / "plotting_sample.parquet"), "--seed", str(stage_seed(seed, STAGE_ANALYZE, i)),
         *flags(analyze_flags)],
        threads=jobs,
    )  # fmt: skip


def write_params(case: dict, seed: int, path: Path) -> Path:
    """Write the params.yaml ``simace analyze`` reads, as ``simace run`` would (simace ``emit_params``)."""
    from simace.simulation.emit_params import emit_params

    sim = {k.replace("-", "_"): v for k, v in case["simulate"].items()}
    params = emit_params(
        seed=stage_seed(seed, STAGE_SIMULATE), rep=1, G_pheno=int(case["cohort"]["G-pheno"]),
        assort1=0.0, assort2=0.0, **sim,
    )  # fmt: skip
    path.write_text(yaml.safe_dump(params, sort_keys=False))
    return path


def main(argv: list[str] | None = None) -> int:
    """Run one seed of the simACE-ts case."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--seed", type=int, required=True, help="a scored seed (seed0 ..) or a pilot seed")
    parser.add_argument("--data", type=Path, required=True, help="directory holding preprocessed_p2/ and SimHumanity/")
    parser.add_argument("--work", type=Path, default=ROOT / ".work" / "simace_ts")
    parser.add_argument("--jobs", type=int, default=int(os.environ.get("SLURM_CPUS_PER_TASK", "1")))
    parser.add_argument(
        "--chroms", type=int, nargs="+", default=list(CHROMOSOMES), help="a subset, for smoke runs only"
    )
    parser.add_argument("--keep-trees", action="store_true")
    parser.add_argument("--case", default=CASE_ID, help="a simace_ts case id in cases.yaml")
    args = parser.parse_args(argv)

    case = load_case(args.case)
    chromosomes = tuple(sorted(set(args.chroms)))
    out = run_seed(case, args.seed, Data(args.data.expanduser()), args.work, args.jobs, chromosomes, args.keep_trees)
    print(out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
