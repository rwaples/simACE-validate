"""Cross-check of the layer-2 engine: msprime against ped-sim on fixed families.

ped-sim (GPL-3.0, run as a binary, built by tools/build_pedsim.sh) and the
msprime drop in genedrop/drop.py each simulate C copies of three family
structures over the same genome, with Poisson crossovers and no
interference (ped-sim ``--pois``). The realised IBD fraction of each
relative pair, (IBD1 + 2 * IBD2 length) / (2 * genome length), is compared
between the two engines with a two-sample KS test and a Welch t-test.
This checks the gene-drop that layer 2 uses, not simACE.
"""

from __future__ import annotations

import subprocess
from collections import defaultdict
from typing import TYPE_CHECKING

import numpy as np
import polars as pl

from genedrop.drop import BP_PER_MORGAN, Genome, drop

if TYPE_CHECKING:
    from pathlib import Path

# name: (ped-sim def body, the printed pair's two sample suffixes, msprime template)
# A template row is (generation, mother, father) with ids local to one copy;
# the pair is the template's last two rows.
STRUCTURES = {
    "FS": ("2 2 1", ("_g2-b1-i1", "_g2-b1-i2"), [(0, -1, -1), (0, -1, -1), (1, 0, 1), (1, 0, 1)]),
    "HS": ("2 1 2 1:1 2:1", ("_g2-b1-i1", "_g2-b2-i1"), [(0, -1, -1), (0, -1, -1), (0, -1, -1), (1, 0, 1), (1, 0, 2)]),
    "1C": (
        "3 1",
        ("_g3-b1-i1", "_g3-b2-i1"),
        [(0, -1, -1), (0, -1, -1), (1, 0, 1), (1, 0, 1), (1, -1, -1), (1, -1, -1), (2, 2, 4), (2, 3, 5)],
    ),
}


def write_inputs(work: Path, genome: Genome, copies: int) -> tuple[Path, Path]:
    """Write the def file and a uniform sex-averaged map matching the msprime genome."""
    def_path = work / "families.def"
    def_path.write_text(
        "".join(f"def {name.lower()} {copies} {body.split()[0]}\n{body}\n" for name, (body, _, _) in STRUCTURES.items())
    )
    length_bp = int(genome.morgans * BP_PER_MORGAN)
    map_path = work / "uniform.map"
    map_path.write_text(
        "".join(f"{c} 1 0\n{c} {length_bp} {genome.morgans * 100}\n" for c in range(1, genome.chromosomes + 1))
    )
    return def_path, map_path


def pedsim_ibd(binary: Path, work: Path, genome: Genome, copies: int, seed: int) -> dict[str, np.ndarray]:
    """Run ped-sim and return each structure's per-pair IBD fraction."""
    work.mkdir(parents=True, exist_ok=True)
    def_path, map_path = write_inputs(work, genome, copies)
    out = work / "pedsim"
    cmd = [str(binary), "-d", str(def_path), "-m", str(map_path), "-o", str(out), "--pois", "--nogz"]
    subprocess.run([*cmd, "--seed", str(seed)], check=True, capture_output=True)
    shared: dict[tuple[str, str], float] = defaultdict(float)
    for line in out.with_suffix(".seg").read_text().splitlines():
        s1, s2, _, _, _, kind, _, _, cm = line.split()
        if kind in ("IBD1", "IBD2"):
            shared[tuple(sorted((s1, s2)))] += float(cm) * (2 if kind == "IBD2" else 1)
    total_cm = 2 * genome.morgans * 100 * genome.chromosomes
    result = {}
    for name, (_, (a, b), _) in STRUCTURES.items():
        pairs = [tuple(sorted((f"{name.lower()}{k}{a}", f"{name.lower()}{k}{b}"))) for k in range(1, copies + 1)]
        result[name] = np.array([shared[p] for p in pairs]) / total_cm
    return result


def msprime_ibd(genome: Genome, copies: int, seed: int) -> dict[str, np.ndarray]:
    """Drop the same structures through msprime and return each one's per-pair IBD fraction."""
    rows, pairs, labels = [], [], []
    for name, (_, _, template) in STRUCTURES.items():
        for _ in range(copies):
            base = len(rows)
            rows += [(g, -1 if m < 0 else base + m, -1 if f < 0 else base + f) for g, m, f in template]
            pairs.append((len(rows) - 2, len(rows) - 1))
            labels.append(name)
    gen, mother, father = (np.array(col, dtype=np.int32) for col in zip(*rows, strict=True))
    ped = pl.DataFrame({"generation": gen, "mother": mother, "father": father, "twin": np.full(len(gen), -1)})
    no_trait = Genome(genome.chromosomes, genome.morgans, causal_per_chromosome=0)
    _, ibd = drop(ped, np.array(pairs), no_trait, A=1.0, seed=seed)
    labels_arr = np.array(labels)
    return {name: ibd[labels_arr == name] for name in STRUCTURES}
