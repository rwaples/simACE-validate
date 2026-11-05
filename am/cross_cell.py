"""Are pilot cells at different r independent? Correlation over seeds of one statistic between r cells.

xftsim draws its E noise and mating noise from one numpy stream, and the
pilot ran every r with the same seed, so cells with the same seed may share
random numbers. Prints, per simulator and statistic, the correlation over
seeds of the per-seed value at r = 0 with r = 0.5 and r = 0.8, at each g.

Usage: python -m am.cross_cell PILOT_DIR
"""

import json
import sys
from pathlib import Path

import numpy as np

d = Path(sys.argv[1])
print("| simulator | stat | g | corr(r=0, r=0.5) | corr(r=0, r=0.8) |")
print("|---|---|---|---|---|")
for name, f in [("simACE", "simace.jsonl"), ("xftsim M=1000", "xft.jsonl")]:
    vals = {}
    for line in (d / f).read_text().splitlines():
        rec = json.loads(line)
        if rec.get("m", 1000) != 1000:
            continue
        for row in rec["rows"]:
            for s in ["r_FS_L", "r_FS_A", "var_L"]:
                if s in row:
                    vals[(rec["r"], row["generation"], s, rec["seed"])] = row[s]
    for s in ["r_FS_L", "r_FS_A", "var_L"]:
        for g in range(2, 7):
            seeds = sorted({k[3] for k in vals if k[1] == g and k[2] == s})
            x0 = np.array([vals[(0.0, g, s, sd)] for sd in seeds])
            c = [np.corrcoef(x0, [vals[(r, g, s, sd)] for sd in seeds])[0, 1] for r in (0.5, 0.8)]
            print(f"| {name} | {s} | {g} | {c[0]:+.2f} | {c[1]:+.2f} |")
