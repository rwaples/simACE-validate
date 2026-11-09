"""Compare layer-5 reference samples across a regeneration: same draws, manifests differ only where expected.

    python tools/compare_reference_records.py OLD_DIR NEW_DIR

For each OLD_DIR/*.jsonl, prints whether NEW_DIR's file holds the same seeds and statistics,
the largest absolute difference of any statistic, and which manifest keys changed (spec keys
listed separately). Another host's CPU rounds reductions differently, so values from the same
draws can differ in the last bits; a difference above TOL means different draws. Exits 1 if any
shared file differs beyond TOL, or a file is missing from NEW_DIR.
"""

import json
import math
import sys
from pathlib import Path

TOL = 1e-12


def diff(a: float, b: float) -> float:
    """Absolute difference; 0 when both are NaN (half-sib rows of a regime without half-sibs), inf when one is."""
    if math.isnan(a) or math.isnan(b):
        return 0.0 if math.isnan(a) and math.isnan(b) else math.inf
    return abs(a - b)


def main() -> int:
    """Print one line per reference and a total."""
    old_dir, new_dir = Path(sys.argv[1]), Path(sys.argv[2])
    bad = 0
    for old in sorted(old_dir.glob("*.jsonl")):
        new = new_dir / old.name
        if not new.exists():
            print(f"{old.name}: MISSING in {new_dir}")
            bad += 1
            continue
        a, b = old.read_text().splitlines(), new.read_text().splitlines()
        ma, mb = json.loads(a[0])["manifest"], json.loads(b[0])["manifest"]
        keys = sorted(k for k in set(ma) | set(mb) if k != "spec" and ma.get(k) != mb.get(k))
        spec = sorted(k for k in set(ma["spec"]) | set(mb["spec"]) if ma["spec"].get(k) != mb["spec"].get(k))
        ra, rb = [json.loads(x) for x in a[1:]], [json.loads(x) for x in b[1:]]
        same_shape = [r["seed"] for r in ra] == [r["seed"] for r in rb] and all(
            set(x["stats"]) == set(y["stats"]) for x, y in zip(ra, rb, strict=True)
        )
        if same_shape:
            gap = max(
                (diff(x["stats"][k], y["stats"][k]) for x, y in zip(ra, rb, strict=True) for k in x["stats"]),
                default=0.0,
            )
        else:
            gap = float("inf")
        ok = gap <= TOL
        bad += not ok
        print(
            f"{old.name}: {'same draws' if ok else 'DIFFER'} (max |diff| {gap:.3g}, byte-identical {a[1:] == b[1:]}); "
            f"manifest keys changed: {', '.join(keys) or 'none'}; spec keys changed: {', '.join(spec) or 'none'}"
        )
    print(f"{bad} reference(s) differ beyond {TOL:g} or are missing")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
