"""Compare layer-5 reference samples across a regeneration: records identical, manifests differ only where expected.

    python tools/compare_reference_records.py OLD_DIR NEW_DIR

For each OLD_DIR/*.jsonl, prints whether NEW_DIR's file holds the same per-seed records
(every line after the manifest, byte for byte) and which manifest keys changed, with the
spec keys that changed listed separately. Exits 1 if any shared file's records differ or
a file is missing from NEW_DIR.
"""

import json
import sys
from pathlib import Path


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
        same = a[1:] == b[1:]
        bad += not same
        print(
            f"{old.name}: records {'identical' if same else 'DIFFER'} ({len(a) - 1} vs {len(b) - 1}); "
            f"manifest keys changed: {', '.join(keys) or 'none'}; spec keys changed: {', '.join(spec) or 'none'}"
        )
    print(f"{bad} reference(s) differ or are missing")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
