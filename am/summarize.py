"""Summarize the issue #4 pilot: simACE vs xftsim under random and assortative mating.

Usage: summarize.py PILOT_DIR   (reads simace.jsonl and xft.jsonl, prints Markdown)

Every number in the design note comes from this output. For each statistic,
generation and r it prints:
- the mean and SE over seeds for simACE and for xftsim (M = 1000);
- the Bulmer/infinitesimal prediction (am/stats.py::bulmer), N -> infinity;
- a two-sample Welch t and KS p for simACE vs xftsim;
- the minimum mean gap a Welch t over R seeds per side would detect with
  80% power at a harness-sized alpha/K.
It also prints the finite-M gap: xftsim Var(A) minus the prediction, by M.
"""

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy import stats as st

sys.path.insert(0, str(Path(__file__).resolve().parent))
from stats import bulmer, bulmer2

S = 0.5  # C + E of the pilot trait
V0 = 0.5
STATS = ["var_A", "var_L", "r_mate_L", "r_mate_A", "r_FS_A", "r_FS_L", "r_PO_A", "r_PO_L", "r_MO_L", "r_FO_L"]
ALPHA_K = 0.05 / 1000  # a harness-sized Bonferroni threshold (status.md: K = 835)


def predict(r: float, g: int) -> dict[str, float]:
    """Bulmer-recursion expectation of every statistic after g matings at mate correlation r."""
    v = bulmer(V0, S, r, g)
    vt, vp = v[g], v[g - 1]
    lt, lp = vt + S, vp + S
    mid = vp / 2 * (1 + r * vp / lp)  # Var(midparent A) = Cov(A) of full sibs
    return {
        "var_A": vt,
        "var_L": lt,
        "r_mate_L": r,
        "r_mate_A": r * vp / lp,
        "r_FS_A": mid / vt,
        "r_FS_L": mid / lt,
        "r_PO_A": mid / np.sqrt(vp * vt),
        "r_PO_L": vp * (1 + r) / 2 / np.sqrt(lp * lt),
        "r_MO_L": vp * (1 + r) / 2 / np.sqrt(lp * lt),
        "r_FO_L": vp * (1 + r) / 2 / np.sqrt(lp * lt),
    }


def load(path: Path) -> list[dict]:
    """Read one JSON record per line."""
    return [json.loads(line) for line in path.read_text().splitlines() if line.strip()]


def table(recs: list[dict], traits: int = 1, stats: list[str] = STATS) -> dict:
    """Map (r, m, generation, stat) to the array of that statistic over seeds, for one trait mode."""
    out = defaultdict(list)
    for rec in recs:
        if rec.get("traits", 1) != traits:
            continue
        for row in rec["rows"]:
            for k, v in row.items():
                if k in stats:
                    out[(rec["r"], rec.get("m"), row["generation"], k)].append(v)
    return {k: np.array(v) for k, v in out.items()}


def mdd(sd1: float, sd2: float, n: int, alpha: float) -> float:
    """Mean gap a two-sided Welch t over n seeds per side detects with 80% power."""
    se = np.sqrt((sd1**2 + sd2**2) / n)
    return float((st.norm.ppf(1 - alpha / 2) + st.norm.ppf(0.8)) * se)


def main() -> None:
    """Print the pilot summary for the directory in argv[1]."""
    d = Path(sys.argv[1])
    sa = table(load(d / "simace.jsonl"))
    xf = table(load(d / "xft.jsonl"))
    rs = sorted({k[0] for k in sa})
    gens = sorted({k[2] for k in sa})
    n_seeds = {len(v) for v in sa.values()} | {len(v) for v in xf.values()}
    print(f"Seeds per cell: {sorted(n_seeds)}. alpha/K for MDD: {ALPHA_K:.2e}.\n")

    for r in rs:
        print(f"## r = {r}\n")
        print(
            "| stat | g | predicted | simACE mean (SE) | xftsim M=1000 mean (SE) | sim−xft | Welch p | KS p "
            "| MDD R=20 | MDD R=50 | MDD R=100 |"
        )
        print("|---|---|---|---|---|---|---|---|---|---|---|")
        for s in STATS:
            for g in gens:
                a = sa.get((r, None, g, s))
                b = xf.get((r, 1000, g, s))
                if a is None or b is None or len(a) < 2:
                    continue
                p = predict(r, g)[s]
                tp = st.ttest_ind(a, b, equal_var=False).pvalue
                kp = st.ks_2samp(a, b).pvalue
                sd1, sd2 = a.std(ddof=1), b.std(ddof=1)
                print(
                    f"| {s} | {g} | {p:.4f} | {a.mean():.4f} ({sd1 / np.sqrt(len(a)):.4f}) "
                    f"| {b.mean():.4f} ({sd2 / np.sqrt(len(b)):.4f}) | {a.mean() - b.mean():+.4f} "
                    f"| {tp:.3g} | {kp:.3g} | {mdd(sd1, sd2, 20, ALPHA_K):.4f} "
                    f"| {mdd(sd1, sd2, 50, ALPHA_K):.4f} | {mdd(sd1, sd2, 100, ALPHA_K):.4f} |"
                )
        print()

    print("## Gap from the prediction: (mean − predicted) / SE, by simulator\n")
    print("| r | stat | g | simACE z | xftsim M=1000 z |")
    print("|---|---|---|---|---|")
    for r in rs:
        for s in ["var_A", "r_mate_L", "r_FS_L", "r_PO_L"]:
            for g in gens:
                a = sa.get((r, None, g, s))
                b = xf.get((r, 1000, g, s))
                if a is None or b is None:
                    continue
                p = predict(r, g)[s]
                za = (a.mean() - p) / (a.std(ddof=1) / np.sqrt(len(a)))
                zb = (b.mean() - p) / (b.std(ddof=1) / np.sqrt(len(b)))
                print(f"| {r} | {s} | {g} | {za:+.2f} | {zb:+.2f} |")
    print()

    print("## Sensitivity: prediction at each simulator's realised mean r_mate_L (g >= 2)\n")
    print(
        "| r | simACE r_mate_L | xftsim r_mate_L | stat | g | pred at r | pred at simACE r | pred at xftsim r "
        "| simACE mean | xftsim mean |"
    )
    print("|---|---|---|---|---|---|---|---|---|---|")
    for r in rs:
        if r == 0:
            continue
        ra = float(np.mean([sa[(r, None, g, "r_mate_L")].mean() for g in gens if (r, None, g, "r_mate_L") in sa]))
        rb = float(np.mean([xf[(r, 1000, g, "r_mate_L")].mean() for g in gens if (r, 1000, g, "r_mate_L") in xf]))
        g = max(gens)
        for s in ["var_A", "r_FS_A", "r_PO_A"]:
            print(
                f"| {r} | {ra:.4f} | {rb:.4f} | {s} | {g} | {predict(r, g)[s]:.4f} | {predict(ra, g)[s]:.4f} "
                f"| {predict(rb, g)[s]:.4f} | {sa[(r, None, g, s)].mean():.4f} | {xf[(r, 1000, g, s)].mean():.4f} |"
            )
    print()

    print("## Spread: SD over seeds, simACE / xftsim M=1000 (median and range over r and g)\n")
    print("| stat | median | min | max |")
    print("|---|---|---|---|")
    for s in STATS:
        ratios = [
            sa[(r, None, g, s)].std(ddof=1) / xf[(r, 1000, g, s)].std(ddof=1)
            for r in rs
            for g in gens
            if (r, None, g, s) in sa and (r, 1000, g, s) in xf and len(sa[(r, None, g, s)]) > 1
        ]
        print(f"| {s} | {np.median(ratios):.2f} | {min(ratios):.2f} | {max(ratios):.2f} |")
    print()

    n_tests = sum(
        1
        for r in rs
        for s in STATS
        for g in gens
        if (r, None, g, s) in sa and (r, 1000, g, s) in xf and len(sa[(r, None, g, s)]) > 1
    )
    print(
        f"Two-sample rows above: {n_tests} (each with a Welch and a KS p); "
        f"Bonferroni at 0.05 over {2 * n_tests} p-values: {0.05 / (2 * n_tests):.2e}.\n"
    )

    rows = []
    for r in rs:
        for s in STATS:
            for g in gens:
                a, b = sa.get((r, None, g, s)), xf.get((r, 1000, g, s))
                if a is None or b is None or len(a) < 2:
                    continue
                sd1, sd2 = a.std(ddof=1), b.std(ddof=1)
                rows.append(
                    (
                        st.ttest_ind(a, b, equal_var=False).pvalue,
                        st.ks_2samp(a, b).pvalue,
                        r,
                        s,
                        g,
                        mdd(sd1, sd2, 50, ALPHA_K),
                        mdd(sd1, sd2, 100, ALPHA_K),
                    )
                )
    best = min(rows)
    print(
        f"Smallest Welch p: {best[0]:.3g} ({best[3]}, r = {best[2]}, g = {best[4]}; KS p {best[1]:.3g}). "
        f"Smallest KS p: {min(x[1] for x in rows):.3g}.\n"
    )
    print("## MDD range by row group (over r and g)\n")
    print("| group | R = 50 | R = 100 |")
    print("|---|---|---|")
    groups = {
        "Var(A)": ["var_A"],
        "Var(L)": ["var_L"],
        "r_mate_L": ["r_mate_L"],
        "r_mate_A": ["r_mate_A"],
        "FS (A, L)": ["r_FS_A", "r_FS_L"],
        "PO (A, L)": ["r_PO_A", "r_PO_L"],
    }
    for name, ss in groups.items():
        sel = [x for x in rows if x[3] in ss]
        print(
            f"| {name} | {min(x[5] for x in sel):.4f} to {max(x[5] for x in sel):.4f} "
            f"| {min(x[6] for x in sel):.4f} to {max(x[6] for x in sel):.4f} |"
        )
    print()

    ms = sorted({k[1] for k in xf if k[1] is not None})

    def ci(x: np.ndarray, p: float) -> tuple[float, float, float]:
        m, se = x.mean() - p, x.std(ddof=1) / np.sqrt(len(x))
        t = st.t.ppf(0.975, len(x) - 1)
        return m, m - t * se, m + t * se

    for r in [r for r in rs if any((r, m, 1, "var_A") in xf for m in ms if m != 1000)]:
        print(f"## Finite M: Var(A) − prediction at r = {r}, mean [95% t CI]\n")
        print("| g | predicted | " + " | ".join(f"xftsim M={m}" for m in ms) + " | simACE |")
        print("|---|---|" + "---|" * (len(ms) + 1))
        worst = defaultdict(float)
        for g in [0, *gens]:
            p = predict(r, g)["var_A"] if g > 0 else V0
            cells = []
            for m in ms:
                b = xf.get((r, m, g, "var_A"))
                if b is None:
                    cells.append("")
                    continue
                c = ci(b, p)
                worst[m] = max(worst[m], abs(c[1]), abs(c[2]))
                cells.append(f"{c[0]:+.4f} [{c[1]:+.4f}, {c[2]:+.4f}]")
            a = sa.get((r, None, g, "var_A"))
            if a is None:
                cells.append("")
            else:
                c = ci(a, p)
                worst["simACE"] = max(worst["simACE"], abs(c[1]), abs(c[2]))
                cells.append(f"{c[0]:+.4f} [{c[1]:+.4f}, {c[2]:+.4f}]")
            print(f"| {g} | {p:.4f} | " + " | ".join(cells) + " |")
        print()
        print(
            "Largest |CI endpoint| over g (a bound on the gap from the recursion that these seeds support): "
            + ", ".join(f"{k}: {v:.4f}" for k, v in worst.items())
            + "\n"
        )
        lo, hi = min(ms), max(ms)
        print(f"Sentinel, xftsim M={lo} − M={hi} Var(A), mean [95% Welch CI] by g: ", end="")
        parts = []
        for g in gens:
            x1, x2 = xf[(r, lo, g, "var_A")], xf[(r, hi, g, "var_A")]
            res = st.ttest_ind(x1, x2, equal_var=False)
            c = res.confidence_interval(0.95)
            parts.append(f"g{g} {x1.mean() - x2.mean():+.4f} [{c.low:+.4f}, {c.high:+.4f}]")
        print("; ".join(parts) + "\n")

    print("## Cost per run (median wall s, max RSS MB)\n")
    for name, f in [("simace", "simace.jsonl"), ("xftsim", "xft.jsonl")]:
        recs = load(d / f)
        by = defaultdict(list)
        for rec in recs:
            by[(rec.get("traits", 1), rec.get("m") or 0)].append(rec)
        for (t, m), rr in sorted(by.items()):
            w = np.median([x["wall_s"] for x in rr])
            rss = max((x.get("maxrss_mb", float("nan")) for x in rr), default=float("nan"))
            print(f"- {name} traits={t} M={m or '-'}: {len(rr)} runs, median wall {w:.1f} s, max RSS {rss:.0f} MB")
    print()
    two_trait(d)


STATS2 = [
    "var_A", "var_A_t2", "var_L", "r_w_A", "r_w_L", "r_w_L_parents",
    "r_mate_L", "r_mate_L_t2", "r_mate_L12", "r_mate_L21",
    "r_mate_A", "r_mate_A12", "r_FS_A", "r_FS_L", "r_PO_A", "r_PO_L",
]  # fmt: skip


def two_trait(d: Path) -> None:
    """Two-trait AM: every mate correlation r; simACE vs xftsim and the two-trait Bulmer recursion."""
    sa = table(load(d / "simace.jsonl"), traits=2, stats=STATS2)
    xf = table(load(d / "xft.jsonl"), traits=2, stats=STATS2)
    if not sa:
        return
    rs = sorted({k[0] for k in sa})
    gens = sorted({k[2] for k in sa})
    rows = []
    for r in rs:
        pred = bulmer2(V0, S, r, max(gens))
        print(f"## Two-trait AM, every mate correlation r = {r}\n")
        print("| stat | g | predicted | simACE mean (SE) | xftsim mean (SE) | sim−xft | Welch p | z sim | z xft |")
        print("|---|---|---|---|---|---|---|---|---|")
        for s_ in STATS2:
            for g in gens:
                a, b = sa.get((r, None, g, s_)), xf.get((r, 1000, g, s_))
                if a is None or b is None or len(a) < 2:
                    continue
                p = pred[g].get(s_, r if s_.startswith("r_mate_L") else None)
                se_a, se_b = a.std(ddof=1) / np.sqrt(len(a)), b.std(ddof=1) / np.sqrt(len(b))
                tp = st.ttest_ind(a, b, equal_var=False).pvalue
                rows.append((tp, r, s_, g, mdd(a.std(ddof=1), b.std(ddof=1), 100, ALPHA_K)))
                zs = "" if p is None else f"{(a.mean() - p) / se_a:+.2f}"
                zx = "" if p is None else f"{(b.mean() - p) / se_b:+.2f}"
                ps = "" if p is None else f"{p:.4f}"
                print(
                    f"| {s_} | {g} | {ps} | {a.mean():.4f} ({se_a:.4f}) | {b.mean():.4f} ({se_b:.4f}) "
                    f"| {a.mean() - b.mean():+.4f} | {tp:.3g} | {zs} | {zx} |"
                )
        print()
    best = min(rows)
    print(
        f"Two-trait rows: {len(rows)}; Bonferroni at 0.05: {0.05 / len(rows):.2e}. "
        f"Smallest Welch p: {best[0]:.3g} ({best[2]}, r = {best[1]}, g = {best[3]}). "
        f"MDD at R = 100: {min(x[4] for x in rows):.4f} to {max(x[4] for x in rows):.4f}.\n"
    )


if __name__ == "__main__":
    main()
