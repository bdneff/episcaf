#!/usr/bin/env python3
"""
plot_fp_reduction.py

Clearer views of the native-volume carve's effect than the twin scatter.

  --mode dist        For designs the antibody tolerates (true clash == 0), the cylinder
                     count distribution before vs after. The high-count tail (the false
                     positives) collapses leftward. One population, one axis, obvious shift.

  --mode selectivity Across true-clash bins, the fraction of designs the cylinder flags
                     heavily (cylinder >= thr), before vs after. The gap is large where
                     there's no real clash and vanishes where the clash is real -- so you
                     can see it's removing false positives, not signal, in one panel.

    python scripts/plot_fp_reduction.py --native metrics_native_cyl.csv --mode dist
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BEFORE = "#9aa7b4"   # muted grey-blue (original)
AFTER = "#c0504d"    # warm red (native-aware)


def load(path):
    d = pd.read_csv(path, low_memory=False)
    tc = pd.to_numeric(d["af3_n_clash_res"], errors="coerce")
    before = pd.to_numeric(d["cylinder_ca_clashes"], errors="coerce")
    after = pd.to_numeric(d["cylinder_native_aware"], errors="coerce")
    m = tc.notna() & before.notna() & after.notna()   # gated (carved) designs
    return tc[m].to_numpy(), before[m].to_numpy(), after[m].to_numpy()


def plot_dist(tc, before, after, out, thr):
    b, a = before[tc == 0], after[tc == 0]
    n = len(b)
    hi = int(max(b.max(), a.max()))
    bins = np.arange(0, hi + 3, 2)
    fig, ax = plt.subplots(figsize=(7.5, 4.6))
    ax.hist(b, bins=bins, color=BEFORE, alpha=0.85, label="original cylinder")
    ax.hist(a, bins=bins, color=AFTER, alpha=0.65,
            label="native antigen volume subtracted")
    nb, na = int((b >= thr).sum()), int((a >= thr).sum())
    ax.axvline(thr, color="#444", lw=1, ls="--")
    ax.annotate(f"heavily penalised (cylinder \u2265 {thr}):\n"
                f"{nb:,}  \u2192  {na:,} designs",
                xy=(thr, ax.get_ylim()[1] * 0.82), xytext=(thr + hi * 0.18,
                ax.get_ylim()[1] * 0.82), fontsize=10, color="#333",
                arrowprops=dict(arrowstyle="->", color="#333"))
    ax.set_xlabel("cylinder count (scaffold residues in the antibody approach path)")
    ax.set_ylabel("number of designs")
    ax.set_title(f"Designs the antibody tolerates (zero AlphaFold clash, n={n:,})\n"
                 "cylinder count collapses after subtracting native antigen volume",
                 fontsize=11)
    ax.legend(frameon=False, fontsize=9.5)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=150, bbox_inches="tight")
    print(f"wrote {out}.png/.pdf (dist; n={n}, >= {thr}: {nb}->{na})")


def plot_selectivity(tc, before, after, out, thr):
    edges = [0, 1, 4, 7, 11, 16, 26, 10**9]
    labels = ["0", "1\u20133", "4\u20136", "7\u201310", "11\u201315", "16\u201325", "26+"]
    fb, fa, xs = [], [], []
    for i in range(len(edges) - 1):
        m = (tc >= edges[i]) & (tc < edges[i + 1])
        if m.sum() < 20:
            fb.append(np.nan); fa.append(np.nan); xs.append(i); continue
        fb.append((before[m] >= thr).mean())
        fa.append((after[m] >= thr).mean())
        xs.append(i)
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(7.8, 4.6))
    ax.plot(x, fb, "-o", color=BEFORE, lw=2, ms=6, label="original cylinder")
    ax.plot(x, fa, "-o", color=AFTER, lw=2, ms=6,
            label="native antigen volume subtracted")
    ax.fill_between(x, fa, fb, color=AFTER, alpha=0.10)
    ax.set_xticks(x); ax.set_xticklabels(labels)
    ax.set_xlabel("AlphaFold antibody clashes (true clash)")
    ax.set_ylabel(f"fraction of designs flagged (cylinder \u2265 {thr})")
    ax.set_title("The carve removes flags where there's no real clash,\n"
                 "and leaves the genuinely-clashing designs alone", fontsize=11)
    ax.annotate("false positives\nremoved here", xy=(0.3, (fb[0] + fa[0]) / 2),
                xytext=(1.2, 0.30), fontsize=9.5, color=AFTER,
                arrowprops=dict(arrowstyle="->", color=AFTER))
    ax.legend(frameon=False, fontsize=9.5, loc="lower right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{out}.{ext}", dpi=150, bbox_inches="tight")
    print(f"wrote {out}.png/.pdf (selectivity; thr={thr})")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--native", required=True)
    ap.add_argument("--mode", choices=["dist", "selectivity"], default="dist")
    ap.add_argument("--out", default=None)
    ap.add_argument("--thr", type=float, default=10, help="heavy-penalty threshold")
    args = ap.parse_args()
    tc, before, after = load(args.native)
    out = args.out or f"fp_{args.mode}"
    if args.mode == "dist":
        plot_dist(tc, before, after, out, int(args.thr))
    else:
        plot_selectivity(tc, before, after, out, int(args.thr))


if __name__ == "__main__":
    main()
