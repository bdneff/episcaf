#!/usr/bin/env python3
"""
plot_cylinder_before_after.py

Before/after scatter of AlphaFold antibody clashes (af3_n_clash_res, x) vs cylinder
count (y). Left panel = original cylinder, right panel = native antigen volume
subtracted. The shaded box is the false-positive zone John pointed at -- designs with
little/no real antibody clash but a high cylinder count -- and each panel is annotated
with how many designs fall in it, so the emptying of that box is the headline.

Runs on the add_native_cylinder output; only the gated (carved) designs are plotted,
which is the same population the original scatter showed.

    python scripts/plot_cylinder_before_after.py \
        --native runs/run_rfd3_mpnn/04_filter/metrics_native_cyl.csv \
        --out cylinder_before_after
"""
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--native", required=True, help="add_native_cylinder output csv")
    ap.add_argument("--out", default="cylinder_before_after")
    ap.add_argument("--clash_max", type=float, default=0,
                    help="x-axis upper limit (0 = auto, 98th pct of true clash)")
    ap.add_argument("--fp_clash", type=float, default=3, help="zone: true clash <=")
    ap.add_argument("--fp_cyl", type=float, default=5, help="zone: cylinder >=")
    args = ap.parse_args()

    d = pd.read_csv(args.native, low_memory=False)
    tc = pd.to_numeric(d["af3_n_clash_res"], errors="coerce")
    before = pd.to_numeric(d["cylinder_ca_clashes"], errors="coerce")
    after = pd.to_numeric(d["cylinder_native_aware"], errors="coerce")
    # keep only designs that were actually carved (native-aware computed = gated)
    m = tc.notna() & before.notna() & after.notna()
    tc, before, after = tc[m].to_numpy(), before[m].to_numpy(), after[m].to_numpy()
    n = len(tc)

    xmax = args.clash_max or float(np.percentile(tc, 98))
    xmax = max(xmax, args.fp_clash + 1)
    ymax = float(max(before.max(), after.max())) * 1.05

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharex=True, sharey=True)
    for i, (ax, title, yv) in enumerate([(axes[0], "Original cylinder", before),
                          (axes[1], "Native antigen volume subtracted", after)]):
        # false-positive zone: shade first so points sit on top
        ax.add_patch(Rectangle((-0.5, args.fp_cyl), args.fp_clash + 0.5,
                               ymax - args.fp_cyl, facecolor="#d24a4a", alpha=0.08,
                               edgecolor="#d24a4a", linewidth=1.0, linestyle="--",
                               zorder=0))
        ax.scatter(tc, yv, s=5, alpha=0.12, color="#3b5b7a",
                   edgecolors="none", rasterized=True, zorder=1)
        nz = int(((tc <= args.fp_clash) & (yv >= args.fp_cyl)).sum())
        ax.text(0.96, 0.95, f"{nz:,} designs\nin the zone", transform=ax.transAxes,
                ha="right", va="top", fontsize=10.5, color="#a82a2a",
                fontweight="bold", linespacing=1.3)
        if i == 0:  # name the zone once
            ax.text(args.fp_clash + 1.2, ymax * 0.93,
                    "false-positive zone:\nno real clash, high cylinder",
                    fontsize=8.5, color="#a82a2a", va="top", linespacing=1.2)
        ax.set_title(title, fontsize=11.5, pad=8)
        ax.set_xlabel("AlphaFold antibody clashes")
        ax.set_xlim(-0.5, xmax)
        ax.set_ylim(-1, ymax)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylabel("cylinder count\n(scaffold residues in the antibody approach path)")
    fig.suptitle("Designs the antibody tolerates (low clash) but the cylinder penalizes "
                 f"(high count)   \u2014   n = {n:,} gated designs",
                 fontsize=11, y=1.02)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(f"{args.out}.{ext}", dpi=150, bbox_inches="tight")
    print(f"wrote {args.out}.png and {args.out}.pdf  (n={n})")


if __name__ == "__main__":
    main()
