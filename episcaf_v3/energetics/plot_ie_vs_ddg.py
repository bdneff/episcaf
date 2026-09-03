#!/usr/bin/env python3
"""Validate per-residue interaction energy against experimental ddG (Decision D2, on 3HFM).

Reads holo_ie_mean.csv (from holo_ie.py) and skempi_3hfm_ddg.csv (the antigen alanine-scan ground
truth). Produces two views and a false-positive report:

  Panel A  per-residue overlay -- residue number on x; favorable interaction energy (-ab_total) on
           the LEFT y; experimental ddG on the RIGHT y (markers at measured residues). Shows whether
           the two track along the sequence.
  Panel B  the actual test -- scatter of favorable IE (x) vs ddG (y) for the measured residues, with
           the FALSE-POSITIVE zone (high IE, low ddG) shaded. Jacob's criterion: that zone must be
           EMPTY. High IE with low measured ddG = a residue we would wrongly lock in design.
           A hot spot with low IE (upper-left) is an acceptable false negative -- it is likely a
           Category-2 structural/orientational residue that a per-residue energy scan cannot see.

The point is an empirical fit: tune the IE method (cutoff, terms, with/without the solvent channel)
to eliminate false positives on 3HFM, then test generalization on other complexes. So the headline
number is PRECISION of the high-IE calls, not correlation.

Run (anywhere with pandas/numpy/scipy/matplotlib):
    python plot_ie_vs_ddg.py --ie <path>/holo_ie_mean.csv --ddg energetics/skempi_3hfm_ddg.csv \
        --out manuscript/figures/ie_vs_ddg_3hfm.png
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

IE_BLUE = "#2a6f97"      # interaction energy
DDG_RED = "#d1495b"      # experimental ddG
FP_SHADE = "#f2c14e"     # false-positive zone


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--ie", required=True, help="holo_ie_mean.csv from holo_ie.py")
    ap.add_argument("--ddg", default="skempi_3hfm_ddg.csv", help="antigen alanine ddG ground truth")
    ap.add_argument("--channel", default="ab_total",
                    help="which energy column to test (ab_total | ab_coul | ab_lj; add sol_* to try solvent)")
    ap.add_argument("--hot", type=float, default=1.0, help="ddG hot-spot threshold (kcal/mol)")
    ap.add_argument("--ie-pct", type=float, default=90.0,
                    help="'high IE' percentile (over all antigen residues) for the false-positive test")
    ap.add_argument("--out", default="ie_vs_ddg_3hfm.png")
    args = ap.parse_args()

    ie = pd.read_csv(args.ie)
    ddg = pd.read_csv(args.ddg, comment="#")
    # favorable interaction energy: flip sign so hot spots point UP (favorable = large positive)
    ie["fav_ie"] = -ie[args.channel]
    ie_thr = np.percentile(ie["fav_ie"], args.ie_pct)

    m = ie.merge(ddg, left_on="ag_res_idx", right_on="resid", how="inner")  # measured residues only
    print(f"antigen residues: {len(ie)}; with experimental ddG: {len(m)}")
    print(f"'high IE' cut = {args.ie_pct:.0f}th pct of -{args.channel} = {ie_thr:.1f} kJ/mol; "
          f"hot-spot cut = ddG > {args.hot}")

    # correlation (context, not the headline)
    rho, prho = stats.spearmanr(m["fav_ie"], m["ddg_kcal_mol"])
    r, pr = stats.pearsonr(m["fav_ie"], m["ddg_kcal_mol"])
    print(f"\nover the {len(m)} measured residues:  Spearman rho={rho:+.2f} (p={prho:.2g})  "
          f"Pearson r={r:+.2f}")

    # the false-positive test: measured residues above the IE cut but NOT hot spots
    fp = m[(m["fav_ie"] > ie_thr) & (m["ddg_kcal_mol"] < args.hot)]
    tp = m[(m["fav_ie"] > ie_thr) & (m["ddg_kcal_mol"] >= args.hot)]
    n_high = len(m[m["fav_ie"] > ie_thr])
    print(f"\nHIGH-IE measured residues: {n_high}  ->  true positives {len(tp)}, "
          f"FALSE POSITIVES {len(fp)}")
    if len(fp):
        print("  !! FALSE POSITIVES (high IE, low ΔΔG -- the thing to eliminate):")
        for _, x in fp.iterrows():
            print(f"     {x.resname_x}{int(x.ag_res_idx)}  fav_ie={x.fav_ie:.1f}  ΔΔG={x.ddg_kcal_mol:+.2f}")
    else:
        print("  OK -- no false positives at this cut (every high-IE measured residue is a hot spot).")
    print(f"  precision of high-IE calls = {len(tp)}/{n_high}" if n_high else "  (no high-IE residues)")

    print("\ntop 10 antigen residues by favorable IE (ΔΔG shown where measured):")
    top = ie.sort_values("fav_ie", ascending=False).head(10).merge(ddg, left_on="ag_res_idx",
                                                                    right_on="resid", how="left")
    for _, x in top.iterrows():
        d = f"{x.ddg_kcal_mol:+.2f}" if pd.notna(x.get("ddg_kcal_mol")) else "  (no value)"
        print(f"    {x.resname_x}{int(x.ag_res_idx):<4d} fav_ie={x.fav_ie:7.1f} kJ/mol   ΔΔG={d}")

    # ---- figure ----
    fig, (axA, axB) = plt.subplots(2, 1, figsize=(11, 9))

    # Panel A: per-residue overlay, dual axis
    axA.fill_between(ie["ag_res_idx"], 0, ie["fav_ie"], color=IE_BLUE, alpha=0.30, step="mid")
    axA.plot(ie["ag_res_idx"], ie["fav_ie"], color=IE_BLUE, lw=1.0)
    axA.axhline(ie_thr, color=IE_BLUE, ls=":", lw=1, alpha=0.7)
    axA.set_xlabel("antigen (lysozyme) residue number")
    axA.set_ylabel(f"favorable interaction energy  ($-${args.channel}, kJ/mol)", color=IE_BLUE)
    axA.tick_params(axis="y", labelcolor=IE_BLUE)
    axA2 = axA.twinx()
    axA2.stem(m["ag_res_idx"], m["ddg_kcal_mol"], linefmt=DDG_RED, markerfmt="o", basefmt=" ")
    axA2.set_ylabel("experimental $\\Delta\\Delta G$ (kcal/mol)", color=DDG_RED)
    axA2.tick_params(axis="y", labelcolor=DDG_RED)
    for _, x in m[m["ddg_kcal_mol"] > 1.5].iterrows():
        axA2.annotate(f"{x.resname_x}{int(x.ag_res_idx)}", (x.ag_res_idx, x.ddg_kcal_mol),
                      textcoords="offset points", xytext=(0, 4), ha="center", fontsize=9, color=DDG_RED)
    axA.set_title("Per-residue: interaction energy vs experimental $\\Delta\\Delta G$")

    # Panel B: the false-positive test
    from matplotlib.patches import Rectangle
    ylo = m["ddg_kcal_mol"].min() - 0.6
    yhi = m["ddg_kcal_mol"].max() + 0.6
    xlo = m["fav_ie"].min() - 15
    xhi = m["fav_ie"].max() * 1.08
    axB.set_xlim(xlo, xhi); axB.set_ylim(ylo, yhi)
    # shade ONLY the lower-right box: high IE (> ie_thr) AND low ddG (< hot) = false positives
    axB.add_patch(Rectangle((ie_thr, ylo), xhi - ie_thr, args.hot - ylo,
                            color=FP_SHADE, alpha=0.18, zorder=0))
    axB.axhline(args.hot, color="0.5", ls="--", lw=1)
    axB.axvline(ie_thr, color=IE_BLUE, ls=":", lw=1)
    axB.scatter(m["fav_ie"], m["ddg_kcal_mol"], s=40, color="0.35", zorder=3)
    if len(fp):
        axB.scatter(fp["fav_ie"], fp["ddg_kcal_mol"], s=80, color=DDG_RED, zorder=4,
                    label="false positive")
    for _, x in m.iterrows():
        axB.annotate(f"{x.resname_x}{int(x.ag_res_idx)}", (x.fav_ie, x.ddg_kcal_mol),
                     textcoords="offset points", xytext=(4, 3), fontsize=8)
    axB.set_xlabel(f"favorable interaction energy  ($-${args.channel}, kJ/mol)")
    axB.set_ylabel("experimental $\\Delta\\Delta G$ (kcal/mol)")
    axB.set_title(f"False-positive test (shaded zone must be empty)   Spearman $\\rho$={rho:+.2f}")
    axB.text(0.02, 0.95, "high IE + low $\\Delta\\Delta G$ =\nwrongly locked in design", transform=axB.transAxes,
             fontsize=9, va="top", color=DDG_RED)
    if len(fp):
        axB.legend(loc="lower right")

    fig.tight_layout()
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=150)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
