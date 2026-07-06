#!/usr/bin/env python3
"""
plot_topn_overlay.py  ->  scratch_figs/topn_overlay.png   (review copy; not in manuscript yet)

John's Figure-7 variant. Same rows-x-metrics layout as passes_overlay, but:
  * highlight the SELECTED designs (top-n per epitope by composite), not the four-filter passes;
  * show n = 1, 10, 20 as OVERLAPPING traces so you can read how the selected distribution
    tightens as you take fewer per epitope;
  * SINGLE axis, DENSITY-normalized (not twin-axis counts) so the shapes are comparable across
    n and against the DP3 reference -- n=1 has ~1 design/epitope and n=20 ~20, so counts would
    make the small-n traces vanish.
  * top row = the DP3 mAb designs that actually succeeded (four-filter passes) -- the reference
    the DP4 top-n traces are meant to be projected against.

Ordered n -> sequential single-hue blue ramp (dark=tightest n=1 .. light=n=20); DP3 reference in
crimson. Validated: adjacent-CVD separation passes; linewidth + legend give the light end relief.

Only the DP4 polyclonal-tiling category is local (1d2k/4wat/6m0j). The other categories John lists
(single-island #1 on the cluster; sampling/mutant/whole-epitope not generated yet) get a row each
once their metrics exist -- the script takes any (label, all_df, scored_df) row.

Inputs (local off-cluster copies):
  DP3 all   : metrics_native_cyl_full.csv     (four-filter pass mask computed here)
  DP4 all   : metrics_12mer.csv               (full population, grey)
  DP4 ranks : composite_12mer_allscored.csv   (rank_in_epitope 1..192 -> top-n selection)
"""
import argparse
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
plt.rcParams.update({"font.size": 14, "axes.titlesize": 16, "axes.labelsize": 15,
                     "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 13,
                     "figure.titlesize": 18})

ap = argparse.ArgumentParser()
ap.add_argument("--ns", default="1,10,20", help="comma list of per-epitope selection depths")
ap.add_argument("--out", default="")
args = ap.parse_args()
NS = [int(x) for x in args.ns.split(",")]                # selection depths (ordered small->large)

LOCAL = Path("/Users/bneff/Desktop/projects/episcaf")
DP3   = LOCAL/"known_antigen/analysis/data/metrics_native_cyl_full.csv"
DP4   = LOCAL/"12mer_tiling/analysis/data/metrics_12mer.csv"
DP4RK = LOCAL/"12mer_tiling/analysis/data/composite_12mer_allscored.csv"
OUT   = Path(args.out) if args.out else (
        Path(__file__).resolve().parents[2]/f"scratch_figs/topn_overlay_{'_'.join(map(str,NS))}.png")

# sequential blue ramp, dark=tightest n .. light=most-relaxed n (ordered encoding of n)
_pos    = np.linspace(0.95, 0.45, len(NS))               # dark->light across NS
N_COLOR = {n: plt.cm.Blues(p) for n, p in zip(NS, _pos)}
N_LW    = {n: lw for n, lw in zip(NS, np.linspace(2.4, 1.4, len(NS)))}  # thin the relaxed traces
DP3_C   = "#c0392b"                                       # reference (designs that succeeded)
POP_C   = "#d9d9d9"

def num(s): return pd.to_numeric(s, errors="coerce")

dp3 = pd.read_csv(DP3, low_memory=False)
ispass = ((num(dp3.epitope_chunk_rmsd)<=1)&(num(dp3.overall_rmsd)<=2)
          &(num(dp3.mean_pae)<5)&(num(dp3.af3_n_clash_res)==0))
m12 = pd.read_csv(DP4, low_memory=False)
rk  = pd.read_csv(DP4RK, low_memory=False)

COLS=[("epitope_chunk_rmsd","Epitope\nRMSD (A)"),("mean_pae","Mean PAE\n(global)"),
      ("epitope_pae","Epitope\nPAE"),("overall_rmsd","Overall\nRMSD (A)"),("ptm","pTM"),
      ("af3_n_clash_res","AF3 clash\n(real Ab)"),
      ("cylinder_ca_clashes","Cylinder\n(plain)"),
      ("cylinder_native_aware","Cylinder\n(native-aware)")]

# rows: (label, full-population df, {n: selected df} or None, is_dp3_reference)
rows=[("DP3 mAb\n(succeeded)", dp3, {"pass": dp3[ispass]}, True)]
for ag in ["1d2k","4wat","6m0j"]:
    sel={n: rk[(rk.antigen==ag)&(num(rk.rank_in_epitope)<=n)] for n in NS}
    rows.append((f"{ag.upper()} 12mer\n(DP4 polyclonal)", m12[m12.antigen==ag], sel, False))

xlim={}
for key,_ in COLS:
    v=[num(a[key]).dropna().values for _,a,_,ref in rows if (key in a and not (key=="af3_n_clash_res" and not ref))]
    v=np.concatenate(v) if v else np.array([0,1])
    lo=max(0.15,np.percentile(v,0.5)) if key=="ptm" else 0
    xlim[key]=(lo,np.percentile(v,99.5))

def edges(key,lo,hi):
    if key in ("cylinder_ca_clashes","cylinder_native_aware","af3_n_clash_res"):
        return np.arange(np.floor(lo),np.ceil(hi)+2,2)
    if key=="ptm": return np.linspace(lo,hi,26)
    return np.linspace(lo,hi,34)

fig,axes=plt.subplots(len(rows),len(COLS),figsize=(2.35*len(COLS),2.15*len(rows)+1.2))
for r,(name,alld,sel,ref) in enumerate(rows):
    for c,(key,label) in enumerate(COLS):
        ax=axes[r,c]
        if key=="af3_n_clash_res" and not ref:
            ax.text(0.5,0.5,"n/a\n(no antibody)",ha="center",va="center",style="italic",
                    color="#999",transform=ax.transAxes,fontsize=13)
            ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]; continue
        lo,hi=xlim[key]; e=edges(key,lo,hi)
        xa=np.clip(num(alld[key]).dropna().values,lo,hi)
        ax.hist(xa,bins=e,color=POP_C,density=True,zorder=1)           # full population, density
        if ref:
            xp=np.clip(num(sel["pass"][key]).dropna().values,lo,hi)
            if len(xp):
                ax.hist(xp,bins=e,density=True,histtype="step",color=DP3_C,lw=2.4,zorder=5)
                ax.axvline(np.median(xp),ls="--",color=DP3_C,lw=1.0,zorder=4)
        else:
            for n in NS:
                d=sel[n]
                if key not in d: continue
                xn=np.clip(num(d[key]).dropna().values,lo,hi)
                if len(xn):
                    ax.hist(xn,bins=e,density=True,histtype="step",
                            color=N_COLOR[n],lw=N_LW[n],zorder=3+ (0 if n==20 else n//10))
        ax.set_xlim(lo,hi); ax.tick_params(axis='x',labelsize=8)
        ax.set_yticks([])                                              # density units carry no meaning to read
        for s in ("top","right","left"): ax.spines[s].set_visible(False)
        if r==0: ax.set_title(label,fontsize=16,fontweight="bold")
        if c==0:
            npass = len(sel["pass"]) if ref else None
            sub = f"succeeded n={npass:,}" if ref else f"pop n={len(alld):,}"
            ax.set_ylabel(f"{name}\n{sub}", color=(DP3_C if ref else "#111"),
                          fontweight="bold",rotation=0,ha="right",va="center",labelpad=52,fontsize=13)

handles=[Patch(fc=POP_C,ec="none",label="full population (density)"),
         Line2D([0],[0],color=DP3_C,lw=2.4,label="DP3: designs that succeeded (four-filter)")]
handles+=[Line2D([0],[0],color=N_COLOR[n],lw=N_LW[n]+0.6,label=f"DP4 selected: top {n} / epitope") for n in NS]
fig.legend(handles=handles,loc="lower center",ncol=5,frameon=False,
           bbox_to_anchor=(0.5,-0.015),fontsize=13)
fig.suptitle("Where the SELECTED designs land, per metric: DP4 top-n per epitope "
             f"(blue, n={'/'.join(map(str,NS))}) vs the DP3 designs that succeeded (red)\n"
             "each trace density-normalized on one axis so shapes are comparable across n and rows; "
             "dashed = median of the DP3 successes",
             fontsize=15,fontweight="bold")
fig.tight_layout(rect=[0.05,0.05,1,0.93])
fig.savefig(OUT,dpi=135,bbox_inches="tight")
print("wrote",OUT)
for r,(name,alld,sel,ref) in enumerate(rows):
    if ref:
        print(f"  {name.splitlines()[0]:10s} succeeded n={len(sel['pass']):,}")
    else:
        print(f"  {name.splitlines()[0]:10s} pop n={len(alld):,}  "+
              "  ".join(f"top{n}={len(sel[n]):,}" for n in NS))
