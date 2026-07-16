#!/usr/bin/env python3
"""
plot_rfd1_vs_rfd3.py  ->  manuscript/figures/rfd1_vs_rfd3.png

Direct comparison of the two backbone generators on the same DP3 (Lawson 59) epitopes,
both followed by ProteinMPNN + AlphaFold3. Overlays the per-metric distributions and
reports the four-filter pass rate, which is nearly identical between the two.

Filters (all four): epitope RMSD <= 1, overall RMSD <= 2, mean PAE < 5, zero AF3 clashes.

Inputs (local copies off-cluster):
  RFD1 (Lawson dp2) : metrics_full_rfd1_mpnn_LAWSON.csv  (epitope col = epitope_chunk_rmsd_vs_mpnn,
                       clash = len(af3_clash_resindices list))
  RFD3              : metrics_cylinder_full.csv           (epitope_chunk_rmsd, af3_n_clash_res)
"""
import numpy as np, pandas as pd
from pathlib import Path
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 14, "axes.titlesize": 16, "axes.labelsize": 15, "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12, "figure.titlesize": 18})  # paper-legible fonts

import sys
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from configs.paths import LOCAL_METRICS  # noqa: E402
RFD1  = LOCAL_METRICS["rfd1_mpnn_lawson"]
RFD3  = LOCAL_METRICS["rfd3_cylinder_full"]
OUT   = ROOT/"manuscript/figures/rfd1_vs_rfd3.png"

def num(s): return pd.to_numeric(s, errors="coerce")
def clash_count(x):
    # real NaN (no AF3 result) -> NaN; "[]" -> 0; "[1 2 3]" -> 3
    if pd.isna(x): return np.nan
    return len(str(x).strip().strip("[]").split())

r1 = pd.read_csv(RFD1, usecols=["epitope_chunk_rmsd_vs_mpnn","overall_rmsd","mean_pae","af3_clash_resindices"], low_memory=False)
r3 = pd.read_csv(RFD3, low_memory=False)

r1d = pd.DataFrame({"epi":num(r1.epitope_chunk_rmsd_vs_mpnn),"overall":num(r1.overall_rmsd),
                    "pae":num(r1.mean_pae),"clash":r1.af3_clash_resindices.map(clash_count)})
r3d = pd.DataFrame({"epi":num(r3.epitope_chunk_rmsd),"overall":num(r3.overall_rmsd),
                    "pae":num(r3.mean_pae),"clash":num(r3.af3_n_clash_res)})

def stats(d, n_total):
    valid = d.epi.notna()&d.overall.notna()&d.pae.notna()&d.clash.notna()
    p = valid&(d.epi<=1)&(d.overall<=2)&(d.pae<5)&(d.clash==0)
    nv = int(valid.sum())
    return int(p.sum()), n_total, nv, 100*p.sum()/nv  # rate over VALID predictions
p1,n1,nv1,rate1 = stats(r1d, len(r1)); p3,n3,nv3,rate3 = stats(r3d, len(r3))

PAN=[("epi","Epitope RMSD (A)",(0,16)),("overall","Overall RMSD (A)",(0,25)),
     ("pae","Mean PAE",(0,20)),("clash","AF3 clashing res",(0,60))]
C1,C3="#1f77b4","#d62728"
fig,axes=plt.subplots(1,4,figsize=(13,4.0))
for ax,(key,label,(lo,hi)) in zip(axes,PAN):
    e=np.arange(lo,hi+2,2) if key=="clash" else np.linspace(lo,hi,45)
    a=np.clip(r1d[key].dropna().values,lo,hi); b=np.clip(r3d[key].dropna().values,lo,hi)
    ax.hist(a,bins=e,density=True,color=C1,alpha=0.45,label=f"RFD1+MPNN (n={nv1:,})")
    ax.hist(b,bins=e,density=True,histtype="step",color=C3,lw=2,label=f"RFD3+MPNN (n={nv3:,})")
    ax.set_title(label,fontsize=17,fontweight="bold"); ax.set_yticks([]); ax.set_xlim(lo,hi)
    for s in ("top","right"): ax.spines[s].set_visible(False)
axes[0].legend(frameon=False,fontsize=13,loc="upper right")
# Pass-rate stats live in the caption; keep the figure title short and legible.
print(f"[rfd1_vs_rfd3] pass per design generated: RFD1 {p1:,}/{n1:,}={100*p1/n1:.2f}%  "
      f"RFD3 {p3:,}/{n3:,}={100*p3/n3:.2f}%  (missing>=1 metric -> non-pass: RFD1 {n1-nv1:,}, "
      f"RFD3 {n3-nv3:,}; over valid {rate1:.2f}% vs {rate3:.2f}%)")
fig.suptitle("RFD1+MPNN vs RFD3+MPNN on the same DP3 epitopes", fontsize=18, fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.92])
fig.savefig(OUT,dpi=140,bbox_inches="tight")
print(f"wrote {OUT}")
print(f"RFD1: valid {nv1:,}/{n1:,}  pass {p1}  per-valid {rate1:.3f}%  per-total {100*p1/n1:.3f}%")
print(f"RFD3: valid {nv3:,}/{n3:,}  pass {p3}  per-valid {rate3:.3f}%  per-total {100*p3/n3:.3f}%")
