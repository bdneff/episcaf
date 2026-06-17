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

LOCAL = Path("/Users/bneff/Desktop/projects/episcaf/known_antigen/analysis")
RFD1  = LOCAL/"full_run/metrics_full_rfd1_mpnn_LAWSON.csv"
RFD3  = LOCAL/"data/metrics_cylinder_full.csv"
OUT   = Path(__file__).resolve().parents[2]/"manuscript/figures/rfd1_vs_rfd3.png"

def num(s): return pd.to_numeric(s, errors="coerce")
def clash_count(s):  # parse "[1 2 3]" / "[]" -> integer count
    s = str(s).strip().strip("[]").strip()
    return 0 if not s or s.lower()=="nan" else len(s.split())

r1 = pd.read_csv(RFD1, usecols=["epitope_chunk_rmsd_vs_mpnn","overall_rmsd","mean_pae","af3_clash_resindices"], low_memory=False)
r3 = pd.read_csv(RFD3, low_memory=False)

r1d = pd.DataFrame({"epi":num(r1.epitope_chunk_rmsd_vs_mpnn),"overall":num(r1.overall_rmsd),
                    "pae":num(r1.mean_pae),"clash":r1.af3_clash_resindices.map(clash_count),
                    "clash_na":r1.af3_clash_resindices.isna()})
r3d = pd.DataFrame({"epi":num(r3.epitope_chunk_rmsd),"overall":num(r3.overall_rmsd),
                    "pae":num(r3.mean_pae),"clash":num(r3.af3_n_clash_res),"clash_na":num(r3.af3_n_clash_res).isna()})

def passrate(d):
    valid = d.epi.notna()&d.overall.notna()&d.pae.notna()&(~d.clash_na)
    p = valid&(d.epi<=1)&(d.overall<=2)&(d.pae<5)&(d.clash==0)
    return int(p.sum()), len(d), 100*p.sum()/len(d)
p1,n1,rate1 = passrate(r1d); p3,n3,rate3 = passrate(r3d)

PAN=[("epi","Epitope RMSD (A)",(0,16)),("overall","Overall RMSD (A)",(0,25)),
     ("pae","Mean PAE",(0,20)),("clash","AF3 clashing res",(0,60))]
C1,C3="#1f77b4","#d62728"
fig,axes=plt.subplots(1,4,figsize=(20,4.6))
for ax,(key,label,(lo,hi)) in zip(axes,PAN):
    e=np.arange(lo,hi+2,2) if key=="clash" else np.linspace(lo,hi,45)
    a=np.clip(r1d[key].dropna().values,lo,hi); b=np.clip(r3d[key].dropna().values,lo,hi)
    ax.hist(a,bins=e,density=True,color=C1,alpha=0.45,label=f"RFD1+MPNN")
    ax.hist(b,bins=e,density=True,histtype="step",color=C3,lw=2,label="RFD3+MPNN")
    ax.set_title(label,fontsize=12,fontweight="bold"); ax.set_yticks([]); ax.set_xlim(lo,hi)
    for s in ("top","right"): ax.spines[s].set_visible(False)
axes[0].legend(frameon=False,fontsize=10,loc="upper right")
fig.suptitle(f"RFD1 vs RFD3 on the same DP3 epitopes (both + ProteinMPNN + AlphaFold3): "
             f"nearly identical four-filter pass rate\n"
             f"RFD1+MPNN {p1:,}/{n1:,} = {rate1:.2f}%      RFD3+MPNN {p3:,}/{n3:,} = {rate3:.2f}%",
             fontsize=13,fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.9])
fig.savefig(OUT,dpi=140,bbox_inches="tight")
print(f"wrote {OUT}")
print(f"RFD1+MPNN pass {p1}/{n1} = {rate1:.3f}%   RFD3+MPNN pass {p3}/{n3} = {rate3:.3f}%")
