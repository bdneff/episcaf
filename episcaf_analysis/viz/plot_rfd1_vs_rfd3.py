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
fig,axes=plt.subplots(1,4,figsize=(20,4.6))
for ax,(key,label,(lo,hi)) in zip(axes,PAN):
    e=np.arange(lo,hi+2,2) if key=="clash" else np.linspace(lo,hi,45)
    a=np.clip(r1d[key].dropna().values,lo,hi); b=np.clip(r3d[key].dropna().values,lo,hi)
    ax.hist(a,bins=e,density=True,color=C1,alpha=0.45,label=f"RFD1+MPNN (n={nv1:,})")
    ax.hist(b,bins=e,density=True,histtype="step",color=C3,lw=2,label=f"RFD3+MPNN (n={nv3:,})")
    ax.set_title(label,fontsize=12,fontweight="bold"); ax.set_yticks([]); ax.set_xlim(lo,hi)
    for s in ("top","right"): ax.spines[s].set_visible(False)
axes[0].legend(frameon=False,fontsize=9.5,loc="upper right")
fig.suptitle(
    "RFD1 vs RFD3 on the same DP3 epitopes (both + ProteinMPNN + AlphaFold3); distributions density-normalized over designs with all four metrics\n"
    f"four-filter pass rate per design generated: RFD1+MPNN {p1:,}/{n1:,} = {100*p1/n1:.2f}%   "
    f"RFD3+MPNN {p3:,}/{n3:,} = {100*p3/n3:.2f}%   "
    f"(designs missing >=1 metric count as non-pass: RFD1 {n1-nv1:,}, RFD3 {n3-nv3:,}; over valid only {rate1:.2f}% vs {rate3:.2f}%)",
    fontsize=11.5,fontweight="bold")
fig.tight_layout(rect=[0,0,1,0.88])
fig.savefig(OUT,dpi=140,bbox_inches="tight")
print(f"wrote {OUT}")
print(f"RFD1: valid {nv1:,}/{n1:,}  pass {p1}  per-valid {rate1:.3f}%  per-total {100*p1/n1:.3f}%")
print(f"RFD3: valid {nv3:,}/{n3:,}  pass {p3}  per-valid {rate3:.3f}%  per-total {100*p3/n3:.3f}%")
