#!/usr/bin/env python3
"""
add_native_cylinder.py

Native-antigen-aware cylinder clash metric + before/after visualization.

For each design:
  - build the locked cylinder over the AF3 epitope (same as add_cylinder_metric);
  - load the design's native complex ({id}.pdb in --native_dir; antigen = chain A,
    antibody = chains B/C, matching the Fab-probe conventions);
  - locate the design epitope inside the native antigen by sequence and Kabsch-align
    the native antigen onto the AF3 epitope;
  - compute three numbers:
        cylinder_ca_clashes     plain scaffold-CA count inside the cylinder
        native_in_cylinder      TEST: native (non-epitope) CAs inside the cylinder
                                -> how much nature-allowed volume the cylinder covers
        cylinder_native_aware   SOLUTION: scaffold CAs inside the cylinder but NOT
                                within --exclude_dist of a native (non-epitope)
                                heavy atom (carve the native footprint out)

With --viz_token TKN (or --viz_first) it also writes a 3D before/after figure for
that design: the cylinder with all flagged scaffold residues, then the same cylinder
with the residues sitting in native volume carved out.

Only the native antigen is needed (the epitope is taken from it), so this works on
no-antibody targets too.

    python scripts/add_native_cylinder.py \
        --metrics_csv runs/run_rfd3_mpnn/04_filter/metrics_cylinder_full.csv \
        --dp2_parquet datasets/dp2.parquet \
        --native_dir  /tgen_labs/.../abdb/complex_pdbfiles/cleaned \
        --out_csv     runs/run_rfd3_mpnn/04_filter/metrics_native_cyl.csv \
        --limit 500 --viz_first
"""
import argparse, difflib, gzip, math
from pathlib import Path
from typing import List, Optional

import gemmi
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection  # noqa

RADIUS, HEIGHT, OFFSET = 16.0, 40.0, -4.0
AB_CHAINS = ("B", "C")

TOKEN_COLS=["token","assay_scaffolded_epitope_id"]; AF3_DIR_COLS=["af3_dir"]
ID_COLS=["id"]; DP2_EPI_COLS=["assay_scaffolded_epitope_chunk_resindices"]
DP2_EPI_MPNN=["scaffolded_epitope_chunk_resindices"]; WS_COLS=["af3_window_start"]
TRUE_CLASH=["af3_n_clash_res"]; EPI_RMSD_COLS=["epitope_chunk_rmsd"]

THREE2ONE = {  # for epitope sequence matching
 'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
 'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
 'THR':'T','TRP':'W','TYR':'Y','VAL':'V'}


# ---------------------------------------------------------------- IO helpers
def parse_index_list(x) -> List[int]:
    if x is None or (isinstance(x, float) and math.isnan(x)): return []
    if isinstance(x,(list,tuple,np.ndarray)): return [int(i) for i in x]
    s=str(x).replace("[","").replace("]","").replace(","," ")
    out=[]
    for t in s.split():
        try: out.append(int(t))
        except ValueError: pass
    return out

def first_present(cols,cands):
    for c in cands:
        if c in cols: return c
    return None

def read_gemmi(p: Path):
    s=str(p)
    if s.endswith(".gz"):
        with gzip.open(p,"rt") as f: doc=gemmi.cif.read_string(f.read())
        return gemmi.make_structure_from_block(doc.sole_block())
    if p.suffix.lower()==".cif":
        return gemmi.make_structure_from_block(gemmi.cif.read(str(p)).sole_block())
    return gemmi.read_structure(str(p))

def chain_by_name(model,name):
    for ch in model:
        if ch.name==name: return ch
    return None

def find_af3_cif(af3_dir: Path):
    for pat in ("*_model.cif","*_model.cif.gz","model.cif","model.cif.gz"):
        h=next(af3_dir.glob(pat),None) or next(af3_dir.rglob(pat),None)
        if h: return h
    return None

def load_af3(af3_dir: Path):
    """chain-A: CA coords, residue index, 1-letter resname."""
    cif=find_af3_cif(Path(af3_dir))
    if cif is None: return None
    try: st=read_gemmi(cif)
    except Exception: return None
    ch=chain_by_name(st[0],"A") or st[0][0]
    ca,ridx,seq=[],[],[]
    for i,res in enumerate(ch):
        a=res.find_atom("CA",altloc="*")
        if a:
            ca.append([a.pos.x,a.pos.y,a.pos.z]); ridx.append(i)
            seq.append(THREE2ONE.get(res.name.upper(),"X"))
    if not ca: return None
    return np.asarray(ca,float),np.asarray(ridx,int),"".join(seq)

def load_native_antigen(pdb: Path):
    """antigen (chain A): CA coords, 1-letter seq, heavy-atom coords + per-atom
    residue position (0..M-1)."""
    try: st=read_gemmi(pdb)
    except Exception: return None
    ch=chain_by_name(st[0],"A")
    if ch is None: return None
    ca,seq,heavy,heavy_resi=[],[],[],[]
    for ri,res in enumerate(ch):
        a=res.find_atom("CA",altloc="*")
        if not a: continue
        ca.append([a.pos.x,a.pos.y,a.pos.z]); seq.append(THREE2ONE.get(res.name.upper(),"X"))
        ridx=len(ca)-1
        for at in res:
            if at.element==gemmi.Element("H"): continue
            heavy.append([at.pos.x,at.pos.y,at.pos.z]); heavy_resi.append(ridx)
    if len(ca)<5: return None
    return (np.asarray(ca,float),"".join(seq),
            np.asarray(heavy,float),np.asarray(heavy_resi,int))


# ---------------------------------------------------------------- geometry
def cylinder_frame(epi_ca,all_ca):
    centroid=epi_ca.mean(0)
    _,_,Vt=np.linalg.svd(epi_ca-centroid); normal=Vt[-1]
    if np.dot(normal,all_ca.mean(0)-centroid)>0: normal=-normal
    return centroid, normal, centroid+OFFSET*normal

def inside_cylinder(pts,base,normal,r=RADIUS,h=HEIGHT):
    v=pts-base; proj=v@normal
    dist=np.linalg.norm(v-np.outer(proj,normal),axis=1)
    return (proj>=0)&(proj<=h)&(dist<=r)

def kabsch(P,Q):
    Pc,Qc=P-P.mean(0),Q-Q.mean(0)
    U,_,Vt=np.linalg.svd(Pc.T@Qc)
    d=np.sign(np.linalg.det(Vt.T@U.T)); D=np.diag([1,1,d])
    R=Vt.T@D@U.T
    return R, Q.mean(0)-R@P.mean(0)

def match_epitope(epi_seq,epi_ca,nseq,nca,min_match=5,max_rmsd=8.0):
    """Locate the design epitope inside the native antigen by sequence and align the
    antigen onto the design. Tolerant to crystal gaps / differing flanks: anchors on
    the longest matching block, Kabsch-aligns native->design on it, and maps the full
    epitope onto the antigen by that offset. Returns (R,t,native_epitope_idx_set,rmsd)
    or None if no acceptable match."""
    sm=difflib.SequenceMatcher(None,epi_seq,nseq,autojunk=False)
    a,b,size=sm.find_longest_match(0,len(epi_seq),0,len(nseq))
    if size<min_match:
        return None
    P=nca[b:b+size]                      # native CAs in the matched block
    Q=epi_ca[a:a+size]                   # design epitope CAs in the matched block
    R,t=kabsch(P,Q)
    rmsd=float(np.sqrt((((R@P.T).T+t-Q)**2).sum(1).mean()))
    if rmsd>max_rmsd:
        return None
    off=b-a                              # antigen_index = off + epitope_index
    nepi={off+k for k in range(len(epi_ca)) if 0<=off+k<len(nca)}
    return R,t,nepi,rmsd


# ---------------------------------------------------------------- viz
def _basis(normal):
    n=normal/np.linalg.norm(normal)
    a=np.array([1.,0,0]) if abs(n[0])<0.9 else np.array([0,1.,0])
    u=np.cross(n,a); u/=np.linalg.norm(u); v=np.cross(n,u)
    return n,u,v

def voxelize_cylinder(base,normal,r,h,step):
    """Voxel centers filling the cylinder interior."""
    n,u,v=_basis(normal)
    ts=np.arange(step/2,h,step)
    g=np.arange(-r+step/2,r,step)
    A,B=np.meshgrid(g,g); disk=(A*A+B*B)<=r*r
    Aa,Bb=A[disk],B[disk]
    return np.vstack([base+t*n+np.outer(Aa,u)+np.outer(Bb,v) for t in ts])

def make_viz(epi_ca,native_heavy,exclude_dist,base,normal,native_in_cyl,
             n_plain,n_aware,out_path,r=RADIUS,h=HEIGHT,step=2.5):
    vox=voxelize_cylinder(base,normal,r,h,step)
    if len(native_heavy):
        d,_=cKDTree(native_heavy).query(vox,k=1)
        vcarved=d<=exclude_dist
    else:
        vcarved=np.zeros(len(vox),bool)
    pct=100.0*vcarved.sum()/max(len(vox),1)

    fig=plt.figure(figsize=(12,5.8))
    panels=[("Cylinder volume (before)",False),
            (f"Native volume carved out (\u2212{pct:.0f}% of volume)",True)]
    for k,(title,after) in enumerate(panels):
        ax=fig.add_subplot(1,2,k+1,projection="3d")
        if not after:
            ax.scatter(*vox.T,s=22,c="#d62728",alpha=0.95,depthshade=False)
        else:
            ax.scatter(*vox[~vcarved].T,s=22,c="#d62728",alpha=0.95,depthshade=False)
            if vcarved.any():
                ax.scatter(*vox[vcarved].T,s=22,c="#d62728",alpha=0.10,
                           depthshade=False,label="carved (native-occupied)")
        ax.scatter(*epi_ca.T,s=40,c="#2ca02c",depthshade=False,label="epitope")
        if len(native_in_cyl):
            ax.scatter(*native_in_cyl.T,s=16,c="#1f77b4",alpha=0.7,
                       depthshade=False,label="native antigen CAs")
        ax.set_title(title,fontsize=11)
        ax.legend(loc="upper left",fontsize=8,frameon=False)
        ax.set_xlabel("x");ax.set_ylabel("y");ax.set_zlabel("z")
        ax.set_box_aspect((1,1,1)); ax.view_init(elev=14,azim=-60)
    fig.suptitle(f"flagged scaffold residues: {n_plain} \u2192 {n_aware}",
                 y=0.04,fontsize=10)
    fig.tight_layout()
    fig.savefig(out_path,dpi=150,bbox_inches="tight")
    plt.close(fig)


# ---------------------------------------------------------------- per design
def process(af3_dir,native_pdb,epi_ris,exclude_dist):
    a=load_af3(af3_dir)
    if a is None: return None
    ca,res_idx,seq=a
    if not epi_ris or max(epi_ris)>=len(ca): return None
    epi_pos=[p for p,ri in enumerate(res_idx) if ri in set(epi_ris)]
    if len(epi_pos)<3: return None
    epi_ca=ca[epi_pos]; epi_seq="".join(seq[p] for p in epi_pos)
    centroid,normal,base=cylinder_frame(epi_ca,ca)

    scaf_pos=[p for p in range(len(ca)) if res_idx[p] not in set(epi_ris)]
    scaf_ca=ca[scaf_pos]
    ins_scaf=inside_cylinder(scaf_ca,base,normal)
    n_plain=int(ins_scaf.sum())

    nat=load_native_antigen(Path(native_pdb)) if native_pdb else None
    native_in=0; n_aware=n_plain; carved=np.zeros(n_plain,bool)
    native_in_cyl_xyz=np.empty((0,3)); inside_xyz=scaf_ca[ins_scaf]
    native_heavy_al=np.empty((0,3))
    if nat is not None:
        nca,nseq,nheavy,nresi=nat
        m=match_epitope(epi_seq,epi_ca,nseq,nca)   # gap-tolerant locate + align
        if m is not None:
            R,t,nepi_set,_=m
            nca_al=(R@nca.T).T+t
            nheavy_al=(R@nheavy.T).T+t
            nonepi=np.ones(len(nca),bool); nonepi[list(nepi_set)]=False
            nin=inside_cylinder(nca_al,base,normal)&nonepi
            native_in=int(nin.sum()); native_in_cyl_xyz=nca_al[nin]
            # carve scaffold residues sitting in native (non-epitope) volume
            heavy_nonepi=nheavy_al[nonepi[nresi]]
            native_heavy_al=heavy_nonepi
            if n_plain>0 and len(heavy_nonepi):
                d,_=cKDTree(heavy_nonepi).query(inside_xyz,k=1)
                carved=d<=exclude_dist
                n_aware=int((~carved).sum())
    return dict(n_plain=n_plain,native_in=native_in,n_aware=n_aware,
                epi_ca=epi_ca,inside=inside_xyz,carved=carved,
                outside=scaf_ca[~ins_scaf],native_in_cyl=native_in_cyl_xyz,
                native_heavy=native_heavy_al,base=base,normal=normal)


def main():
    ap=argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics_csv",required=True)
    ap.add_argument("--dp2_parquet",required=True)
    ap.add_argument("--native_dir",required=True)
    ap.add_argument("--out_csv",required=True)
    ap.add_argument("--limit",type=int,default=0)
    ap.add_argument("--exclude_dist",type=float,default=4.0)
    ap.add_argument("--gate",type=float,default=0.0,
                    help="only process designs with epitope_chunk_rmsd < gate (0 = all)")
    ap.add_argument("--viz_token",default=None)
    ap.add_argument("--viz_first",action="store_true")
    ap.add_argument("--viz_out",default="native_cylinder_viz.png")
    args=ap.parse_args()

    df=pd.read_csv(args.metrics_csv,low_memory=False)
    if args.limit>0: df=df.head(args.limit).copy()
    tok=first_present(df.columns,TOKEN_COLS); dirc=first_present(df.columns,AF3_DIR_COLS)
    idc=first_present(df.columns,ID_COLS); wsc=first_present(df.columns,WS_COLS)

    dp2=pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"]=dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()
    epic=first_present(dp2.columns,DP2_EPI_COLS); epim=first_present(dp2.columns,DP2_EPI_MPNN)
    look=dp2.drop_duplicates("assay_scaffolded_epitope_id").set_index("assay_scaffolded_epitope_id")

    native_index={p.stem.lower(): p for p in Path(args.native_dir).glob("*.pdb")}
    erc=first_present(df.columns,EPI_RMSD_COLS)
    # preserve any existing cylinder_ca_clashes (so gate-skipped rows keep their value)
    plain=(pd.to_numeric(df["cylinder_ca_clashes"],errors="coerce").to_numpy()
           if "cylinder_ca_clashes" in df.columns else np.full(len(df),np.nan))
    natin=np.full(len(df),np.nan); aware=np.full(len(df),np.nan)
    did_viz=False; n_ok=n_fail=n_skip=0
    for pos,(_,row) in enumerate(df.iterrows()):
        t=str(row[tok]).lower()
        if t not in look.index or pd.isna(row.get(dirc)): n_fail+=1; continue
        if args.gate>0:
            erm=pd.to_numeric(row.get(erc),errors="coerce")
            if not (pd.notna(erm) and erm<args.gate): n_skip+=1; continue
        dr=look.loc[t]
        if epic is not None: epi=parse_index_list(dr[epic])
        else:
            ws=int(row[wsc]) if (wsc and pd.notna(row.get(wsc))) else 0
            epi=[ws+i for i in parse_index_list(dr[epim])]
        npdb=native_index.get(str(row[idc]).lower()) if idc else None
        res=process(str(row[dirc]),npdb,epi,args.exclude_dist)
        if res is None: n_fail+=1; continue
        plain[pos],natin[pos],aware[pos]=res["n_plain"],res["native_in"],res["n_aware"]
        n_ok+=1
        want=(args.viz_token and t==args.viz_token.lower()) or (args.viz_first and not did_viz and npdb is not None and res["native_in"]>=0)
        if want and not did_viz:
            make_viz(res["epi_ca"],res["native_heavy"],args.exclude_dist,
                     res["base"],res["normal"],res["native_in_cyl"],
                     res["n_plain"],res["n_aware"],args.viz_out)
            did_viz=True; print(f"  wrote viz for {t} -> {args.viz_out}")
        if (pos+1)%100==0: print(f"  {pos+1}  ok={n_ok} fail={n_fail}")

    df["cylinder_ca_clashes"]=plain
    df["native_in_cylinder"]=natin
    df["cylinder_native_aware"]=aware
    out=Path(args.out_csv); out.parent.mkdir(parents=True,exist_ok=True)
    df.to_csv(out,index=False)
    print(f"\nWrote {out}  (ok={n_ok} fail={n_fail} gate-skipped={n_skip})")

    nv=pd.to_numeric(df["native_in_cylinder"],errors="coerce").dropna()
    if len(nv):
        print(f"\nnative_in_cylinder over {len(nv)} designs: "
              f"mean {nv.mean():.1f}  median {nv.median():.0f}  "
              f"frac>0 {(nv>0).mean():.2f}  max {nv.max():.0f}")
    tc=first_present(df.columns,TRUE_CLASH)
    if tc is not None:
        truth=pd.to_numeric(df[tc],errors="coerce")
        print(f"\nCorrelation vs true clash ({tc}):")
        for c in ("cylinder_ca_clashes","cylinder_native_aware"):
            p=pd.to_numeric(df[c],errors="coerce"); m=p.notna()&truth.notna()
            if m.sum()>=10:
                print(f"  {c:24s} Pearson {p[m].corr(truth[m]):+.3f}  "
                      f"Spearman {p[m].corr(truth[m],method='spearman'):+.3f}")


if __name__=="__main__":
    main()
