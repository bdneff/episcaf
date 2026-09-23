#!/usr/bin/env python3
"""Copy a LOW-RMSD and HIGH-RMSD C1 design PDB for one epitope, for local rendering (render_design.tcl).

Run ON GEMINI (the design PDBs live on the cluster), from the repo root, after `git pull`:
    conda activate ~/rfd3/env/rfd3_py312     # or any env with pandas
    python scripts/fetch_design_pair.py 8pww

It picks the pair from the local superset (by epitope RMSD), then finds each design's AF3 output PDB
under the C1 run and copies both into data/dp4_binding/designs/. Then `git add` + push (or just scp the
two PDBs back to your laptop) and I render them. `find` is scoped to the C1 run dir, so it does not walk
all of /tgen_labs.
"""
import sys, re, shutil, subprocess
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from configs.paths import DATA
RUN = DATA["antibody_runs"] / "whole_epitope_rfd3"          # the C1 run tree (on $WS)
epi = (sys.argv[1] if len(sys.argv) > 1 else "8pww").lower()
OUT = ROOT / "data/dp4_binding/designs"; OUT.mkdir(parents=True, exist_ok=True)

sup = ROOT / "data/libraries/dp4_superset.csv"
if not sup.exists():
    sup = ROOT / "data/libraries/dp4_superset.csv.gz"   # the committed artifact on the cluster
s = pd.read_csv(sup, low_memory=False)
s = s[s.component == "C1"].copy()
s["r"] = pd.to_numeric(s.epitope_rmsd, errors="coerce")
s["epi"] = s.target.astype(str).str.split("_").str[0].str.lower()
sub = s[(s.epi == epi) & s.r.notna()]
if sub.empty:
    sys.exit(f"no C1 designs for epitope '{epi}'")
lo = sub.loc[sub.r.idxmin()]
hiset = sub[(sub.r > 3) & (sub.r < 6)].sort_values("r")
hi = hiset.iloc[0] if len(hiset) else sub.loc[sub.r.idxmax()]

def fetch(row, tag):
    pid = str(row.predID)
    m = re.search(r"(contig\d+).*?(dldesign_\d+)", pid)
    dl = m.group(2) if m else pid.split("_")[-1]
    contig = m.group(1) if m else ""
    # fast name search scoped to the run dir; filter to this contig + a real coordinate PDB
    try:
        found = subprocess.run(["find", str(RUN), "-name", f"*{dl}*.pdb"],
                               capture_output=True, text=True, timeout=600).stdout.split()
    except Exception as e:
        found = []; print(f"  find failed: {e}")
    found = [f for f in found if contig in f]
    print(f"{tag}: rmsd={row.r:.2f}  match tokens=({contig},{dl})  candidates={len(found)}")
    for f in found[:6]:
        print("   ", f)
    if found:
        dst = OUT / f"design_{epi}_{tag}_rmsd{row.r:.2f}.pdb"
        shutil.copy(found[0], dst); print("   copied ->", dst.relative_to(ROOT))
    else:
        print("   no PDB found automatically; predID:", pid)

print(f"epitope {epi}: low rmsd {lo.r:.2f}, high rmsd {hi.r:.2f}\n")
fetch(lo, "low"); fetch(hi, "high")
print(f"\nDone. Copy data/dp4_binding/designs/*.pdb back to the laptop (or commit) for rendering.")
