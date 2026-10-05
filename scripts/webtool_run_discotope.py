#!/usr/bin/env python3
"""Run DiscoTope-3.0 on the four web-tool antigens via BioLib and bake the scores into the viewer.

The ESM-IF1 compute runs on BioLib's cloud (anonymous, no local torch/GPU), so the only dependency
here is the lightweight client:  pip install pybiolib. It submits each antigen-only PDB, downloads the
per-residue result, then calls scripts/webtool_discotope_bake.py to write DiscoTope-3.0 scores into the
B-factor of docs/pdb/<id>_scored.pdb (what predict mode colors). One command, fully reproducible.

Run (from repo root, in a python with pybiolib):  python scripts/webtool_run_discotope.py
"""
import os, sys, glob, subprocess
import biolib

TARGETS = {"8PWW": "A", "7OX3": "C", "5FHX": "A", "8VDL": "C"}

def main():
    biolib.utils.STREAM_STDOUT = True
    app = biolib.load("DTU/DiscoTope-3")
    pairs = []
    for t in TARGETS:
        print(f"=== DiscoTope-3 on {t} ===", flush=True)
        job = app.cli(args=f"-f docs/pdb/{t}_scored.pdb --struc_type solved --out_dir output")
        out = f"/tmp/dt_out/{t}"; os.makedirs(out, exist_ok=True); job.save_files(out)
        csv = glob.glob(out + "/**/*discotope3.csv", recursive=True)[0]
        pairs.append(f"{t}={csv}"); print(f"  -> {csv}", flush=True)
    subprocess.run([sys.executable, "scripts/webtool_discotope_bake.py", *pairs,
                    "--resid-col", "res_id", "--score-col", "DiscoTope-3.0_score"], check=True)

if __name__ == "__main__":
    main()
