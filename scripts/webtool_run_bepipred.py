#!/usr/bin/env python3
"""Run BepiPred-3.0 on the four web-tool antigens via BioLib and report the output layout.

Companion to webtool_run_discotope.py. BepiPred-3.0 is *sequence*-based: it takes a FASTA of the
antigen chains and returns a per-residue B-cell-epitope score, so its scores are mapped onto the
served antigen PDB by residue ORDER (sequence position i -> the i-th CA residue of the chain), via
webtool_discotope_bake.py --by-order. This differs from DiscoTope-3.0, which is structure-based and
maps by res_id.

Every biolib app.cli() call submits a cloud container (anonymous, no local torch), so a run takes a
few minutes. This script submits once on the combined FASTA, saves all output under /tmp/bp_out, and
prints the file tree + the head of any CSV so the exact score/accession columns can be read before
baking.

Run (from repo root, in a python with pybiolib):  python scripts/webtool_run_bepipred.py
"""
import os, glob, sys
import biolib

FASTA = "/tmp/ag_fasta/all.fasta"
OUT = "/tmp/bp_out"

def main():
    biolib.utils.STREAM_STDOUT = True
    if not os.path.exists(FASTA):
        sys.exit(f"missing {FASTA}; build it first (antigen chains, one record per target)")
    app = biolib.load("DTU/BepiPred-3")
    print("app:", app, flush=True)
    job = app.cli(args=f"-i {FASTA} -o output -pred mjv_pred")
    os.makedirs(OUT, exist_ok=True)
    job.save_files(OUT)
    files = sorted(glob.glob(OUT + "/**/*", recursive=True))
    print("=== FILES ===", flush=True)
    for f in files:
        print(" ", f, flush=True)
    for c in [f for f in files if f.lower().endswith(".csv")]:
        print(f"=== HEAD {c} ===", flush=True)
        with open(c) as fh:
            for i, line in enumerate(fh):
                if i >= 8:
                    break
                print("   " + line.rstrip(), flush=True)

if __name__ == "__main__":
    main()
