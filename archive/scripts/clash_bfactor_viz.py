#!/usr/bin/env python3
"""
clash_bfactor_viz.py

Write per-residue clash frequencies as B-factors into the true complex PDB
for visualization in VMD.

B-factor encoding:
  - Antigen (chain A), non-epitope residues : 0.0
  - Antigen (chain A), epitope residues     : 50.0  (fixed, for highlighting)
  - Antibody (chains B+C), each residue     : (clash_count / n_designs) * 100

Note: af3_clash_resindices are indices into the antibody selection (segid B or C),
not the full universe. Index 0 = first antibody residue.

Usage:
    python scripts/clash_bfactor_viz.py \
        --metrics_csv  runs/run_rfd3_mpnn/04_filter/metrics_partial.csv \
        --true_pdb     /tgen_labs/altin/.../cleaned/4xwo_5P.pdb \
        --epitope_id   4xwo_5P \
        --dp2_parquet  datasets/dp2.parquet \
        --out_pdb      runs/run_rfd3_mpnn/clash_density/4xwo_5P/4xwo_5P_clash_bfactor.pdb
"""

import argparse
import ast
from collections import Counter
from pathlib import Path

import MDAnalysis as mda
import numpy as np
import pandas as pd


def parse_resindices(x) -> list:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return []
    if isinstance(x, (list, np.ndarray)):
        return [int(i) for i in x]
    s = str(x).strip()
    if not s or s in ("[]", "nan"):
        return []
    try:
        return [int(i) for i in ast.literal_eval(s)]
    except Exception:
        return []


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metrics_csv",  required=True)
    parser.add_argument("--true_pdb",     required=True)
    parser.add_argument("--epitope_id",   required=True)
    parser.add_argument("--dp2_parquet",  required=True)
    parser.add_argument("--out_pdb",      required=True)
    args = parser.parse_args()

    # --- Load data ---
    print(f"Loading metrics for {args.epitope_id} ...")
    df = pd.read_csv(args.metrics_csv, low_memory=False)
    sub = df[df["id"] == args.epitope_id].copy()
    n_designs = len(sub)
    print(f"  {n_designs} designs")

    print(f"Loading dp2 for epitope residue indices ...")
    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2_row = dp2[dp2["id"] == args.epitope_id].iloc[0]
    true_epi_ris = set(int(i) for i in dp2_row["epitope_chunk_resindices"])
    print(f"  {len(true_epi_ris)} epitope residues (indices into full universe)")

    print(f"Loading true complex PDB ...")
    u = mda.Universe(args.true_pdb)
    print(f"  {len(u.residues)} total residues, chains: {set(u.atoms.segids)}")

    # antibody residues as a separate selection — clash indices are into this
    ab_residues = u.select_atoms("segid B or segid C").residues
    print(f"  {len(ab_residues)} antibody residues")

    # --- Count clash frequency per antibody residue ---
    print(f"Counting clash frequencies ...")
    clash_counter = Counter()
    n_with_clash = 0

    for _, row in sub.iterrows():
        ris = parse_resindices(row["af3_clash_resindices"])
        if ris:
            n_with_clash += 1
            for ri in ris:
                clash_counter[ri] += 1

    print(f"  {n_with_clash} / {n_designs} designs had at least one clash")
    print(f"  {len(clash_counter)} unique antibody residues ever clashed")

    if clash_counter:
        print(f"  Top 5 most frequent clash residues:")
        for ri, count in clash_counter.most_common(5):
            res = ab_residues[ri]
            print(f"    ab_resindex={ri}  resname={res.resname}  resid={res.resid}"
                  f"  segid={res.segid}  count={count} ({100*count/n_designs:.1f}%)")

    # --- Assign B-factors ---
    # reset everything to 0
    u.atoms.tempfactors = 0.0

    # epitope residues on antigen: fixed value of 50
    for ri in true_epi_ris:
        if ri < len(u.residues):
            u.residues[ri].atoms.tempfactors = 50.0

    # antibody residues: clash frequency as % of designs
    max_count = max(clash_counter.values()) if clash_counter else 1
    for ri, count in clash_counter.items():
        if ri < len(ab_residues):
            bfactor = (count / n_designs) * 100
            ab_residues[ri].atoms.tempfactors = bfactor

    # --- Write output PDB ---
    out_pdb = Path(args.out_pdb)
    out_pdb.parent.mkdir(parents=True, exist_ok=True)
    u.atoms.write(str(out_pdb))
    print(f"\nWrote: {out_pdb}")

    max_rate = 100 * max_count / n_designs
    print(f"""
=== VMD VISUALIZATION INSTRUCTIONS ===

vmd -e clash_viz.tcl -args {out_pdb.name}

B-factor scale:
  Antigen non-epitope : 0.0
  Epitope patch       : 50.0
  Antibody            : 0.0 (no clash) to {max_rate:.1f} (max clash rate %)

=== STATS ===
  Total designs:              {n_designs}
  Designs with clash:         {n_with_clash} ({100*n_with_clash/n_designs:.1f}%)
  Epitope residues:           {len(true_epi_ris)}
  Antibody residues clashed:  {len(clash_counter)} / {len(ab_residues)}
  Max clash rate:             {max_rate:.1f}%
""")


if __name__ == "__main__":
    main()
