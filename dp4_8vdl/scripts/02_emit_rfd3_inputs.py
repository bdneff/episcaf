#!/usr/bin/env python3
"""
02_emit_rfd3_inputs.py -- one RFD3 input JSON per contig for an 8VDL scaffolding run.

Reads a contigs CSV from 01_generate_contigs.py (multi-island aware: fixes every residue listed in
`fixed_resids`, so both the epitope20 single island and the hotspots two-island motif work). Emits
the RFD3 JSON spec + an inputs_manifest.csv that 03_rfd3_array.sbatch consumes.

Usage:
  python scripts/02_emit_rfd3_inputs.py --contigs_csv 01_contigs/epitope20.csv --out_dir 02_rfd3/epitope20/inputs
"""
from __future__ import annotations
import argparse
import csv
import json
from pathlib import Path

import pandas as pd


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--contigs_csv", required=True)
    ap.add_argument("--out_dir", required=True)
    args = ap.parse_args()

    out_dir = Path(args.out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(args.contigs_csv)
    manifest_path = out_dir.parent / "inputs_manifest.csv"

    n = 0
    with open(manifest_path, "w", newline="") as mf:
        w = csv.writer(mf)
        w.writerow(["design_id", "json_path"])
        for _, row in df.iterrows():
            target = str(row["target"])
            chain = str(row["chain"])
            resids = [int(x) for x in str(row["fixed_resids"]).split(",")]
            design_id = f"{target}_contig{int(row['contig_id']):04d}"
            fixed_atoms = {f"{chain}{r}": "BKBN" for r in resids}
            spec = {
                "input":                   str(Path(row["input_pdb"]).resolve()),  # portable: resolve here
                "contig":                  str(row["contig_string"]),
                "length":                  f"{int(row['total_len'])}-{int(row['total_len'])}",
                "select_fixed_atoms":      fixed_atoms,
                "select_unfixed_sequence": False,
            }
            json_path = out_dir / f"{design_id}.json"
            json_path.write_text(json.dumps({design_id: spec}, indent=2))
            w.writerow([design_id, str(json_path)])
            n += 1

    print(f"wrote {n:,} RFD3 input JSONs -> {out_dir}")
    print(f"wrote manifest -> {manifest_path}")
    print(f"\nsubmit RFD3 (one task per contig, 8 backbones each):")
    print(f"  sbatch --array=1-{n}%200 scripts/03_rfd3_array.sbatch {out_dir.parent}")


if __name__ == "__main__":
    main()
