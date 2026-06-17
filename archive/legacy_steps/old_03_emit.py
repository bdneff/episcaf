#!/usr/bin/env python3
"""
Emit one RFD3 DesignInputSpecification JSON per row in contigs.parquet.

Assumptions (validated on your data):
- parquet column `id` is the ABDB complex stem, and cleaned PDB is {id}.pdb
- `contig_string` and `contig_length` are already correct and should be passed through
"""

import argparse
import json
import os
from pathlib import Path

import pandas as pd


def dump_json(path: Path, obj: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2)
        f.write("\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs_parquet", required=True)
    ap.add_argument("--out_dir", required=True, help="e.g. runs/.../02_rfd3/inputs")
    ap.add_argument("--cleaned_pdb_dir", required=True)
    ap.add_argument("--limit", type=int, default=0, help="If >0, only emit first N JSONs")
    ap.add_argument("--dump_trajectory", action="store_true")
    args = ap.parse_args()

    df = pd.read_parquet(args.contigs_parquet)

    required = ["design_id", "id", "contig_string", "contig_length", "seed"]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Missing required columns in contigs parquet: {missing}")

    out_dir = Path(args.out_dir)
    cleaned = Path(args.cleaned_pdb_dir)

    if args.limit and args.limit > 0:
        df = df.head(args.limit).copy()

    manifest = []
    for _, r in df.iterrows():
        abdb_id = r["id"]
        input_pdb = cleaned / f"{abdb_id}.pdb"
        if not input_pdb.exists():
            raise FileNotFoundError(f"Missing cleaned PDB for id={abdb_id}: {input_pdb}")

        design_id = r["design_id"]
        js = {
            "design_id": design_id,
            "input_structure": str(input_pdb),
            "contig": r["contig_string"].replace("/", ","),
            "length": r["contig_length"],
            "dump_trajectory": bool(args.dump_trajectory),
            "seed": int(r["seed"]),
            # keep this explicit (you've been using it successfully)
            "prevalidate_inputs": True,
        }

        out_path = out_dir / f"{design_id}.json"
        dump_json(out_path, js)

        manifest.append({
            "design_id": design_id,
            "json": str(out_path),
            "input_structure": str(input_pdb),
            "contig": r["contig_string"].replace("/",","),
            "length": r["contig_length"],
            "seed": int(r["seed"]),
        })

    man_path = out_dir.parent / "inputs_manifest.csv"
    pd.DataFrame(manifest).to_csv(man_path, index=False)

    print(f"Wrote {len(manifest)} JSONs to: {out_dir}")
    print(f"Wrote manifest: {man_path}")


if __name__ == "__main__":
    main()
