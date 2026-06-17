#!/usr/bin/env python3
from pathlib import Path
import argparse
import gemmi

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_root", required=True, help="AF3 outputs dir, e.g. runs/.../03_af3/outputs")
    ap.add_argument("--pattern", default="*_model.cif", help="Glob pattern for AF3 model cif files")
    ap.add_argument("--overwrite", action="store_true")
    args = ap.parse_args()

    root = Path(args.in_root)
    n_in = 0
    n_written = 0
    n_skip = 0

    for cif in root.rglob(args.pattern):
        n_in += 1
        pdb = cif.with_suffix(".pdb")
        if pdb.exists() and not args.overwrite:
            n_skip += 1
            continue
        try:
            st = gemmi.read_structure(str(cif))
            st.write_pdb(str(pdb))
            n_written += 1
        except Exception as e:
            print("FAIL:", cif, "->", e)

    print(f"Found {n_in} CIFs, wrote {n_written} PDBs, skipped {n_skip}")

if __name__ == "__main__":
    main()
