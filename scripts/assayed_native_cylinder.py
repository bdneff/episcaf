#!/usr/bin/env python3
"""Compute the native-aware cylinder for the ASSAYED scaffold designs (DP3 binding set).

Runs on Gemini. The carve is IDENTICAL to scripts/dp3_native_cylinder.py -- this script imports
process()/load_af3()/parse_index_list() from it, so the geometry cannot drift. The only new
piece is locating each assayed design's AF3 structure, since the assayed (dp2.parquet) designs
are a different molecule set than the RFD3 run and dp2.parquet stores no path.

Inputs:
  --worklist   results/assayed_cylinder_worklist.csv (from build_assayed_cylinder_worklist.py):
               one row per assayed scaffold with designedSequence, id, assay_scaffolded_epitope_id.
  --dp2_parquet  dp2.parquet -- epitope residue indices, looked up by assay_scaffolded_epitope_id.
  --native_dir   AbDb cleaned complex dir (ABDB_CLEANED_PDB_DIR); antigen = chain A of <id>.pdb.
  --af3_root     root under which the assayed designs' AF3 outputs live, one dir per design
                 named by the assay hash: <af3_root>/<assay_scaffolded_epitope_id>/<hash>_model.cif.gz
                 (verified 403/403 under sourced_antibody_v1/no_antibody/af3_predictions). The dir
                 is resolved directly; --key_col overrides which worklist column names the dir.
  --out_csv      output; one row per design with cylinder_ca_clashes, native_in_cylinder,
                 cylinder_native_aware, keyed by designedSequence (joins to dp3_binding_metrics.csv).

Example (defaults in scripts/assayed_native_cylinder.sbatch already point at the right folder):
  conda activate ~/rfd3/env/rfd3_py312
  B=/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/sourced_antibody_v1/no_antibody
  python3 scripts/assayed_native_cylinder.py \
      --worklist    results/assayed_cylinder_worklist.csv \
      --dp2_parquet $B/assay_scaffold_simple_metrics_403.parquet \
      --native_dir  /tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned \
      --af3_root    $B/af3_predictions \
      --out_csv     results/assayed_native_cyl.csv \
      --exclude_dist 1.0
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
# identical carve + helpers; do NOT re-implement
from dp3_native_cylinder import process, parse_index_list, first_present, DP2_EPI_COLS  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--worklist", required=True)
    ap.add_argument("--dp2_parquet", required=True)
    ap.add_argument("--native_dir", required=True)
    ap.add_argument("--af3_root", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--key_col", default="assay_scaffolded_epitope_id",
                    help="worklist column whose value identifies each AF3 output under --af3_root")
    ap.add_argument("--exclude_dist", type=float, default=1.0)
    args = ap.parse_args()

    wl = pd.read_csv(args.worklist)
    wl[args.key_col] = wl[args.key_col].astype(str).str.lower()

    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"] = dp2.assay_scaffolded_epitope_id.astype(str).str.lower()
    epic = first_present(dp2.columns, DP2_EPI_COLS)
    epi_lookup = (dp2.drop_duplicates("assay_scaffolded_epitope_id")
                     .set_index("assay_scaffolded_epitope_id")[epic])

    native_index = {p.stem.lower(): p for p in Path(args.native_dir).glob("*.pdb")}
    missing_native = sorted(set(wl.id.str.lower()) - set(native_index))
    if missing_native:
        print(f"WARNING: no native PDB for ids: {missing_native}")

    # The assayed-design AF3 outputs are stored one dir per design, named exactly by the
    # assay hash: <af3_root>/<key>/<key>_model.cif.gz (verified 403/403, sourced_antibody_v1/
    # no_antibody). So resolve the dir directly -- no tree walk. find_af3_cif (in the shared
    # carve) handles the .cif.gz.
    af3_root = Path(args.af3_root)

    plain = np.full(len(wl), np.nan); natin = np.full(len(wl), np.nan)
    aware = np.full(len(wl), np.nan); af3_found = []
    n_ok = n_nodir = n_fail = 0
    for pos, row in enumerate(wl.itertuples(index=False)):
        h = getattr(row, "assay_scaffolded_epitope_id")
        key = str(getattr(row, args.key_col)).lower()
        d = af3_root / key
        af3 = d if d.is_dir() else None
        af3_found.append(str(af3) if af3 else "")
        if af3 is None:
            n_nodir += 1; continue
        epi = parse_index_list(epi_lookup.get(h))
        npdb = native_index.get(str(row.id).lower())
        res = process(str(af3), npdb, epi, args.exclude_dist)
        if res is None:
            n_fail += 1; continue
        plain[pos], natin[pos], aware[pos] = res
        n_ok += 1
        if (pos + 1) % 50 == 0:
            print(f"  {pos+1}/{len(wl)} ok={n_ok} no_af3_dir={n_nodir} fail={n_fail}", flush=True)

    out = wl.copy()
    out["af3_dir"] = af3_found
    out["cylinder_ca_clashes"] = plain
    out["native_in_cylinder"] = natin
    out["cylinder_native_aware"] = aware
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out_csv, index=False)
    print(f"\nwrote {args.out_csv}  ok={n_ok} no_af3_dir={n_nodir} fail={n_fail} of {len(wl)}")
    if n_nodir:
        ex = out.loc[out.cylinder_native_aware.isna(), args.key_col].head(3).tolist()
        print(f"  no AF3 dir found for {n_nodir} (e.g. {ex}). If 0 found, the assayed outputs are "
              f"keyed differently -- check --af3_root layout and set --key_col accordingly.")


if __name__ == "__main__":
    main()
