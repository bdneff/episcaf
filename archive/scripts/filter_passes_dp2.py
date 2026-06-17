#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd

REQ = ["epitope_chunk_rmsd_vs_mpnn", "rmsd_ca_all", "mean_pae", "af3_n_clash_res"]

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True, help="metrics_by_pred_with_clash.csv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--epi_rmsd_max", type=float, default=1.0)
    ap.add_argument("--rmsd_ca_all_max", type=float, default=2.0)
    ap.add_argument("--mean_pae_max", type=float, default=5.0)
    ap.add_argument("--clash_max", type=int, default=0)
    ap.add_argument("--top", type=int, default=25)
    args = ap.parse_args()

    df = pd.read_csv(args.csv)

    missing = [c for c in REQ if c not in df.columns]
    if missing:
        raise SystemExit(f"Missing required columns: {missing}\nHave: {list(df.columns)}")

    finished = df[df[REQ].notna().all(axis=1)].copy()

    passed = finished[
        (finished["epitope_chunk_rmsd_vs_mpnn"] < args.epi_rmsd_max) &
        (finished["rmsd_ca_all"] < args.rmsd_ca_all_max) &
        (finished["mean_pae"] < args.mean_pae_max) &
        (finished["af3_n_clash_res"] <= args.clash_max)
    ].copy()

    passed = passed.sort_values(["mean_pae", "overall_rmsd", "epitope_chunk_rmsd_vs_mpnn"], ascending=True)

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    passed.to_csv(out, index=False)

    print(f"Total rows: {len(df):,}")
    print(f"Finished rows (have all filter metrics): {len(finished):,}")
    print(f"Passed: {len(passed):,} / {len(finished):,} ({(100*len(passed)/len(finished) if len(finished) else 0):.4f}%)")
    print(f"Wrote: {out}")

    if args.top > 0:
        cols = [c for c in [
            "id","assay_scaffolded_epitope_id","pred",
            "epitope_chunk_rmsd_vs_mpnn","overall_rmsd","mean_pae",
            "af3_n_clash_res","pae_mean_all","rmsd_ca_all","ptm"
        ] if c in passed.columns]
        print("\nTop passes:")
        print(passed[cols].head(args.top).to_string(index=False))
if __name__ == "__main__":
    main()
