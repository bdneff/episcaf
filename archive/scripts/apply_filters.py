#!/usr/bin/env python3
import pandas as pd

df = pd.read_csv("runs/run_test_rfd3_nompmn/04_filter/metrics_rfd3_af3_with_clash_LIMIT200.csv")

mask = (
    df["rmsd_bb_epitope"].notna() &
    (df["rmsd_bb_epitope"] < 2.0) &
    (df["rmsd_bb_all"] < 1.0) &
    (df["pae_mean_all"] < 5.0) &
    (df["af3_n_clash_res"] <= 0)
)

passed = df[mask].copy()

print("Total rows:", len(df))
print("Passing rows:", len(passed))

if len(passed) > 0:
    print("\nTop passing rows:")
    print(
        passed.sort_values("pae_mean_all")
        [["id","pred","rmsd_bb_all","rmsd_bb_epitope","pae_mean_all","af3_n_clash_res"]]
        .head(100)
        .to_string(index=False)
    )
