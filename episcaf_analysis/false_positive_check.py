#!/usr/bin/env python3
"""
false_positive_check.py

Direct test of John's concern. A false positive is a design with little/no TRUE
antibody clash that the cylinder nonetheless flags. The 1 A native-volume carve should
shrink that population. This reports, per true-clash band, how many designs the
cylinder flags before (plain) vs after (native-aware) -- and, as a selectivity check,
whether it leaves the genuinely-clashing designs mostly alone.

Run on the add_native_cylinder output, so plain (cylinder_ca_clashes) and carved
(cylinder_native_aware) come from the SAME cylinder and differ only by the carve.

    python scripts/false_positive_check.py \
        --native runs/run_rfd3_mpnn/04_filter/metrics_native_cyl.csv
"""
import argparse
import numpy as np
import pandas as pd


def report(name, plain, aware):
    n = len(plain)
    dm = (aware.mean() - plain.mean()) / plain.mean() * 100 if plain.mean() else 0
    print(f"\n{name}   (n={n})")
    print(f"   mean cylinder:   {plain.mean():5.1f}  ->  {aware.mean():5.1f}   ({dm:+.0f}%)")
    print(f"   {'flagged at':<14}{'before':>9}{'after':>9}{'removed':>10}")
    for thr in (1, 5, 10):
        b = int((plain >= thr).sum()); a = int((aware >= thr).sum())
        pct = (b - a) / b * 100 if b else 0
        print(f"   cylinder >= {thr:<3d}{b:>9d}{a:>9d}{b-a:>7d} ({pct:+.0f}%)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--native", required=True, help="add_native_cylinder output csv")
    args = ap.parse_args()

    d = pd.read_csv(args.native, low_memory=False)
    tc = pd.to_numeric(d["af3_n_clash_res"], errors="coerce")
    plain = pd.to_numeric(d["cylinder_ca_clashes"], errors="coerce")
    aware = pd.to_numeric(d["cylinder_native_aware"], errors="coerce").fillna(plain)
    ok = tc.notna() & plain.notna() & aware.notna()
    tc, plain, aware = tc[ok], plain[ok], aware[ok]
    print(f"designs analyzed: {len(tc)}")

    # false-positive bands (low/no true clash) -- where the carve should help
    report("true clash == 0  (false positives)", plain[tc == 0], aware[tc == 0])
    report("true clash <= 3", plain[tc <= 3], aware[tc <= 3])
    # selectivity check: genuinely clashing designs -- carve should barely touch these
    report("true clash >= 10 (real clashes)", plain[tc >= 10], aware[tc >= 10])

    carved = aware < plain
    if carved.any():
        print(f"\noverall: {int(carved.sum())}/{len(plain)} designs carved; "
              f"mean residues removed where carved: {(plain - aware)[carved].mean():.1f}")
    # selectivity ratio: fractional reduction on clash-free vs on real-clash designs
    def frac_drop(mask):
        p = plain[mask]; a = aware[mask]
        return 1 - a.sum() / p.sum() if p.sum() else 0
    print(f"fractional cylinder reduction:  clash-free {100*frac_drop(tc==0):.0f}%   "
          f"vs  real-clash(>=10) {100*frac_drop(tc>=10):.0f}%")
    print("(bigger gap = the carve is selectively removing false positives, "
          "not real signal)")


if __name__ == "__main__":
    main()
