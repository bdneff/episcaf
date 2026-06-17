#!/usr/bin/env python3
"""
composite_swap_validate.py

Validation SWAP: does using the native-aware cylinder (cylinder_native_aware) as the
composite's clash term change recovery of true passers, vs the plain cylinder
(cylinder_ca_clashes)?

It reuses all_scored.csv (the gated, scored set -- already has the four non-clash
percentile score_ columns, score_cylinder, is_pass, and composite) plus the
cylinder_native_aware column from add_native_cylinder.py. The composite is rebuilt by
swapping ONLY the clash term, so the other four terms are byte-identical to your
pipeline. The plain-clash run therefore reproduces your composite exactly -- it should
land on your known baseline (per-epitope recall 0.359, coverage 10/13, precision
0.080), which is the sanity check that the rebuild is faithful.

    python scripts/composite_swap_validate.py \
        --scored runs/run_rfd3_mpnn/04_filter/all_scored.csv \
        --native runs/run_rfd3_mpnn/04_filter/metrics_native_cyl.csv \
        --w_cylinder 0.15 --topk 15
"""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd


def oriented_pct(x, method):
    """Lower value -> higher percentile (lower clash is better)."""
    return 1.0 - pd.Series(np.asarray(x, float)).rank(pct=True, method=method).to_numpy()


def select_and_validate(df, comp, topk):
    """Dedupe by token (keep best composite), top-k per id, validate vs is_pass."""
    d = df.assign(_c=comp)
    sel = (d.sort_values("_c", ascending=False)
             .drop_duplicates("token", keep="first")
             .sort_values(["id", "_c"], ascending=[True, False]))
    sel["_r"] = sel.groupby("id").cumcount() + 1
    sel = sel[sel["_r"] <= topk]
    truepass = d[d["is_pass"] == 1]
    idx = set(sel.index)
    tp = int(truepass.index.isin(idx).sum())
    grp_all = set(truepass["id"].astype(str))
    grp_kept = set(sel.loc[sel.index.isin(truepass.index), "id"].astype(str))
    per = [gg.index.isin(idx).sum() / min(len(gg), topk)
           for _, gg in truepass.groupby(truepass["id"].astype(str))]
    return dict(per_ep=float(np.mean(per)) if per else float("nan"),
                coverage=f"{len(grp_kept)}/{len(grp_all)}",
                precision=tp / len(sel) if len(sel) else float("nan"),
                recall=tp / len(truepass) if len(truepass) else float("nan"),
                shortlist=len(sel))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scored", required=True, help="all_scored.csv")
    ap.add_argument("--native", required=True, help="add_native_cylinder output csv")
    ap.add_argument("--w_cylinder", type=float, default=0.15)
    ap.add_argument("--topk", type=int, default=15)
    args = ap.parse_args()

    d = pd.read_csv(args.scored, low_memory=False)
    need = ("token", "id", "is_pass", "composite", "score_cylinder",
            "cylinder_ca_clashes", "af3_dir")
    miss = [c for c in need if c not in d.columns]
    if miss:
        raise SystemExit(f"all_scored is missing columns: {miss}")
    d["_base"] = d["af3_dir"].astype(str).map(lambda p: Path(p).name)

    nat = pd.read_csv(args.native, low_memory=False)
    if "cylinder_native_aware" not in nat.columns or "af3_dir" not in nat.columns:
        raise SystemExit("native file needs af3_dir and cylinder_native_aware columns")
    nat["_base"] = nat["af3_dir"].astype(str).map(lambda p: Path(p).name)
    na = (nat.dropna(subset=["cylinder_native_aware"])
             .drop_duplicates("_base").set_index("_base")["cylinder_native_aware"])
    d["cylinder_native_aware"] = d["_base"].map(na)
    n_have = int(d["cylinder_native_aware"].notna().sum())
    # graceful fallback: missing native-aware -> plain cylinder
    d["native_aware_filled"] = d["cylinder_native_aware"].fillna(d["cylinder_ca_clashes"])
    plain_v = pd.to_numeric(d["cylinder_ca_clashes"], errors="coerce").to_numpy()
    aware_v = pd.to_numeric(d["native_aware_filled"], errors="coerce").to_numpy()
    n_carved = int(np.nansum(aware_v < plain_v))
    print(f"merged native-aware for {n_have}/{len(d)} scored rows; "
          f"{n_carved} rows carved to a lower value")

    # --- calibrate the oriented-percentile method against stored score_cylinder ----
    target = (pd.to_numeric(d["score_cylinder"], errors="coerce")
              / args.w_cylinder).to_numpy()
    best = None
    for meth in ("average", "min", "max", "dense"):
        err = np.nanmax(np.abs(oriented_pct(plain_v, meth) - target))
        if best is None or err < best[1]:
            best = (meth, err)
    method, err = best
    print(f"percentile calibration: method='{method}'  "
          f"max|err| vs score_cylinder/w = {err:.3g}")
    if err > 1e-3:
        print("  WARNING: calibration imperfect; native-aware composite is approximate")

    # --- composites: plain = stored composite; aware = swap only the clash term ----
    comp_plain = pd.to_numeric(d["composite"], errors="coerce").to_numpy()
    score_cyl = pd.to_numeric(d["score_cylinder"], errors="coerce").to_numpy()
    comp_aware = comp_plain - score_cyl + args.w_cylinder * oriented_pct(aware_v, method)

    r_plain = select_and_validate(d, comp_plain, args.topk)
    r_aware = select_and_validate(d, comp_aware, args.topk)

    print("\n=== composite validation: plain cylinder vs native-aware (clash-term swap) ===")
    print(f"{'metric':<26}{'plain':>14}{'native-aware':>16}")
    for k, label in [("per_ep", "per-epitope recall"), ("coverage", "coverage"),
                     ("precision", "precision"), ("recall", "raw recall"),
                     ("shortlist", "shortlist size")]:
        pv, av = r_plain[k], r_aware[k]
        if isinstance(pv, float):
            print(f"{label:<26}{pv:>14.3f}{av:>16.3f}")
        else:
            print(f"{label:<26}{str(pv):>14}{str(av):>16}")
    print("\n(plain column should match your known baseline: 0.359 / 10-13 / 0.080)")


if __name__ == "__main__":
    main()
