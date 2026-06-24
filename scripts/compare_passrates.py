#!/usr/bin/env python3
"""
compare_passrates.py -- per-filter four-filter pass rates across the three design sets, to
disentangle the backbone generator (RFD1 vs RFD3) from the epitope structure (Lawson's
mostly two-island DP3 vs this run's single islands). Answers John's question (2026-06-24):
how much of the per-island run's different rates is RFD1-vs-RFD3 vs two-island-vs-single?

Three sets (same 59 mAb epitopes throughout):
  1. RFD1+MPNN on Lawson's DP3 contigs        dp2.parquet            (two-island dominated)
  2. RFD3+MPNN on the SAME DP3 epitopes        metrics_cylinder_full.csv  (isolates the generator)
  3. RFD3+MPNN on single islands (THIS run)    metrics_dual_island.parquet (isolates island count)

Each cell is the percentage of that set's VALID designs (all four filter metrics present)
that clear the filter; "valid" counts a design only if it has an AF3 result and a computable
clash. Run wherever the files live (sets 1-2 are local off-cluster; set 3 lives on the cluster
run dir, so pass --perisland there).

Usage:
  python3 scripts/compare_passrates.py \
      --rfd1 ../known_antigen/analysis/full_run/dp2.parquet \
      --rfd3 ../known_antigen/analysis/data/metrics_cylinder_full.csv \
      --perisland runs/dual_island_rfd3/05_analysis/metrics_dual_island.parquet
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import pandas as pd

# (dataset column names) per source; all map to the same four filters.
SCHEMAS = {
    "rfd1": dict(epi="epitope_chunk_rmsd_vs_mpnn", ov="overall_rmsd", pae="mean_pae",
                 clash="af3_clash_resindices", clash_is_list=True),
    "rfd3": dict(epi="epitope_chunk_rmsd", ov="overall_rmsd", pae="mean_pae",
                 clash="af3_n_clash_res", clash_is_list=False),
    "perisland": dict(epi="epitope_chunk_rmsd", ov="overall_rmsd", pae="mean_pae",
                      clash="af3_n_clash_res", clash_is_list=False),
}


def clashlen(x):
    if x is None:
        return np.nan
    if isinstance(x, float) and np.isnan(x):
        return np.nan
    if isinstance(x, (list, np.ndarray)):
        return len(x)
    s = str(x).strip().strip("[]").split()
    return len(s) if s else 0


def load(path: Path, schema: dict) -> pd.DataFrame:
    df = pd.read_parquet(path) if path.suffix == ".parquet" else pd.read_csv(path, low_memory=False)
    cl = df[schema["clash"]].apply(clashlen) if schema["clash_is_list"] \
        else pd.to_numeric(df[schema["clash"]], errors="coerce")
    return pd.DataFrame({
        "epi": pd.to_numeric(df[schema["epi"]], errors="coerce"),
        "ov":  pd.to_numeric(df[schema["ov"]], errors="coerce"),
        "pae": pd.to_numeric(df[schema["pae"]], errors="coerce"),
        "cl":  cl,
    })


def rates(d: pd.DataFrame) -> dict:
    """Each filter as (over ALL designs, over VALID designs). VALID = all four metrics present;
    a design counts toward ALL whether or not it has every metric (a missing metric = non-pass).
    The two diverge when a set has many designs missing a metric (e.g. RFD1's missing AF3)."""
    nt = len(d)
    v = d.epi.notna() & d.ov.notna() & d.pae.notna() & d.cl.notna()
    nv = int(v.sum())
    out = {"n_total": nt, "n_valid": nv}
    masks = {"epitope RMSD <1": d.epi <= 1, "overall RMSD <2": d.ov <= 2,
             "clash = 0": d.cl == 0, "mean PAE <5": d.pae < 5,
             "all four": v & (d.epi <= 1) & (d.ov <= 2) & (d.pae < 5) & (d.cl == 0)}
    for name, m in masks.items():
        npass = int((v & m).sum()) if name != "all four" else int(m.sum())
        out[name] = (100 * npass / nt, 100 * npass / nv)   # (over all, over valid)
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--rfd1", type=Path,
                    default=Path("../known_antigen/analysis/full_run/dp2.parquet"))
    ap.add_argument("--rfd3", type=Path,
                    default=Path("../known_antigen/analysis/data/metrics_cylinder_full.csv"))
    ap.add_argument("--perisland", type=Path, default=None,
                    help="metrics_dual_island.parquet (on the cluster run dir)")
    args = ap.parse_args()

    cols = [("RFD1 / Lawson DP3 contigs", args.rfd1, "rfd1"),
            ("RFD3 / same DP3 epitopes", args.rfd3, "rfd3"),
            ("RFD3 / single island", args.perisland, "perisland")]
    res = {}
    for name, path, key in cols:
        if path is None or not Path(path).exists():
            print(f"[skip] {name}: {path} not available")
            continue
        res[name] = rates(load(Path(path), SCHEMAS[key]))

    filters = ["epitope RMSD <1", "overall RMSD <2", "clash = 0", "mean PAE <5", "all four"]
    names = list(res)
    w = max(len(f) for f in filters) + 2
    print("\ncell = % over ALL designs (% over VALID designs)")
    print(f"{'filter':<{w}}" + "".join(f"{n:>34}" for n in names))
    print(f"{'  n total / valid':<{w}}"
          + "".join(f"{format(res[n]['n_total'],',')+' / '+format(res[n]['n_valid'],','):>34}" for n in names))
    for f in filters:
        print(f"{f:<{w}}" + "".join(f"{res[n][f][0]:>10.1f} ({res[n][f][1]:.1f})  ".rjust(34) for n in names))


if __name__ == "__main__":
    main()
