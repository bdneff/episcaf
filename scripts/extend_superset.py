#!/usr/bin/env python3
"""extend_superset.py -- make dp4_superset.csv a TRUE superset: every arm's rows + every column.

`build_superset.sbatch` builds the episcaf candidate pool (C1/C2/C3, ~334,750 rows, 20 cols). John
(2026-07-20) wants the superset to actually be a superset -- the shipped library a strict subset of it --
so this appends the other two arms and unions the columns:

  + 8VDL   : all 1,280 candidates (results/dp4_8vdl_top10_allmetrics.csv), `selected` = the shipped 20.
  + LX      : the 21,759 PASSING minibinders (pulled from dp4_library.csv, which already carries them with
              their lx_* columns), `selected` = True (all passing shipped). (The full ~484k LX pool stays
              in the raw LX file; John chose candidates-plus-passing-minibinders, not the whole LX pool.)

Columns = union of the superset's 20 and the library's minibinder columns: the 20 + design_ID, model,
designedSequenceLength + the 13 lx_*. Episcaf/8VDL rows are blank in lx_*; minibinder rows are blank in
the episcaf metric + scoring columns (never scored on our axes). Runs LOCALLY -- every input is committed
or local (the C1/C2/C3 superset is read from the committed .gz), no cluster needed.

Usage:
  python scripts/extend_superset.py            # -> data/libraries/dp4_superset.csv (+ gzip it to commit)
"""
from __future__ import annotations
import argparse
import gzip
from pathlib import Path
import pandas as pd

SUPERSET_COLS = ["component", "category", "target", "predID", "island_index",
                 "epitope_rmsd", "overall_rmsd", "epitope_pae", "scaffold_pae", "mean_pae", "ptm",
                 "af3_clashes", "cylinder_clashes",
                 "composite", "rank_in_group", "is_global_pass", "selected", "library_member",
                 "sequence", "designedSequence"]
LX_COLS = ["lx_batch_uuid", "lx_uuid", "lx_binder_file", "lx_sequence", "lx_passes", "lx_plddt",
           "lx_ipae", "lx_pae", "lx_rmsd", "lx_plddt_binder", "lx_iptm", "lx_target", "lx_hotspots"]
OUT_COLS = SUPERSET_COLS + ["design_ID", "model", "designedSequenceLength"] + LX_COLS


def read_any(path: Path) -> pd.DataFrame:
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt") as fh:
            return pd.read_csv(fh, low_memory=False)
    return pd.read_csv(path, low_memory=False)


def vdl_block(metrics_csv: Path, shipped_csv: Path) -> pd.DataFrame:
    """8VDL candidates -> superset rows. rank within run by composite; selected = the shipped design_IDs."""
    a = pd.read_csv(metrics_csv, low_memory=False)
    shipped = set(pd.read_csv(shipped_csv)["design_ID"].astype(str))
    a["rank_in_group"] = (a.groupby("run")["composite"].rank(ascending=False, method="first").astype(int))
    ds = a["designedSequence"].astype(str)
    out = pd.DataFrame({
        "component": "8VDL",
        "category": "scaffolded8VDL",
        "target": "8VDL_" + a["run"].astype(str),
        "predID": a["design_ID"].astype(str),
        "design_ID": a["design_ID"].astype(str),
        "island_index": pd.NA,
        "epitope_rmsd": a["epitope_chunk_rmsd"],
        "overall_rmsd": a["overall_rmsd"],
        "epitope_pae": a["epitope_pae"],
        "scaffold_pae": pd.NA,
        "mean_pae": a["mean_pae"],
        "ptm": pd.NA,
        "af3_clashes": a["af3_n_clash_res"],       # real H/L clash, no cylinder for 8VDL
        "cylinder_clashes": pd.NA,
        "composite": a["composite"],
        "rank_in_group": a["rank_in_group"],
        "is_global_pass": a["pass_indicator"] > 0.5,
        "selected": a["design_ID"].astype(str).isin(shipped),
        "library_member": pd.NA,
        "sequence": ds.str.upper().str[:103],
        "designedSequence": ds,
        "model": "RFD",
        "designedSequenceLength": ds.str.len(),
    })
    return out


def mini_block(library_csv: Path) -> pd.DataFrame:
    """Passing minibinders, pulled straight from the library (already mapped + carrying lx_*)."""
    lib = read_any(library_csv)
    mb = lib[lib["category"] == "minibinder"].copy()
    mb["component"] = "LX"
    mb["predID"] = mb["design_ID"].astype(str)
    mb["selected"] = True                          # every passing minibinder shipped
    return mb


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--superset", default="data/libraries/dp4_superset.csv.gz",
                    help="the C1/C2/C3 episcaf superset from build_superset.sbatch (.csv or .csv.gz)")
    ap.add_argument("--vdl-metrics", default="results/dp4_8vdl_top10_allmetrics.csv")
    ap.add_argument("--vdl-shipped", default="results/dp4_8vdl_top10.csv")
    ap.add_argument("--library", default="data/libraries/dp4_library.csv")
    ap.add_argument("--out", default="data/libraries/dp4_superset.csv")
    args = ap.parse_args()

    base = read_any(Path(args.superset))           # C1/C2/C3 (or an already-extended file)
    # idempotent: keep only the episcaf pool from base, so re-running on the full superset re-appends
    # the 8VDL/LX blocks fresh rather than duplicating them.
    base = base[~base["component"].isin(["8VDL", "LX"])].copy()
    vdl = vdl_block(Path(args.vdl_metrics), Path(args.vdl_shipped))
    mini = mini_block(Path(args.library))

    full = pd.concat([base, vdl, mini], ignore_index=True)
    for c in OUT_COLS:                              # ensure every union column exists, fill gaps
        if c not in full.columns:
            full[c] = pd.NA
    full = full[OUT_COLS]

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    full.to_csv(args.out, index=False)
    n = full["component"].value_counts()
    print(f"[extend] {len(full):,} rows, {len(OUT_COLS)} cols -> {args.out}")
    print(f"[extend] by arm: C1/C2/C3={int(n.reindex(['C1','C2','C3']).sum())} "
          f"8VDL={int(n.get('8VDL',0))} LX={int(n.get('LX',0))}")
    print(f"[extend] selected total: {int(full['selected'].fillna(False).sum())}")


if __name__ == "__main__":
    main()
