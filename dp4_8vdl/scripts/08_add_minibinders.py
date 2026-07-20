#!/usr/bin/env python3
"""08_add_minibinders.py -- fold the LX PfEMP1/EPCR minibinders into dp4_library.csv.

John's ask (2026-07-20): keep ONE dp4_library.csv with everything in it, so the overall numbers add up
at a glance. The LX designs are the de-novo minibinder arm of the same PfEMP1/EPCR project as the 8VDL
scaffolds -- hence this step lives in the 8vdl subdir.

What we add: only the LX designs that PASS their own filters (`passes == True`) -- ~21.8k of ~484k
generations -- mapped into the library's column shape. We carry NO episcaf metrics for them: they were
never scored on our epitope-RMSD / PAE / clash / cylinder axes (John: "we dont need metrics for those"),
and the LX file's own metrics (plddt/ipae/pae/rmsd/iptm) are a different measurement, so those five
library columns stay BLANK rather than being cross-populated from unlike quantities.

Column map (LX -> library):
  sequence           -> sequence           (103-mer, already uppercase = synthesizable form)
  sequence           -> designedSequence   (no epitope case-encoding: a minibinder presents no grafted
                                             island, so designedSequence == the plain sequence)
  len(sequence)      -> designedSequenceLength
  uuid               -> design_ID
  target             -> target              (kept verbatim, e.g. fold_pfemp1_epcr_model_0, for traceability)
  (constant "LX")    -> model
  (constant CATEGORY)-> category
  epitope_rmsd, overall_rmsd, epitope_pae, af3_clashes, cylinder_clashes -> BLANK (our axes, never scored)
  EVERY LX column     -> lx_<name>  (all native LatentX columns kept for post-hoc analysis: plddt, pae,
                                     rmsd, ipae, iptm, plddt_binder, hotspots, uuid, ... ; episcaf rows NaN)

Idempotent: strips any existing rows of this category first, appends fresh, then renumbers
library_member so the episcaf block keeps DP4_1..DP4_<n> (unchanged order) and the minibinders continue
DP4_<n+1>.. . Run AFTER scripts/stage06_assemble.py; re-running is safe.

Usage:
  python dp4_8vdl/scripts/08_add_minibinders.py --lx /path/to/LX_YYYYMMDD.csv
  # writes data/libraries/dp4_library.csv in place (use --out to write elsewhere)
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

CATEGORY = "minibinder"
# The standard library schema (kept in sync with scripts/stage06_assemble.py). 8 identity columns, then
# the metric columns, then the scoring/identity extras. Minibinders fill only the 8 identity columns;
# everything after `target` is blank for them (never scored on our axes).
LIB_COLS = ["library_member", "sequence", "category", "model", "designedSequence",
            "designedSequenceLength", "design_ID", "target",
            "epitope_rmsd", "overall_rmsd", "epitope_pae", "scaffold_pae", "mean_pae", "ptm",
            "af3_clashes", "cylinder_clashes",
            "composite", "rank_in_group", "is_global_pass", "island_index"]
METRIC_COLS = LIB_COLS[8:]      # everything after the 8 identity columns -> blank for minibinders


def read_lx(path: Path) -> pd.DataFrame:
    """Load the LX generations and keep only the filter-passing designs."""
    lx = pd.read_csv(path, low_memory=False)
    for need in ("passes", "sequence", "uuid", "target"):
        if need not in lx.columns:
            raise SystemExit(f"[minibinders] LX file missing column {need!r}; have {list(lx.columns)}")
    p = lx["passes"]
    keep = p if p.dtype == bool else p.astype(str).str.strip().str.lower().isin(("true", "1"))
    passing = lx[keep].copy()
    if passing.empty:
        raise SystemExit(f"[minibinders] 0 rows with passes==True in {path} -- nothing to add")
    return passing


def to_library_rows(passing: pd.DataFrame, category: str) -> pd.DataFrame:
    passing = passing.reset_index(drop=True)
    seq = passing["sequence"].astype(str).str.upper()
    out = pd.DataFrame({
        "library_member": pd.NA,                    # assigned after concat
        "sequence": seq,
        "category": category,
        "model": "LX",
        "designedSequence": seq,                    # no epitope casing for a de-novo binder
        "designedSequenceLength": seq.str.len(),
        "design_ID": passing["uuid"].astype(str),
        "target": passing["target"].astype(str),
    })
    for c in METRIC_COLS:
        out[c] = pd.NA                              # never scored on OUR axes -> blank, not imputed
    out = out[LIB_COLS]
    # Carry EVERY native LatentX column (John, 2026-07-20: keep them for post-hoc analysis), prefixed
    # `lx_` so they never collide with the library schema and their provenance is explicit. Episcaf rows
    # get NaN here (they have no LX metrics); minibinder rows carry the full LX record (plddt, pae, rmsd,
    # ipae, iptm, plddt_binder, hotspots, uuid, batch_uuid, binder_file, ...).
    lx = passing.add_prefix("lx_")
    return pd.concat([out, lx], axis=1)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lx", required=True, help="LX_YYYYMMDD.csv (the ~484k generations; passing subset is added)")
    ap.add_argument("--library", default="data/libraries/dp4_library.csv")
    ap.add_argument("--out", default=None, help="default: overwrite --library in place")
    ap.add_argument("--category", default=CATEGORY)
    args = ap.parse_args()

    lib = pd.read_csv(args.library, low_memory=False)
    missing = [c for c in LIB_COLS if c not in lib.columns]
    if missing:
        raise SystemExit(f"[minibinders] library missing standard columns {missing}\n  found {list(lib.columns)}")

    # idempotent: drop any prior minibinder block, and keep only the STANDARD columns for the episcaf
    # rows (strips any lx_* left from a previous run) so the episcaf side never carries stale LX data.
    episcaf = lib[lib["category"] != args.category][LIB_COLS].copy()
    mini = to_library_rows(read_lx(Path(args.lx)), args.category)

    combined = pd.concat([episcaf, mini], ignore_index=True)          # episcaf lx_* -> NaN
    combined["library_member"] = [f"DP4_{i}" for i in range(1, len(combined) + 1)]
    lxcols = [c for c in combined.columns if c.startswith("lx_")]     # standard 13 first, then lx_*
    combined = combined[LIB_COLS + lxcols]

    out = Path(args.out or args.library)
    combined.to_csv(out, index=False)
    print(f"[minibinders] episcaf {len(episcaf)} + {args.category} {len(mini)} = {len(combined)} rows, "
          f"{len(lxcols)} lx_ columns carried -> {out}")
    print("[minibinders] category counts:")
    print(combined["category"].value_counts().to_string())


if __name__ == "__main__":
    main()
