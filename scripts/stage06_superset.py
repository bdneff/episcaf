#!/usr/bin/env python3
"""Build the ALL-DESIGNS superset of the DP4 scaffolded arms (John, 2026-07-16).

The shipped dp4_library.csv holds only the SELECTED designs. This emits EVERY candidate design
in the scaffolded pools, in dp4_library.csv's own column shape, so the overall distributions can
be poked at and the selected rows compared against the pool they came from:

  component            C1 / C2 / C3
  category, target, predID, island_index          identity (as in dp4_library.csv)
  epitope_rmsd, overall_rmsd,
  epitope_pae, scaffold_pae, mean_pae             the PAE decomposition (see PAE below)
  af3_clashes, cylinder_clashes, ptm              the metrics (blank where a component lacks one)
  composite, rank_in_group                        soft-gate composite + rank within selection group
  is_global_pass                                  clears ALL four Lawson filters (C1/C2 only)
  selected, library_member                        did this design ship, and as which member
  sequence, designedSequence                      see SEQUENCES below

Ranked under `antibody_softgate` (C1/C2) -- the preset that ACTUALLY picked the shipped library --
so `composite`/`rank_in_group` reconcile with what shipped. C3 uses `twelvemer` (no antibody).

PAE. `mean_pae` is the OVERALL PAE -- the whole AF3 PAE matrix averaged -- and it is the metric the
published RFD/Lawson four-filter thresholds at < 5 (kept here as the point of comparison, not gated).
`epitope_pae` is the epitope block only (what the soft-gate weights, 0.10); `scaffold_pae` is the
rest, so overall reads as its two parts. All are PAE (inter-residue error), NOT pLDDT (per-residue
confidence, the AF3-viewer coloring) -- we do not currently extract pLDDT.

SEQUENCES. Filled for SELECTED designs (copied verbatim from dp4_library.csv, so the superset
agrees with what shipped by construction) and, with --sequences, for GLOBAL-PASSING designs (read
from each design's AF3 chain A). Blank otherwise: filling all ~335k would mean reading every design
PDB, and distributions live in the metrics, not the sequences. `designedSequence` (the case-encoded
form) is only available for selected designs -- case-encoding was only ever run on those.

JOIN NOTE. dp4_library.csv's `design_ID` is `C<n>_<predID>_r<N>`: a component prefix plus an AF3
replicate suffix that `predID` does not carry. Stripping both is 1:1 onto predID (verified: 1120/1120
C1, 1660/1660 C2, 4390/4390 C3), which is how `selected` is joined.

C1 and C3 pools are local; C2 (metrics_dual_island.parquet) lives on the cluster -- run this there
for C2 and concat with --append. Usage (one component at a time):

  python scripts/stage06_superset.py --component C1 \
    --metrics-csv $D/known_antigen/analysis/data/metrics_whole_epitope_103.csv \
    --library data/libraries/dp4_library.csv \
    --out data/libraries/dp4_superset.csv
  python scripts/stage06_superset.py --component C3 \
    --metrics-csv $D/12mer_tiling/analysis/data/metrics_12mer.csv --append \
    --library data/libraries/dp4_library.csv \
    --out data/libraries/dp4_superset.csv

NOT INCLUDED: the 8VDL arm (20 shipped). Its run is separate and small; add it as a COMPONENTS
entry once its metrics file is pointed at.
"""
from __future__ import annotations

import argparse
import os
import re
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from episcaf_analysis.presets import PRESETS
from episcaf_analysis.score import score

# component -> how it is scored, grouped, and identified
COMPONENTS = {
    "C1": dict(category="scaffoldedAbEpitope",    preset="antibody_softgate",
               group=["id"],               target="id",      has_ab=True),
    "C2": dict(category="scaffoldedSingleIsland", preset="antibody_softgate",
               group=["id", "island_index"], target="id",    has_ab=True),
    "C3": dict(category="scaffoldedPolyclonal",   preset="twelvemer",
               group=["antigen", "id"],    target="antigen", has_ab=False, seq_col="design_seq"),
}

# metrics name -> dp4_library.csv name
RENAME = {
    "epitope_chunk_rmsd": "epitope_rmsd",
    "af3_n_clash_res": "af3_clashes",
    "cylinder_native_aware": "cylinder_clashes",
}
# PAE is carried as the full three-way decomposition AF3 gives us, none of it filtered on here:
#   mean_pae     -- OVERALL, the whole PAE matrix averaged. This is the metric the published RFD
#                   filters use (four-filter: < 5); kept as the point of comparison, not a gate.
#   epitope_pae  -- the epitope block only (what the soft-gate scorer actually weights, at 0.10)
#   scaffold_pae -- the rest of the chain, so overall can be read as its two parts
OUT_COLS = ["component", "category", "target", "predID", "island_index",
            "epitope_rmsd", "overall_rmsd", "epitope_pae", "scaffold_pae", "mean_pae", "ptm",
            "af3_clashes", "cylinder_clashes",
            "composite", "rank_in_group", "is_global_pass",
            "selected", "library_member", "sequence", "designedSequence"]
# the four-filter (global mean PAE, exactly as sec:fourfilter)
FOUR_FILTER = {"epitope_chunk_rmsd": 1.0, "overall_rmsd": 2.0, "mean_pae": 5.0, "af3_n_clash_res": 0.0}


# collapse the doubled RFD3-name block in a predID (same rule as stage06_assemble.py -- keep in sync so
# the superset's predID matches the library's collapsed design_ID for the `selected` join).
_DOUBLE = re.compile(r"(.+?)_\1(_\d+_model)")
def collapse_id(s: str) -> str:
    return _DOUBLE.sub(r"\1\2", str(s))


def lib_key(design_id: str) -> str:
    """`C1_<predID>_r3` -> `<predID>` (design_ID already collapsed upstream). See JOIN NOTE."""
    return re.sub(r"_r\d+$", "", re.sub(r"^C\d+_", "", str(design_id)))


def add_sequences(scored: pd.DataFrame, mask: pd.Series, remap: tuple[str, str] | None = None) -> int:
    """Read AF3 chain-A sequence for `mask` rows. Cluster-only (needs gemmi + the AF3 outputs).

    `af3_dir` holds ABSOLUTE paths baked in when the metrics were built -- for C1/C2 that was
    /scratch, which is swept on ~30 days. The runs now live under $WS, so `remap=(old, new)`
    rewrites that prefix when the recorded directory is gone. See --af3-remap.
    """
    from pathlib import Path
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                    "episcaf_analysis"))
    import compute_metrics as CM  # noqa: E402

    n = 0
    for i in scored.index[mask]:
        d = scored.at[i, "af3_dir"]
        if not isinstance(d, str) or not d:
            continue
        if remap and not os.path.isdir(d) and d.startswith(remap[0]):
            d = remap[1] + d[len(remap[0]):]
        cif, _, _ = CM.find_af3_files(Path(d))
        if cif is None:
            continue
        scored.at[i, "sequence"] = CM.chain_seq(CM.get_chain(CM.read_structure(cif), "A"))
        n += 1
    return n


def build(component: str, path: str, library: str | None, want_seqs: bool,
          remap: tuple[str, str] | None = None) -> pd.DataFrame:
    spec = COMPONENTS[component]
    df = pd.read_parquet(path) if path.endswith(".parquet") else pd.read_csv(path, low_memory=False)
    if "status" in df.columns:
        df = df[df["status"].astype(str).str.lower().eq("ok")].copy()
    print(f"[{component}] {len(df)} designs from {path}")

    # composite under the preset that actually picked the library (rank ALL rows; no top-k cut)
    preset = {k: (v.copy() if isinstance(v, dict) else v) for k, v in PRESETS[spec["preset"]].items()}
    preset["select"] = None
    scored = score(df.copy(), preset)

    # rank within the shipped selection group
    gcols = [g for g in spec["group"] if g in scored.columns] or [spec["group"][0]]
    if len(gcols) > 1:
        scored["_grp"] = scored[gcols].astype(str).agg("|".join, axis=1)
        gkey = "_grp"
    else:
        gkey = gcols[0]
    scored["rank_in_group"] = scored.groupby(gkey)["composite"].rank(ascending=False, method="first")

    # global four-filter pass (known-antibody arms only)
    if spec["has_ab"] and all(c in scored.columns for c in FOUR_FILTER):
        ok = pd.Series(True, index=scored.index)
        for c, t in FOUR_FILTER.items():
            x = pd.to_numeric(scored[c], errors="coerce")
            ok &= (x == 0) if t == 0.0 else (x <= t)
        scored["is_global_pass"] = ok
    else:
        scored["is_global_pass"] = pd.NA

    # design identity. C1/C2 carry `predID`; C3 does not, but every component's af3_dir basename IS
    # the design id (and is what dp4_library.csv's design_ID wraps -- see JOIN NOTE).
    if "predID" not in scored.columns:
        if "af3_dir" not in scored.columns:
            raise SystemExit(f"[{component}] metrics has neither predID nor af3_dir; cannot identify designs")
        scored["predID"] = scored["af3_dir"].astype(str).map(os.path.basename)
    scored["predID"] = scored["predID"].astype(str).map(collapse_id)   # match the library's collapsed IDs

    scored["component"] = component
    scored["category"] = spec["category"]
    scored["target"] = scored[spec["target"]] if spec["target"] in scored.columns else pd.NA
    scored = scored.rename(columns=RENAME)

    # selected + the shipped sequences, joined from dp4_library.csv (see JOIN NOTE)
    scored["selected"] = False
    for c in ("library_member", "sequence", "designedSequence"):
        scored[c] = pd.NA
    if library:
        lib = pd.read_csv(library, low_memory=False)
        lib = lib[lib["category"].eq(spec["category"])].copy()
        lib["_k"] = lib["design_ID"].map(lib_key)
        keyed = lib.set_index("_k")
        hit = scored["predID"].astype(str).map(lambda p: p in keyed.index)
        scored["selected"] = hit
        for c in ("library_member", "sequence", "designedSequence"):
            scored.loc[hit, c] = scored.loc[hit, "predID"].astype(str).map(keyed[c])
        print(f"[{component}] selected: {int(hit.sum())} / {len(lib)} shipped members matched")
        if int(hit.sum()) != len(lib):
            print(f"[{component}] WARNING: {len(lib) - int(hit.sum())} shipped members did not "
                  f"match a design in this metrics file -- check the join", file=sys.stderr)

    # some components already carry the design sequence in their metrics (C3's design_seq) -- free,
    # so take it for every row rather than only the passing ones.
    sc = spec.get("seq_col")
    if sc and sc in scored.columns:
        miss = scored["sequence"].isna()
        scored.loc[miss, "sequence"] = scored.loc[miss, sc]
        print(f"[{component}] sequence from metrics `{sc}` for {int(miss.sum())} rows")

    # sequences for global-passing designs that did not ship (cluster only)
    if want_seqs:
        need = scored["is_global_pass"].fillna(False).astype(bool) & scored["sequence"].isna()
        if "af3_dir" not in scored.columns:
            print(f"[{component}] WARNING: no af3_dir column; cannot fill sequences", file=sys.stderr)
        elif need.any():
            got = add_sequences(scored, need, remap)
            print(f"[{component}] read AF3 sequences for {got} / {int(need.sum())} passing designs")
            # Resolving nothing means the AF3 outputs moved or were swept -- not "no passers".
            # Fail loudly rather than writing a superset whose `sequence` column is quietly empty.
            if got == 0:
                raise SystemExit(
                    f"[{component}] asked for {int(need.sum())} AF3 sequences and resolved 0.\n"
                    f"  af3_dir example: {scored.loc[need, 'af3_dir'].iloc[0]}\n"
                    f"  If the run moved (e.g. /scratch -> $WS), pass --af3-remap OLD:NEW.")

    for c in OUT_COLS:
        if c not in scored.columns:
            scored[c] = pd.NA
    return scored[OUT_COLS]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--component", required=True, choices=sorted(COMPONENTS))
    ap.add_argument("--metrics-csv", required=True)
    ap.add_argument("--library", help="dp4_library.csv -- adds `selected` + the shipped sequences")
    ap.add_argument("--sequences", action="store_true",
                    help="also read AF3 chain-A sequence for global-passing designs (cluster only)")
    ap.add_argument("--af3-remap", metavar="OLD:NEW", default="",
                    help="rewrite this af3_dir path prefix when the recorded directory is gone "
                         "(the metrics bake in absolute paths; C1/C2 moved /scratch -> $WS). "
                         "Same OLD:NEW form as case_encode_c2.py --af3-remap / "
                         "case_encode_selected.py --pdb-remap")
    ap.add_argument("--out", required=True)
    ap.add_argument("--append", action="store_true", help="append to --out instead of overwriting")
    args = ap.parse_args()

    remap = None
    if args.af3_remap:
        if ":" not in args.af3_remap:
            ap.error("--af3-remap wants OLD:NEW")
        old, new = args.af3_remap.split(":", 1)
        remap = (old, new)

    out = build(args.component, args.metrics_csv, args.library, args.sequences, remap)
    header = not (args.append and os.path.exists(args.out))
    out.to_csv(args.out, mode="a" if args.append else "w", header=header, index=False)
    gp = out["is_global_pass"]
    print(f"[{args.component}] wrote {len(out)} rows -> {args.out}"
          f"  (global-pass: {int(gp.sum()) if gp.notna().any() else 'n/a'},"
          f" selected: {int(out['selected'].sum())},"
          f" with sequence: {int(out['sequence'].notna().sum())})")


if __name__ == "__main__":
    main()
