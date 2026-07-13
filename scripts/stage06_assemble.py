#!/usr/bin/env python3
"""
stage06_assemble.py -> the DP4 synthesis library (06_library), one 8-column DP2 annotated file.

Concatenates the six components into the approved DP2 schema
(library_member, sequence, category, model, designedSequence, designedSequenceLength, design_ID, target)
with global `library_member` numbering. Along the way it applies, per the decisions in
docs/DP4_LIBRARY.md:

  * 56-mAb EXCLUSION (drop 2h32/4xwo/7a3t) -- only C1/C2 need it here (C3 polyclonal is unaffected;
    C4/C5/C6 were already rebuilt on the 56-set).
  * DEPTH cut (top-n per group) on the ranked components C1/C2/C3 (--depth, default 20). C4 exhaustive,
    C5 fixed sample, C6 derived from C1's depth (rebuild C6 if --depth != its build depth).
  * 104->103 epitope-PRESERVING trim on the 104-mer designs (C1/C5/C6): drop one residue from whichever
    terminus is scaffold (lowercase); if a design has epitope (uppercase) at BOTH termini it can't be
    trimmed -- it is dropped and logged (e.g. 3ux9_1P rank 9; its rank-(n+1) replacement is only available
    if the case-encode went deeper than n). C2/C3/C4 are natively 103.

`sequence` = the synthesized 103-mer (uppercase of the trimmed case-encoded string);
`designedSequence` = the case-encoded (epitope UPPER / scaffold lower) trimmed string, for visualization.
C4 is already in DP2 format (its `sequence`/`designedSequence` are used verbatim).

Usage:
  python scripts/stage06_assemble.py --depth 20 --out data/libraries/dp4_library.csv
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

R = Path(__file__).resolve().parents[1]
EXCLUDE = ("2h32", "4xwo", "7a3t")   # canonical 56-mAb exclusion (id-prefix)


def trim_103(se: str):
    """Epitope-preserving 104->103 trim. Returns (trimmed_str, note) or (None, reason) if untrimmable."""
    if len(se) <= 103:
        return se, "already103"
    if se[-1].islower():
        return se[:-1], "cterm"          # drop scaffold C-terminus (default)
    if se[0].islower():
        return se[1:], "nterm"           # drop scaffold N-terminus (C-term is epitope)
    return None, "both_termini_epitope"  # can't trim -> drop


def attach_rank(seq_df, ranked_csv):
    """The scaffoldEPITOPE file is row-aligned to its ranked table (case-encode iterated it in order);
    attach rank_in_group by position (asserted)."""
    rk = pd.read_csv(ranked_csv, low_memory=False)
    assert len(rk) == len(seq_df), f"row count mismatch {ranked_csv}: {len(rk)} vs {len(seq_df)}"
    out = seq_df.copy()
    out["rank_in_group"] = rk["rank_in_group"].to_numpy()
    return out


def scaffold_rows(comp, seq_csv, ranked_csv, category, *, trim, exclude, depth):
    """Build DP2 rows for a scaffold component from its scaffoldEPITOPE file."""
    d = pd.read_csv(seq_csv, low_memory=False)
    d = d[d.get("status", "ok").eq("ok")] if "status" in d else d
    if ranked_csv:
        d = attach_rank(d, ranked_csv)
        d = d[d.rank_in_group <= depth]
    if exclude:
        tgt = d["target"].astype(str).str.lower()
        d = d[~tgt.str.startswith(EXCLUDE)]
    rows, dropped = [], []
    for i, r in enumerate(d.itertuples(index=False)):
        se = str(getattr(r, "scaffoldEPITOPE"))
        tse, note = trim_103(se) if trim else (se, "native103")
        if tse is None:
            dropped.append((getattr(r, "target", "?"), note)); continue
        did = f"{comp}_{getattr(r,'token',getattr(r,'predID',i))}"
        if hasattr(r, "rank_in_group"):
            did += f"_r{int(r.rank_in_group)}"
        rows.append(dict(sequence=tse.upper(), category=category, model="RFD",
                         designedSequence=tse, designedSequenceLength=len(tse),
                         design_ID=did, target=getattr(r, "target", "")))
    return pd.DataFrame(rows), dropped


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--depth", type=int, default=20, help="top-n per group for the ranked components C1/C2/C3")
    ap.add_argument("--out", default="data/libraries/dp4_library.csv")
    args = ap.parse_args()
    res = R / "results"; lib = R / "data/libraries"
    parts, all_dropped = [], {}

    # C1/C2/C3 -- ranked scaffolds (depth-cut); C1 also 104->trim + exclusion; C2 exclusion; C3 neither
    for comp, seqf, rankf, cat, trim, exc in [
        ("C1", res/"dp4_C1_scaffoldEPITOPE.csv", res/"dp4_C1_whole_epitope_ranked.top20.csv", "scaffoldedAbEpitope", True, True),
        ("C2", res/"dp4_C2_scaffoldEPITOPE.csv", res/"dp4_C2_single_island_ranked.top20.csv", "scaffoldedSingleIsland", False, True),
        ("C3", res/"dp4_C3_scaffoldEPITOPE.csv", res/"dp4_C3_12mer_ranked.top20.csv", "scaffoldedPolyclonal", False, False),
    ]:
        df, dropped = scaffold_rows(comp, seqf, rankf, cat, trim=trim, exclude=exc, depth=args.depth)
        parts.append(df); all_dropped[comp] = dropped

    # C5 -- fixed sample (no depth cut), 104->trim, already 56
    df, dropped = scaffold_rows("C5", res/"dp4_C5_scaffoldEPITOPE.csv", None, "metricSpaceTitration",
                                trim=True, exclude=False, depth=args.depth)
    parts.append(df); all_dropped["C5"] = dropped

    # C6 -- controls, 104->trim, already 56 (built at top-20; rebuild if depth!=20). Uses its own DP2 cols.
    c6 = pd.read_csv(res/"dp4_C6_controls.csv", low_memory=False)
    c6rows, c6drop = [], []
    for r in c6.itertuples(index=False):
        se = str(getattr(r, "scaffoldEPITOPE")); tse, note = trim_103(se)
        if tse is None: c6drop.append((getattr(r, "design_ID", "?"), note)); continue
        c6rows.append(dict(sequence=tse.upper(), category="scaffoldedEpitopeControl", model="RFD",
                           designedSequence=tse, designedSequenceLength=len(tse),
                           design_ID=getattr(r, "design_ID"), target=getattr(r, "target", "")))
    parts.append(pd.DataFrame(c6rows)); all_dropped["C6"] = c6drop

    # C4 -- already 8-column format, native 103, already 56
    keep8 = ["sequence", "category", "model", "designedSequence",
             "designedSequenceLength", "design_ID", "target"]
    c4 = pd.read_csv(lib/"dp4_tiled30mers_fasta.csv", low_memory=False)
    parts.append(c4[keep8])

    # 8VDL arm -- already 8-column (07_consolidate), native 103, fixed top-10 per definition
    v8 = pd.read_csv(res/"dp4_8vdl_top10.csv", low_memory=False)
    parts.append(v8[keep8])

    lib_df = pd.concat(parts, ignore_index=True)
    lib_df.insert(0, "library_member", [f"DP4_{i}" for i in range(1, len(lib_df) + 1)])
    assert (lib_df.sequence.str.len() <= 103).all(), "a sequence exceeds 103 residues!"
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    lib_df.to_csv(args.out, index=False)

    print(f"[assemble] depth={args.depth} -> {len(lib_df):,} library members -> {args.out}")
    print("[assemble] per-component:")
    for p, name in zip(parts, ["C1", "C2", "C3", "C5", "C6", "C4", "8VDL"]):
        print(f"    {name}: {len(p):,}")
    print(f"[assemble] sequence lengths: {lib_df.sequence.str.len().value_counts().to_dict()}")
    for comp, dr in all_dropped.items():
        if dr:
            print(f"[assemble] {comp}: {len(dr)} dropped (untrimmable, both-termini epitope): "
                  f"{dr[:3]}{'...' if len(dr)>3 else ''}")


if __name__ == "__main__":
    main()
