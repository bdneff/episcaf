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
import re
from pathlib import Path
import pandas as pd

R = Path(__file__).resolve().parents[1]
EXCLUDE = ("2h32", "4xwo", "7a3t")   # canonical 56-mAb exclusion (id-prefix)

# The predID doubles the RFD3 backbone name (an AF3-output-dir naming artifact): <name>_<name>_<d>_model
# ... for every design-producing arm (C1/C2/C5 and 8VDL, which double differently). Collapse the block
# that repeats immediately before the `_<digit>_model` structure tag. Cosmetic; sequences/metrics
# untouched. Kept identical in stage06_superset.py / extend_superset.py so the superset join still matches.
_DOUBLE = re.compile(r"(.+?)_\1(_\d+_model)")
def collapse_id(s: str) -> str:
    return _DOUBLE.sub(r"\1\2", str(s))

# The metric fields shipped alongside the 8 standard columns. Maps the internal metric name -> the
# shipped column name. Left EMPTY where a design has no such value: C3 has no af3_clashes (no antibody);
# C4 (linear tiles) and C6 (mutants) were never folded, so all are blank; 8VDL has no cylinder; only
# C1/C2/C3 carry the PAE decomposition + ptm. pandas writes NaN as an empty CSV cell. (John, 2026-07-20:
# don't condense -- carry the full metric set the superset has, not just the lean 5.)
METRICS = {
    "epitope_chunk_rmsd":    "epitope_rmsd",
    "overall_rmsd":          "overall_rmsd",
    "epitope_pae":           "epitope_pae",
    "scaffold_pae":          "scaffold_pae",
    "mean_pae":              "mean_pae",
    "ptm":                   "ptm",
    "af3_n_clash_res":       "af3_clashes",
    "cylinder_native_aware": "cylinder_clashes",
}
# Scoring / identity columns carried through from the ranked tables (C1/C2/C3), blank elsewhere.
EXTRAS = ["composite", "rank_in_group", "is_global_pass", "island_index"]


def trim_103(se: str):
    """Epitope-preserving 104->103 trim. Returns (trimmed_str, note) or (None, reason) if untrimmable."""
    if len(se) <= 103:
        return se, "already103"
    if se[-1].islower():
        return se[:-1], "cterm"          # drop scaffold C-terminus (default)
    if se[0].islower():
        return se[1:], "nterm"           # drop scaffold N-terminus (C-term is epitope)
    return None, "both_termini_epitope"  # can't trim -> drop


def attach_source(seq_df, source_csv):
    """The scaffoldEPITOPE file is row-aligned to its ranked/metrics table (case-encode iterated it in
    order); attach the metric columns + the scoring/identity extras by position (asserted)."""
    src = pd.read_csv(source_csv, low_memory=False)
    assert len(src) == len(seq_df), f"row count mismatch {source_csv}: {len(src)} vs {len(seq_df)}"
    out = seq_df.copy()
    for internal, shipped in METRICS.items():
        if internal in src:
            out[shipped] = src[internal].to_numpy()
    for c in ("composite", "rank_in_group", "island_index", "af3_clash_status"):
        if c in src:
            out[c] = src[c].to_numpy()
    if "pass_indicator" in src:      # -> is_global_pass, the four-filter soft-AND crossing 0.5
        out["is_global_pass"] = (pd.to_numeric(src["pass_indicator"], errors="coerce") > 0.5).to_numpy()
    return out


def scaffold_rows(comp, seq_csv, source_csv, category, *, trim, exclude, depth):
    """Build DP2 rows for a scaffold component from its scaffoldEPITOPE file, carrying the 5 scoring
    metrics from `source_csv` (its ranked table, or for C5 its titration table). `depth` is applied
    only when the source has rank_in_group (ranked components); pass None to keep all rows (C5)."""
    d = pd.read_csv(seq_csv, low_memory=False)
    d = d[d.get("status", "ok").eq("ok")] if "status" in d else d
    if source_csv:
        d = attach_source(d, source_csv)
        if "rank_in_group" in d and depth is not None:
            d = d[d.rank_in_group <= depth]
        # Drop designs whose accessibility could not be computed: a <3-residue island gives too few
        # epitope CA pairs to define the design->native rigid-body fit, so both the real clash and the
        # cylinder come out blank (stage05: too_few_epitope_pairs). John, 2026-07-21: cull these -- a
        # 2-residue island scaffolded alone is a marginal target and ships with no accessibility metric.
        if "af3_clash_status" in d:
            n0 = len(d)
            d = d[d["af3_clash_status"].astype(str) != "too_few_epitope_pairs"]
            if n0 - len(d):
                print(f"[{comp}] dropped {n0-len(d)} designs with uncomputable accessibility "
                      f"(too_few_epitope_pairs)")
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
        row = dict(sequence=tse.upper(), category=category, model="RFD",
                   designedSequence=tse, designedSequenceLength=len(tse),
                   design_ID=did, target=getattr(r, "target", ""))
        for shipped in METRICS.values():
            row[shipped] = getattr(r, shipped, float("nan"))
        for c in EXTRAS:
            row[c] = getattr(r, c, float("nan"))
        rows.append(row)
    return pd.DataFrame(rows), dropped


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--depth", type=int, default=20, help="top-n per group (mAb/island) for C1/C2 (and C6 via C1)")
    ap.add_argument("--c3-depth", type=int, default=10,
                    help="top-n per 12-mer window for C3 (default 10). C3 tiles are 12-mers stepping by "
                         "2 residues, so windows overlap heavily; kept shallow, but set to 10 to maximise "
                         "polyclonal coverage given the weaker clash distribution (John, 2026-07-14)")
    ap.add_argument("--out", default="data/libraries/dp4_library.csv")
    args = ap.parse_args()
    res = R / "results"; lib = R / "data/libraries"
    parts, all_dropped = [], {}

    # C1/C2/C3 -- ranked scaffolds (depth-cut); C1 also 104->trim + exclusion; C2 exclusion; C3 neither.
    # C1/C2 take the shipping --depth; C3 is fixed at --c3-depth (3) because neighbouring tiles overlap.
    for comp, seqf, rankf, cat, trim, exc, dep in [
        ("C1", res/"dp4_C1_scaffoldEPITOPE.csv", res/"dp4_C1_whole_epitope_ranked.top20.csv", "scaffoldedAbEpitope", True, True, args.depth),
        ("C2", res/"dp4_C2_scaffoldEPITOPE.csv", res/"dp4_C2_single_island_ranked.top20.csv", "scaffoldedSingleIsland", False, True, args.depth),
        ("C3", res/"dp4_C3_scaffoldEPITOPE.csv", res/"dp4_C3_12mer_ranked.top20.csv", "scaffoldedPolyclonal", False, False, args.c3_depth),
    ]:
        df, dropped = scaffold_rows(comp, seqf, rankf, cat, trim=trim, exclude=exc, depth=dep)
        parts.append(df); all_dropped[comp] = dropped

    # C5 -- fixed sample (no depth cut), native 103. Metrics come from the titration table (row-aligned).
    df, dropped = scaffold_rows("C5", res/"dp4_C5_scaffoldEPITOPE.csv", res/"dp4_C5_titration.csv",
                                "metricSpaceTitration", trim=True, exclude=False, depth=None)
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

    # C4 -- linear tiles, native 103, already 56. Two touch-ups: (1) designedSequence is made the FULL
    # 103-mer construct in EPITOPEscaffold format -- the 30-mer tile (last 30 residues) uppercase (it is
    # the epitope), the GSGA.. filler + ENLYFQGA TEV lowercase (scaffold) -- so every row's
    # designedSequence is a 103-mer in the same casing (John, 2026-07-14); (2) design_ID is a per-antigen
    # tile-start index (1,7,13,...) repeating across antigens, so namespace it C4_<target>_t<pos>.
    keep8 = ["sequence", "category", "model", "designedSequence",
             "designedSequenceLength", "design_ID", "target"]
    c4 = pd.read_csv(lib/"dp4_tiled30mers_fasta.csv", low_memory=False).copy()
    c4["designedSequence"] = c4["sequence"].map(lambda s: s[:-30].lower() + s[-30:].upper())
    c4["designedSequenceLength"] = c4["designedSequence"].str.len()
    c4["design_ID"] = "C4_" + c4["target"].astype(str) + "_t" + c4["design_ID"].astype(str)
    parts.append(c4[keep8])                     # metrics blank (linear controls, never folded)

    # 8VDL arm -- 8-column + (after 07_consolidate rerun) the 5 metric columns; native 103, top-10/def
    v8 = pd.read_csv(res/"dp4_8vdl_top10.csv", low_memory=False)
    v8_cols = keep8 + [c for c in list(METRICS.values()) + EXTRAS if c in v8.columns]
    parts.append(v8[v8_cols])

    lib_df = pd.concat(parts, ignore_index=True)

    # Dedup by peptide sequence, keeping the FIRST occurrence. Parts are concatenated C1,C2,C3,C5,C6,C4,
    # 8VDL, so a peptide picked by both C1 (top ranking) and C5 (metric-space sample of the same pool) is
    # kept as its C1 row and the C5 duplicate is dropped -- no point ordering the same peptide twice
    # (John, 2026-07-21). The C5 point is still covered: the same peptide is assayed under its C1 barcode.
    n0 = len(lib_df)
    lib_df = lib_df.drop_duplicates(subset="sequence", keep="first").reset_index(drop=True)
    print(f"[assemble] dropped {n0-len(lib_df)} duplicate-sequence rows (picked by >1 component)")

    lib_df["design_ID"] = lib_df["design_ID"].map(collapse_id)   # un-double the AF3-artifact predIDs
    lib_df.insert(0, "library_member", [f"DP4_{i}" for i in range(1, len(lib_df) + 1)])
    # 8 standard columns, then the metric columns, then the scoring/identity extras; missing cells blank.
    lib_df = lib_df.reindex(columns=["library_member"] + keep8 + list(METRICS.values()) + EXTRAS)
    assert (lib_df.sequence.str.len() <= 103).all(), "a sequence exceeds 103 residues!"
    assert lib_df.sequence.is_unique, "duplicate sequences remain after dedup"
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
