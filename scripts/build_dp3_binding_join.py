#!/usr/bin/env python3
"""Join DP3 assay binding data to the design metrics in dp2.parquet.

John shipped DP3 PepSeq binding intensities for 8 antibodies, split across two assay
runs (2026-06): IM226 (6 Abs: 6o9i 8cz8 8jnk 6xxv 5fhx 7ox3) and IM229 (2 Abs: 8db4
8pww). 4xwo and 7a3t were dropped (low yield / epitope too small). The two CSVs hold the
SAME 1000 DP2 library members and identical design columns; they differ only in their
NoAb baseline and per-Ab intensity columns. See data/dp3_binding/README.md.

This is the bridge from "what binds" to "what the design looks like": every
`scaffoldedAbEpitope` member carries the exact scaffold sequence it was synthesized from,
which matches `scaffolded_epitope_seq` in the parquet 1-to-1 (verified 403/403). We attach
the AlphaFold3 design metrics (overall_rmsd, epitope_chunk_rmsd_vs_mpnn, mean_pae,
af3_n_clash_res, and the inherited 4-filter is_pass) so binding can be regressed on them.

For each design we also resolve its COGNATE antibody: the Ab raised against the same PDB
as the design's Target (e.g. Target 7ox3_0P -> Ab 7ox3), pulling that Ab's intensity and
its run's NoAb baseline, plus John's log10(1+x) transform and the log-enrichment
(log10(1+Ab) - log10(1+NoAb)) he uses as the binding readout.

Output: results/dp3_binding_metrics.csv (one row per library member; all raw assay columns
kept, design metrics filled in for the scaffolded designs).

Run (local dev, off-cluster):
    python scripts/build_dp3_binding_join.py
"""
from __future__ import annotations
import re, sys
from pathlib import Path
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from configs.paths import DP2_PARQUET_LOCAL  # noqa: E402

BINDING_DIR = ROOT / "data" / "dp3_binding"
IM226 = BINDING_DIR / "20260114_IM226_DP2_mAbs.csv"
IM229 = BINDING_DIR / "20260202_IM229_DP2.csv"
OUT = ROOT / "results" / "dp3_binding_metrics.csv"

# shared design columns (identical across the two runs)
DESIGN_COLS = ["library_member", "sequence", "category", "designedSequence",
               "designedSequenceLength", "antibodySequence", "Model", "Design_ID",
               "Target", "Kd"]

# 4-filter thresholds (CLAUDE.md "Fixed facts")
def is_pass(d: pd.DataFrame) -> pd.Series:
    return ((d.overall_rmsd <= 2) & (d.epitope_chunk_rmsd_vs_mpnn <= 1)
            & (d.mean_pae < 5) & (d.af3_n_clash_res == 0))

log1p10 = lambda x: np.log10(1.0 + x)  # John's transform


def ab_columns(df: pd.DataFrame) -> dict[str, str]:
    """Map antibody PDB id -> its intensity column. Ab cols look like
    'X7ox3_1X_0.1pmol_beadsFirst'; strip the leading X, take the PDB before the first _."""
    out = {}
    for c in df.columns:
        m = re.match(r"^X([0-9a-z]{4})_", c)
        if m:
            out[m.group(1)] = c
    return out


def noab_column(df: pd.DataFrame) -> str:
    cols = [c for c in df.columns if c.startswith("NoAb")]
    assert len(cols) == 1, f"expected one NoAb col, got {cols}"
    return cols[0]


def main() -> None:
    a = pd.read_csv(IM226)
    b = pd.read_csv(IM229)
    assert set(a.library_member) == set(b.library_member), "library_member sets differ"

    # per-run assay columns (NoAb + Abs), keyed back to library_member
    a_assay = [noab_column(a)] + list(ab_columns(a).values())
    b_assay = [noab_column(b)] + list(ab_columns(b).values())
    merged = a[DESIGN_COLS + a_assay].merge(
        b[["library_member"] + b_assay], on="library_member", how="outer")
    assert len(merged) == len(a), "merge changed row count"

    # antibody -> (intensity col, this run's NoAb col)
    ab_map: dict[str, tuple[str, str]] = {}
    for src in (a, b):
        noab = noab_column(src)
        for pdb, col in ab_columns(src).items():
            ab_map[pdb] = (col, noab)
    print(f"antibodies: {sorted(ab_map)}")

    # --- design metrics from the parquet, joined on scaffold sequence ---
    df = pd.read_parquet(DP2_PARQUET_LOCAL)
    df["af3_n_clash_res"] = df.af3_clash_resindices.apply(
        lambda x: 0 if x is None else len(x))
    metric_cols = ["overall_rmsd", "epitope_chunk_rmsd_vs_mpnn", "mean_pae", "af3_n_clash_res"]
    pq = (df[["scaffolded_epitope_seq"] + metric_cols]
          .dropna(subset=["scaffolded_epitope_seq"])
          .drop_duplicates("scaffolded_epitope_seq")
          .rename(columns={"scaffolded_epitope_seq": "designedSequence"}))

    scaf_mask = merged.category == "scaffoldedAbEpitope"
    n_scaf = int(scaf_mask.sum())
    out = merged.merge(pq, on="designedSequence", how="left")
    matched = out.loc[scaf_mask, "overall_rmsd"].notna().sum()
    print(f"scaffoldedAbEpitope: {matched}/{n_scaf} joined to parquet metrics")
    assert matched == n_scaf, "not every scaffold design matched a parquet metric row"
    out["is_pass"] = is_pass(out)

    # --- cognate antibody per design ---
    def cognate_pdb(target: str) -> str | None:
        pdb = str(target).split("_")[0].lower()
        return pdb if pdb in ab_map else None

    out["cognate_ab"] = out.Target.map(cognate_pdb)
    cog_signal, cog_noab = [], []
    for _, r in out.iterrows():
        ab = r.cognate_ab
        if ab is None:
            cog_signal.append(np.nan); cog_noab.append(np.nan); continue
        ab_col, noab_col = ab_map[ab]
        cog_signal.append(r[ab_col]); cog_noab.append(r[noab_col])
    out["cognate_ab_signal"] = cog_signal
    out["cognate_noab"] = cog_noab
    out["cognate_log_ab"] = log1p10(out.cognate_ab_signal)
    out["cognate_log_noab"] = log1p10(out.cognate_noab)
    out["cognate_log_enrichment"] = out.cognate_log_ab - out.cognate_log_noab

    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    print(f"\nwrote {OUT}  ({len(out)} rows, {len(out.columns)} cols)")

    # quick orientation: cognate enrichment of scaffold designs by Ab
    s = out[scaf_mask & out.cognate_ab.notna()]
    print(f"\nscaffold designs with a cognate Ab: {len(s)}")
    g = (s.groupby("cognate_ab")
           .agg(n=("library_member", "size"),
                n_pass=("is_pass", "sum"),
                med_enrich=("cognate_log_enrichment", "median"))
           .sort_values("n", ascending=False))
    print(g.to_string())


if __name__ == "__main__":
    main()
