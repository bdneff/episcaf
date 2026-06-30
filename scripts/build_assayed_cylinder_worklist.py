#!/usr/bin/env python3
"""Build the worklist of assayed scaffold designs to compute the native-aware cylinder on.

The cylinder values we have on disk are from the RFD3 reproduction run (different molecules);
the synthesized/assayed designs live in dp2.parquet and have no cylinder column. This emits the
list of designs whose cylinder we actually need -- the scaffoldedAbEpitope members that were
assayed -- keyed so the cluster run (scripts/assayed_native_cylinder.py) can find each design's
AF3 structure and so the result joins straight back to results/dp3_binding_metrics.csv.

Per design we carry: designedSequence (the join key to binding), the epitope id (-> native
<id>.pdb), and assay_scaffolded_epitope_id (the dp2 hash -> epitope residue indices + the key
under which the assayed AF3 outputs should be findable).

Run (local): python scripts/build_assayed_cylinder_worklist.py
"""
from __future__ import annotations
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from configs.paths import DP2_PARQUET_LOCAL  # noqa: E402

BIND = ROOT / "results" / "dp3_binding_metrics.csv"
OUT = ROOT / "results" / "assayed_cylinder_worklist.csv"


def main() -> None:
    b = pd.read_csv(BIND)
    scaf = b.loc[b.category == "scaffoldedAbEpitope",
                 ["designedSequence", "library_member", "Target", "cognate_ab", "is_pass"]].copy()

    dp2 = pd.read_parquet(DP2_PARQUET_LOCAL)
    key = (dp2[["scaffolded_epitope_seq", "assay_scaffolded_epitope_id", "id"]]
           .dropna(subset=["scaffolded_epitope_seq"])
           .drop_duplicates("scaffolded_epitope_seq"))
    key["assay_scaffolded_epitope_id"] = key.assay_scaffolded_epitope_id.astype(str).str.lower()

    m = scaf.merge(key, left_on="designedSequence", right_on="scaffolded_epitope_seq", how="left")
    missing = m.assay_scaffolded_epitope_id.isna().sum()
    assert missing == 0, f"{missing} assayed scaffolds have no dp2 assay hash"
    m = m.drop(columns="scaffolded_epitope_seq")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    m.to_csv(OUT, index=False)
    print(f"wrote {OUT}  ({len(m)} assayed scaffold designs)")
    print(f"  epitope ids (need <id>.pdb in native_dir): {sorted(m.id.unique())}")
    print(f"  with a usable cognate antibody: {m.cognate_ab.notna().sum()}")


if __name__ == "__main__":
    main()
