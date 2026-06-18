#!/usr/bin/env python3
"""
build_dp4_tiled30mers.py  --  regenerate the DP4 tiled-30mer linear-control library.

This wraps tile_fasta for the specific DP4 control set John asked for: all mAb-target
antigens plus the 3 tiled antigens (1D2K, 4WAT, 6M0J), chopped into overlapping 30mers at
step 6, in the constant assay construct (model "(none)").

Inputs are the in-repo sequence files, so a clone reproduces the library with no external
dependency (the full dp2.parquet is NOT needed here -- the 59 mAb antigen sequences were
extracted from it once into data/sequences/dp3_mab_antigens.fasta; see data/README.md):
  data/sequences/dp3_mab_antigens.fasta   59 mAb antigens (target = dp2 patch id)
  data/sequences/control_antigens.fasta   3 tiled antigens (1d2k/4wat/6m0j)

Output (DP2 annotated format, library_member numbered continuously DP4_1..N across both
sources, mAb first), plus a per-antigen summary for spot-checking:
  data/libraries/dp4_tiled30mers.csv
  data/libraries/dp4_tiled30mers_summary.csv

Decisions baked in (flag to John if they change):
  - 59 mAb antigens (the 2 duplicate antigens are kept; collapse if 57 is wanted).
  - target = dp2 patch id (e.g. 2h32_0P); remap to protein names downstream if needed.
  - 6M0J is the 223-aa PDB-resolved region, matching the scaffolded 12mer set.

Usage:
  python -m episcaf_pipeline.build_dp4_tiled30mers          # uses in-repo defaults
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

from episcaf_pipeline.tile_fasta import read_fasta, tile_sequence, CONSTRUCT_PREFIX

MER, STEP, CATEGORY, ID_PREFIX = 30, 6, "tiled30mer", "DP4"
TILED_ANTIGENS = {"1d2k", "4wat", "6m0j"}


def build(fastas):
    """Tile each FASTA in order; library_member numbered continuously across all."""
    rows, seqlen = [], {}
    for fa in fastas:
        for rid, seq in read_fasta(fa):
            seqlen[rid] = len(seq)
            for start, tile in tile_sequence(seq, MER, STEP, include_cterm=True):
                rows.append(dict(
                    sequence=CONSTRUCT_PREFIX + tile, category=CATEGORY, model="(none)",
                    designedSequence=tile, designedSequenceLength=len(tile),
                    design_ID=start, target=rid,
                ))

    df = pd.DataFrame(rows)
    df.insert(0, "library_member", [f"{ID_PREFIX}_{i}" for i in range(1, len(df) + 1)])
    df = df[["library_member", "sequence", "category", "model", "designedSequence",
             "designedSequenceLength", "design_ID", "target"]]

    summary = df.groupby("target").size().rename("tiles").reset_index()
    summary["len"] = summary["target"].map(seqlen)
    summary["mer"], summary["step"] = MER, STEP
    summary = summary[["target", "len", "mer", "step", "tiles"]].sort_values("target")
    return df, summary


def main():
    data = Path(__file__).resolve().parent.parent / "data"      # episcaf_v2/data
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mab-fasta", type=Path,
                    default=data / "sequences/dp3_mab_antigens.fasta")
    ap.add_argument("--antigens", type=Path,
                    default=data / "sequences/control_antigens.fasta")
    ap.add_argument("--out", type=Path, default=data / "libraries/dp4_tiled30mers.csv")
    args = ap.parse_args()

    df, summary = build([args.mab_fasta, args.antigens])      # mAb first, then tiled antigens
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    summ_path = args.out.with_name(args.out.stem + "_summary.csv")
    summary.to_csv(summ_path, index=False)

    n_mab = summary[~summary.target.isin(TILED_ANTIGENS)].shape[0]
    print(f"{summary.shape[0]} antigens ({n_mab} mAb + {summary.shape[0]-n_mab} tiled) "
          f"-> {len(df)} tiled-30mers")
    print(f"  library -> {args.out}")
    print(f"  summary -> {summ_path}")


if __name__ == "__main__":
    main()
