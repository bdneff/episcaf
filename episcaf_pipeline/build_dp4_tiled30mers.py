#!/usr/bin/env python3
"""
build_dp4_tiled30mers.py  --  regenerate the DP4 tiled-30mer linear-control library.

This wraps tile_fasta for the specific DP4 control set John asked for: all mAb-target
antigens (from dp2.parquet) plus the 3 tiled antigens (1D2K, 4WAT, 6M0J), chopped into
overlapping 30mers at step 6, in the constant assay construct (model "(none)").

It writes one combined library in DP2 annotated format, with library_member numbered
continuously (DP4_1..N) across both sources, and a per-antigen summary CSV
(target, len, mer, step, tiles) for spot-checking.

Notes / decisions baked in (flag to John if they change):
  - mAb antigens deduped by sequence -> 59 distinct antigens (the 2 duplicate antigens
    are kept; collapse with a name map if 57 is wanted).
  - target = dp2 id (a patch id like 2h32_0P); remap to protein names downstream if needed.
  - 6M0J is the 223-aa PDB-resolved region (control_antigens.fasta), matching the
    scaffolded 12mer set, not the full deposited FASTA.

Usage:
  python -m episcaf_pipeline.build_dp4_tiled30mers \
      --dp2 ../known_antigen/analysis/full_run/dp2.parquet \
      --antigens ../12mer_tiling/data/control_antigens.fasta \
      --out ../12mer_tiling/data/dp4_tiled30mers.csv
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

from episcaf_pipeline.tile_fasta import (
    read_parquet_records, read_fasta, tile_sequence, CONSTRUCT_PREFIX,
)

MER, STEP, CATEGORY, ID_PREFIX = 30, 6, "tiled30mer", "DP4"


def build(dp2: Path, antigens: Path):
    rows, seqlen = [], {}

    def add(rid: str, seq: str):
        seqlen[rid] = len(seq)
        for start, tile in tile_sequence(seq, MER, STEP, include_cterm=True):
            rows.append(dict(
                sequence=CONSTRUCT_PREFIX + tile, category=CATEGORY, model="(none)",
                designedSequence=tile, designedSequenceLength=len(tile),
                design_ID=start, target=rid,
            ))

    for rid, seq in read_parquet_records(dp2, "id", "antigen_seq", dedupe=True):
        add(rid, seq)
    for rid, seq in read_fasta(antigens):
        add(rid, seq)

    df = pd.DataFrame(rows)
    df.insert(0, "library_member", [f"{ID_PREFIX}_{i}" for i in range(1, len(df) + 1)])
    df = df[["library_member", "sequence", "category", "model", "designedSequence",
             "designedSequenceLength", "design_ID", "target"]]

    summary = (df.groupby("target").size().rename("tiles").reset_index())
    summary["len"] = summary["target"].map(seqlen)
    summary["mer"], summary["step"] = MER, STEP
    summary = summary[["target", "len", "mer", "step", "tiles"]].sort_values("target")
    return df, summary


def main():
    here = Path(__file__).resolve().parent.parent          # episcaf_v2/
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dp2", type=Path,
                    default=here.parent / "known_antigen/analysis/full_run/dp2.parquet")
    ap.add_argument("--antigens", type=Path,
                    default=here.parent / "12mer_tiling/data/control_antigens.fasta")
    ap.add_argument("--out", type=Path,
                    default=here.parent / "12mer_tiling/data/dp4_tiled30mers.csv")
    args = ap.parse_args()

    df, summary = build(args.dp2, args.antigens)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)
    summ_path = args.out.with_name(args.out.stem + "_summary.csv")
    summary.to_csv(summ_path, index=False)

    n_mab = summary[~summary.target.isin(["1d2k", "4wat", "6m0j"])].shape[0]
    print(f"{summary.shape[0]} antigens ({n_mab} mAb + 3 tiled) -> {len(df)} tiled-30mers")
    print(f"  library -> {args.out}")
    print(f"  summary -> {summ_path}")


if __name__ == "__main__":
    main()
