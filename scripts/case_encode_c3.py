#!/usr/bin/env python3
"""
case_encode_c3.py -- case-encoded scaffoldEPITOPE for C3 (polyclonal 12-mer tiles). LOCAL.

Unlike C1/C5 (whole-epitope, positions from dp2) and C2 (single-island, positions from the ledger),
the C3 selection already carries both `design_seq` (the 103-mer) and `epitope_seq` (the 12-mer), and
the 12-mer is a unique contiguous substring of the design (verified: 1 match in all 8,780). So we just
locate it and uppercase it -> case-encoded sequence (epitope UPPER, scaffold lower). C3 designs are
natively 103-mers (no 104->103 trim needed). Feeds assembly's designedSequence column.

Usage:
  python scripts/case_encode_c3.py \
      --selection results/dp4_C3_12mer_ranked.top20.csv \
      --out results/dp4_C3_scaffoldEPITOPE.csv
"""
from __future__ import annotations
import argparse, re
from pathlib import Path
import pandas as pd


def case_encode(seq: str, start: int, length: int) -> str:
    chars = [c.lower() for c in seq]
    for i in range(start, start + length):
        if 0 <= i < len(chars):
            chars[i] = chars[i].upper()
    return "".join(chars)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selection", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    d = pd.read_csv(args.selection, low_memory=False)
    rows, n_ok, n_bad = [], 0, 0
    for r in d.itertuples(index=False):
        seq, epi = str(getattr(r, "design_seq")), str(getattr(r, "epitope_seq"))
        if seq == "nan" or epi == "nan":
            n_bad += 1; continue
        i = seq.find(epi)
        if i < 0 or seq.count(epi) != 1:      # not found, or ambiguous (multiple matches)
            n_bad += 1
            rows.append(dict(token=getattr(r, "token", ""), id=getattr(r, "id", ""),
                             antigen=getattr(r, "antigen", ""), status="epitope_not_unique",
                             scaffoldEPITOPE="")); continue
        se = case_encode(seq, i, len(epi))
        n_ok += 1
        rows.append(dict(token=getattr(r, "token", ""), id=getattr(r, "id", ""),
                         antigen=getattr(r, "antigen", ""), target=getattr(r, "antigen", ""),
                         design_seq=seq, scaffoldEPITOPE=se,
                         n_islands=len(re.findall(r"[A-Z]+", se)), status="ok"))

    out = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(args.out, index=False)
    print(f"[c3] {len(d)} designs -> {n_ok} encoded | {n_bad} skipped (not-found/ambiguous)")
    if n_ok:
        ok = out[out.status == "ok"]
        print(f"[c3] island-count: {ok.n_islands.value_counts().to_dict()} (expect all 1); "
              f"design_seq len: {ok.design_seq.str.len().value_counts().to_dict()}")
    print(f"[c3] wrote {args.out}")


if __name__ == "__main__":
    main()
