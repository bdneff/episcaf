#!/usr/bin/env python3
"""
stage07_named_peptides.py -- export the LadnerLab oligo-encoder input from the assembled DP4 library.

The encoder's step-1 `main` takes a `-i` file of `name,seq` with NO header (one peptide per line,
max line length 128; see episcaf_pipeline/oligo_encoding/README.md). This slices the two columns
`library_member,sequence` out of the 8-column library (scripts/stage06_assemble.py) and writes that
file, validating the format the encoder assumes (no header, unique names, no stray commas, every line
<= 128 chars, every sequence a clean 103-mer of standard amino acids).

Usage:
  python scripts/stage07_named_peptides.py \
      --library data/libraries/dp4_library.csv \
      --out data/libraries/dp4_named_peptides.csv
"""
from __future__ import annotations
import argparse
from pathlib import Path
import pandas as pd

AA = set("ACDEFGHIKLMNPQRSTVWY")
MAX_LINE = 128


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--library", default="data/libraries/dp4_library.csv")
    ap.add_argument("--out", default="data/libraries/dp4_named_peptides.csv")
    ap.add_argument("--sample", type=int, default=0,
                    help="if >0, emit only N evenly-spaced rows (deterministic) as a smoke-test input "
                         "that still spans every component; 0 = the whole library")
    args = ap.parse_args()

    d = pd.read_csv(args.library, low_memory=False)
    if args.sample > 0 and args.sample < len(d):
        step = len(d) / args.sample
        idx = [int(i * step) for i in range(args.sample)]   # evenly spaced -> touches all categories
        d = d.iloc[idx].reset_index(drop=True)
    name = d["library_member"].astype(str)
    seq = d["sequence"].astype(str)

    # validate the encoder's assumptions before writing (fail loudly, don't ship a bad input)
    assert name.is_unique, "library_member is not unique"
    assert not name.str.contains("[,\n]").any(), "a name contains a comma/newline"
    bad = seq[~seq.map(lambda s: set(s) <= AA)]
    assert bad.empty, f"{len(bad)} sequences have non-standard residues, e.g. {bad.iloc[0]!r}"
    lines = name + "," + seq
    too_long = lines[lines.str.len() > MAX_LINE]
    assert too_long.empty, f"{len(too_long)} lines exceed {MAX_LINE} chars (max {lines.str.len().max()})"

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    # header=False: the encoder reads name,seq positionally with no header row
    d[["library_member", "sequence"]].to_csv(args.out, index=False, header=False)

    print(f"[stage07] wrote {len(d):,} named peptides -> {args.out}")
    print(f"[stage07] name+seq line length: min {lines.str.len().min()}, max {lines.str.len().max()} "
          f"(limit {MAX_LINE})")
    print(f"[stage07] sequence length: {sorted(seq.str.len().unique())}")


if __name__ == "__main__":
    main()
