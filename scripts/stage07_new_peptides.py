#!/usr/bin/env python3
"""Stage 07b2 -- emit the peptides that still need oligo encoding (a top-up encoder input).

When the library grows (a deeper per-group top-k, a new arm), the peptides already encoded do NOT need
re-encoding -- only the new ones do. This writes exactly those, in the encoder's input format, so a
top-up run costs 38 peptides instead of 36,000.

`--encoded` accepts either form of "what is already encoded":
  * a step-2 encodings table (has an `AA Peptide` column)  -- the authoritative source, or
  * a `name,seq` peptide file (no header)                  -- e.g. a previous dp4_named_peptides.csv.
Matching is by SEQUENCE, never by name: a top-up renumbers `library_member` for everything after the
insertion point (adding 38 8VDL rows shifted all 21,759 minibinders), so names are not stable but the
molecules are.

LINE ENDINGS ARE LOAD-BEARING. The LadnerLab encoder splits on newline and takes the rest of the line
as the peptide, so a CRLF file leaves a trailing \\r on every sequence and the encoder rejects all of
them ("Warning: Line N: ... is invalid and will be skipped", then "Processed 0 lines"). That is exactly
how the 2026-07-23 top-up failed -- a csv.writer wrote \\r\\n (the default `excel` dialect). This script
writes bare \\n and hard-fails if a CR ever reaches the output.

Usage
-----
    python scripts/stage07_new_peptides.py \
        --peptides data/libraries/dp4_named_peptides.csv \
        --encoded  $WS/runs/dp4_encoding_full/DP4_best_encodings \
        --out      data/libraries/dp4_named_peptides.new.csv
"""
from __future__ import annotations
import argparse
import csv
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_peptides(path: Path):
    """`name,seq` with no header -> [(name, seq)]."""
    out = []
    with path.open(newline="") as fh:
        for ln, rec in enumerate(csv.reader(fh), start=1):
            if not rec or not rec[0].strip():
                continue
            if len(rec) < 2:
                sys.exit(f"ERROR: {path}:{ln} is not `name,seq`: {rec!r}")
            out.append((rec[0].strip(), rec[1].strip()))
    return out


def read_encoded(path: Path):
    """Already-encoded peptide SEQUENCES, from a step-2 encodings table or a name,seq file."""
    with path.open(newline="") as fh:
        first = fh.readline()
    if "AA Peptide" in first:
        with path.open(newline="") as fh:
            rdr = csv.DictReader(fh)
            seqs = {(r.get("AA Peptide") or "").strip() for r in rdr}
        seqs.discard("")
        return seqs
    return {s for _, s in read_peptides(path)}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--peptides", default="data/libraries/dp4_named_peptides.csv",
                    help="the CURRENT full encoder input (name,seq, no header)")
    ap.add_argument("--encoded", required=True,
                    help="what is already encoded: a step-2 encodings table (`AA Peptide` column) "
                         "or a previous name,seq peptide file")
    ap.add_argument("--out", default="data/libraries/dp4_named_peptides.new.csv")
    ap.add_argument("--length", type=int, default=103, help="expected peptide length (DP4: 103)")
    args = ap.parse_args()

    res = lambda p: Path(p) if Path(p).is_absolute() else _ROOT / p       # noqa: E731
    peptides = read_peptides(res(args.peptides))
    encoded = read_encoded(res(args.encoded))
    out = res(args.out)

    new = [(n, s) for n, s in peptides if s not in encoded]

    bad_len = [n for n, s in new if len(s) != args.length]
    if bad_len:
        sys.exit(f"ERROR: {len(bad_len)} new peptides are not {args.length} aa, e.g. {bad_len[:3]}")
    bad_aa = [n for n, s in new if set(s) - AA]
    if bad_aa:
        sys.exit(f"ERROR: {len(bad_aa)} new peptides carry non-standard residues, e.g. {bad_aa[:3]}")
    dupes = len(new) - len({s for _, s in new})
    if dupes:
        sys.exit(f"ERROR: {dupes} duplicate sequences among the new peptides")

    out.parent.mkdir(parents=True, exist_ok=True)
    # Written by hand, NOT csv.writer: the encoder needs bare \n (see the module docstring).
    with out.open("w", newline="\n") as fh:
        for n, s in new:
            fh.write(f"{n},{s}\n")

    raw = out.read_bytes()
    if b"\r" in raw:
        sys.exit(f"ERROR: {out} contains CR -- the encoder would reject every line")

    print(f"[new-peptides] {len(peptides):,} current - {len(encoded):,} already encoded "
          f"-> {len(new):,} need encoding -> {out}")
    if new:
        w = max(len(n) for n, _ in new)
        print(f"[new-peptides] line length {len(new[0][0])+1+args.length} "
              f"(name width {w}), LF endings verified, all {args.length} aa")


if __name__ == "__main__":
    main()
