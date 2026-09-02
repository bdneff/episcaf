#!/usr/bin/env python3
"""Flag oligos where the 5' primer's 3' end recurs INTERNALLY (mispriming risk).

Background (John Altin, 2026-08-05, post-order QC of DP4). The first primer is the 20-mer 5'
Twist adapter `ACCTATACTTCCAAGGCGCA`. If the tail of that primer also appears *inside* an oligo,
the primer can anneal internally instead of terminally during amplification -> a truncated
product. In DP4 these internal sites arise because the **tiled30mer (C4)** arm encodes an internal
TEV site, and the 5' primer also encodes TEV, so a stochastic codon choice can make an internal
stretch match the primer's 3' end.

This is NOT a fix for a shipped library (DP4 is ordered; John rates the impact low -- 10 bp is
likely too short to prime, 15 bp is borderline). It is the reproducible check behind that finding
and the seed for a **primer-similarity filter at the encoding step** for future designs: reject /
re-draw any candidate encoding whose payload contains a >=K-bp suffix of either primer.

Verified 2026-08-05 against `data/libraries/dp4_order_file.csv` -- reproduces John's counts exactly:
  internal last-10bp (CCAAGGCGCA):      52 / 36000  (34 tiled30mer, 17 minibinder, 1 epitopeControl)
  internal last-15bp (TACTTCCAAGGCGCA): 29 / 36000  (all 29 tiled30mer)  <- the higher-risk set
The 15-bp (riskier) sites are entirely the TEV-encoding C4 arm, confirming the mechanism.

Usage
-----
    python scripts/check_primer_internal_sites.py                 # summary + per-category counts
    python scripts/check_primer_internal_sites.py --min-len 12    # sweep the risk threshold
    python scripts/check_primer_internal_sites.py --list          # print flagged Seq IDs
"""
from __future__ import annotations
import argparse
import collections
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ORDER = ROOT / "data" / "libraries" / "dp4_order_file.csv"
LIBRARY = ROOT / "data" / "libraries" / "dp4_library.csv"

# The two Twist adapters / primers (see stage07_order_file.py, memory `oligo-adapter-trap`).
PRIMER5 = "ACCTATACTTCCAAGGCGCA"   # first primer -- the one John flagged
PRIMER3 = "GGTGACTCTCTGTCTTGGCT"   # second primer (3'); checked too, for completeness


def total_occ(s: str, sub: str) -> int:
    n = start = 0
    while (i := s.find(sub, start)) >= 0:
        n += 1
        start = i + 1
    return n


def category_map() -> dict[str, str]:
    cat: dict[str, str] = {}
    with LIBRARY.open() as f:
        for row in csv.DictReader(f):
            cat[row["library_member"]] = row["category"]
    return cat


def scan(min_len: int, primer: str, primer_at_5prime: bool):
    """Return (flagged Seq IDs, category counter) for oligos whose payload contains the
    last `min_len` bp of `primer` as an INTERNAL occurrence (beyond the terminal adapter copy)."""
    probe = primer[-min_len:]
    cat = category_map()
    flagged, by_cat = [], collections.Counter()
    n = 0
    with ORDER.open() as f:
        r = csv.reader(f)
        next(r)
        for seqid, oligo in r:
            n += 1
            # The adapter itself contains exactly one terminal copy of any of its own suffixes;
            # every occurrence beyond that first one is an internal (mispriming) site.
            if total_occ(oligo, probe) > 1:
                lm = seqid.rsplit("_", 1)[0]          # DP4_x_00046 -> DP4_x
                c = cat.get(lm, "?")
                flagged.append((seqid, c))
                by_cat[c] += 1
    return n, probe, flagged, by_cat


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--min-len", type=int, default=None,
                    help="single suffix length to check; default reports both 10 and 15 (John's)")
    ap.add_argument("--primer", choices=["5", "3"], default="5",
                    help="which primer's suffix to look for internally (default 5', the flagged one)")
    ap.add_argument("--list", action="store_true", help="print each flagged Seq ID")
    args = ap.parse_args()

    primer = PRIMER5 if args.primer == "5" else PRIMER3
    lengths = [args.min_len] if args.min_len else [10, 15]

    for L in lengths:
        n, probe, flagged, by_cat = scan(L, primer, args.primer == "5")
        print(f"\nprimer {args.primer}' last-{L}bp ({probe}): "
              f"{len(flagged)} / {n} oligos have an internal site")
        for c, v in by_cat.most_common():
            print(f"    {v:4d}  {c}")
        if args.list:
            for seqid, c in flagged:
                print(f"      {seqid}\t{c}")


if __name__ == "__main__":
    main()
