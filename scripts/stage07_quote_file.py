#!/usr/bin/env python3
"""Stage 07d -- emit the 2-column vendor QUOTE file (name,349mer) from the order file.

Twist (and Agilent) will quote off a plain two-column list of the oligos to be synthesized; the
sequences do NOT have to be final to get a price, so this can be sent while the library is still being
checked (John, 2026-07-23). This is a pure re-shaping of `dp4_order_file.csv` -- it never re-encodes and
never re-derives a sequence, so the quote and the order are the same molecules by construction.

The only real work is the NAME. The order file's `Seq ID` is `<library_member>_<encoding_id>`
(e.g. `DP4_1_00046`: library member DP4_1, encoding 00046). The requested quote format is a flat
zero-padded name (`DP4_00001`), so each row is renamed `DP4_%05d` from its library-member number --
a 1:1, order-preserving relabel, so DP4_00001 is library member DP4_1.

Checks before writing (all hard-fail): every oligo is `--length` nt, carries the 5'/3' Twist adapters,
names are unique and contiguous 1..N, and the row count matches the library.

Usage
-----
    python scripts/stage07_quote_file.py \
        --order-file data/libraries/dp4_order_file.csv \
        --out data/libraries/dp4_quote_file.csv

    # raw two columns with no header row (some vendor templates want this):
    python scripts/stage07_quote_file.py --no-header ...
"""
from __future__ import annotations
import argparse
import csv
import re
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]

# The DP4 20-mer Twist adapters -- the same pair stage07_order_file.py pins (see memory
# `oligo-adapter-trap`: the encoder master's own default is a 19-mer and is WRONG for DP4).
STD_PREFIX = "ACCTATACTTCCAAGGCGCA"
STD_SUFFIX = "GGTGACTCTCTGTCTTGGCT"

_SEQID = re.compile(r"^(?P<lib>DP4_(?P<n>\d+))_(?P<enc>\d+)$")


def load_order(path: Path):
    """(library_member_number, oligo) per row, in file order."""
    with path.open() as fh:
        rdr = csv.reader(fh)
        header = next(rdr)
        if len(header) < 2:
            sys.exit(f"ERROR: {path} is not a 2-column order file (header={header})")
        rows = []
        for ln, rec in enumerate(rdr, start=2):
            if not rec or not rec[0].strip():
                continue
            m = _SEQID.match(rec[0].strip())
            if m is None:
                sys.exit(f"ERROR: {path}:{ln} Seq ID {rec[0]!r} is not <library_member>_<encoding_id>")
            rows.append((int(m.group("n")), rec[1].strip()))
    return header, rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--order-file", default="data/libraries/dp4_order_file.csv")
    ap.add_argument("--out", default="data/libraries/dp4_quote_file.csv")
    ap.add_argument("--library", default="data/libraries/dp4_library.csv",
                    help="cross-check the row count against the shipped library")
    ap.add_argument("--length", type=int, default=349, help="expected oligo length (DP4: 349 nt)")
    ap.add_argument("--prefix", default=STD_PREFIX)
    ap.add_argument("--suffix", default=STD_SUFFIX)
    ap.add_argument("--no-header", action="store_true",
                    help="write bare `name,sequence` rows with no header line")
    args = ap.parse_args()

    order = _ROOT / args.order_file if not Path(args.order_file).is_absolute() else Path(args.order_file)
    out = _ROOT / args.out if not Path(args.out).is_absolute() else Path(args.out)

    _, rows = load_order(order)
    if not rows:
        sys.exit(f"ERROR: no data rows in {order}")

    # --- checks (hard-fail; a quote that disagrees with the order is worse than no quote) ---
    bad_len = [(n, len(s)) for n, s in rows if len(s) != args.length]
    if bad_len:
        sys.exit(f"ERROR: {len(bad_len)} oligos are not {args.length} nt, e.g. {bad_len[:3]}")
    bad_ad = [n for n, s in rows if not (s.startswith(args.prefix) and s.endswith(args.suffix))]
    if bad_ad:
        sys.exit(f"ERROR: {len(bad_ad)} oligos lack the 20-mer Twist adapters, e.g. {bad_ad[:3]}")

    nums = [n for n, _ in rows]
    if len(set(nums)) != len(nums):
        sys.exit(f"ERROR: duplicate library-member numbers in {order}")
    if sorted(nums) != list(range(1, len(nums) + 1)):
        sys.exit(f"ERROR: library-member numbers are not contiguous 1..{len(nums)} "
                 f"(min {min(nums)}, max {max(nums)}, n {len(nums)})")

    lib = _ROOT / args.library if not Path(args.library).is_absolute() else Path(args.library)
    if lib.exists():
        with lib.open() as fh:
            n_lib = sum(1 for _ in fh) - 1
        if n_lib != len(rows):
            sys.exit(f"ERROR: order file has {len(rows)} rows but {lib.name} has {n_lib}")

    width = max(5, len(str(len(rows))))
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        if not args.no_header:
            w.writerow(["name", "sequence"])
        for n, s in rows:
            w.writerow([f"DP4_{n:0{width}d}", s])

    print(f"[quote] wrote {len(rows):,} oligos ({args.length} nt, 20-mer adapters verified) -> {out}")
    print(f"[quote] names DP4_{1:0{width}d}..DP4_{len(rows):0{width}d}"
          f"{'  (no header row)' if args.no_header else '  (header: name,sequence)'}")


if __name__ == "__main__":
    main()
