#!/usr/bin/env python3
"""
01_generate_contigs.py -- RFD3 contigs for the 8VDL PfEMP1 conserved-epitope scaffolding.

Two targets on 8VDL chain C (see ../README.md for the biology):
  epitope20 : the whole contiguous epitope C651-670 (20 aa, one island)
  hotspots  : only the functional hotspots F655, F656, E666 (two islands: 655-656 and 666)

Each target wraps its fixed motif in a constant TOTAL_LEN (=103) construct. The scaffold budget
(103 - motif residues) is split -- randomly, seeded -- across the k+1 scaffold segments (N-flank,
any inter-island gaps, C-flank), each >= MIN_SCAFFOLD. The fixed atoms carry the residues' native
crystal coordinates, so the motif's 3-D geometry is preserved regardless of construct spacing.

Output: one contigs CSV consumed by 02_emit_rfd3_inputs.py.

Usage:
  python scripts/01_generate_contigs.py --target epitope20 --n-contigs 500 --out 01_contigs/epitope20.csv
  python scripts/01_generate_contigs.py --target hotspots  --n-contigs 500 --out 01_contigs/hotspots.csv
"""
from __future__ import annotations
import argparse
import random
from pathlib import Path

import pandas as pd

TOTAL_LEN    = 103          # PepSeq maximum construct length
MIN_SCAFFOLD = 10           # minimum residues per scaffold segment (flanks and gaps)
CHAIN        = "C"
PDB_PATH     = "data/8VDL.pdb"

# Fixed motif per target: a list of (start_resid, end_resid) islands on chain C, in native order.
TARGETS = {
    "epitope20": [(651, 670)],            # the whole conserved EPCR-binding epitope
    "hotspots":  [(655, 656), (666, 666)],  # F655/F656 (adjacent) + E666
    "contact":   [(652, 653), (655, 657), (659, 661), (666, 667), (669, 670), (673, 673)],
    # ^ the AbDb/IEDB 4A contact epitope (13 residues, 6 islands; contact_epitope.py). Natively
    #   close-packed, so generate with --native-gaps (hold native inter-island spacing, vary flanks).
}


def random_composition(total: int, parts: int, min_each: int, rng: random.Random) -> list[int]:
    """`parts` positive ints, each >= min_each, summing to `total` (uniform via stars-and-bars)."""
    r = total - min_each * parts
    if r < 0:
        raise ValueError(f"cannot place {parts} segments of >= {min_each} in budget {total}")
    cuts = sorted(rng.randint(0, r) for _ in range(parts - 1))
    segs, prev = [], 0
    for c in cuts:
        segs.append(c - prev); prev = c
    segs.append(r - prev)
    return [s + min_each for s in segs]


def build_contig(islands: list[tuple[int, int]], segs: list[int]) -> str:
    """Interleave scaffold segments with fixed islands -> RFD3 comma-style contig string."""
    parts = [f"{segs[0]}-{segs[0]}"]
    for i, (a, b) in enumerate(islands):
        parts.append(f"{CHAIN}{a}-{b}")
        parts.append(f"{segs[i + 1]}-{segs[i + 1]}")
    return ",".join(parts)


def contig_len(islands: list[tuple[int, int]], segs: list[int]) -> int:
    return sum(segs) + sum(b - a + 1 for a, b in islands)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--target", required=True, choices=sorted(TARGETS))
    ap.add_argument("--n-contigs", type=int, default=500, help="John's 'n designs' knob (contigs)")
    ap.add_argument("--total-len", type=int, default=TOTAL_LEN)
    ap.add_argument("--min-scaffold", type=int, default=MIN_SCAFFOLD)
    ap.add_argument("--native-gaps", action="store_true",
                    help="hold the NATIVE inter-island spacing and randomize only the N/C flanks "
                         "(for close-packed multi-island motifs like contact, where forcing large "
                         "random gaps between natively-adjacent islands makes strained scaffolds)")
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()

    islands = TARGETS[args.target]
    motif_len = sum(b - a + 1 for a, b in islands)
    fixed_resids = sorted(r for a, b in islands for r in range(a, b + 1))
    rng = random.Random(args.seed)

    # native inter-island gaps (residues between consecutive islands in the crystal)
    native_gaps = [islands[i + 1][0] - islands[i][1] - 1 for i in range(len(islands) - 1)]
    if args.native_gaps:
        flank_budget = args.total_len - motif_len - sum(native_gaps)   # split across N/C flanks only
    else:
        budget = args.total_len - motif_len
        n_seg = len(islands) + 1                    # N-flank, gaps, C-flank (all randomized)

    seen, rows = set(), []
    tries = 0
    while len(rows) < args.n_contigs and tries < args.n_contigs * 50:
        tries += 1
        if args.native_gaps:
            n_flank, c_flank = random_composition(flank_budget, 2, args.min_scaffold, rng)
            segs = [n_flank, *native_gaps, c_flank]
        else:
            segs = random_composition(budget, n_seg, args.min_scaffold, rng)
        contig = build_contig(islands, segs)
        if contig in seen:
            continue
        seen.add(contig)
        assert contig_len(islands, segs) == args.total_len
        rows.append({
            "contig_id":     len(rows),
            "target":        args.target,
            "chain":         CHAIN,
            "islands":       ";".join(f"{a}-{b}" for a, b in islands),
            "fixed_resids":  ",".join(map(str, fixed_resids)),
            "scaffold_segs": ",".join(map(str, segs)),
            "contig_string": contig,
            "total_len":     args.total_len,
            "input_pdb":     PDB_PATH,   # relative (data/8VDL.pdb); resolved at emit time, portable
        })

    df = pd.DataFrame(rows)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out, index=False)

    gap_mode = (f"native gaps {native_gaps}, flanks randomized" if args.native_gaps
                else "all segments randomized")
    print(f"target={args.target}  islands={islands}  motif={motif_len} aa  ({len(islands)} islands)")
    print(f"  spacing: {gap_mode}")
    print(f"wrote {len(df)} unique contigs -> {args.out}"
          + ("" if len(df) == args.n_contigs else f"  (requested {args.n_contigs}; budget-limited)"))
    print(f"  designs = {len(df)} contigs x 8 RFD3 x 8 MPNN = {len(df) * 64:,}")
    print("  sample:")
    print(df[["contig_id", "scaffold_segs", "contig_string"]].head(4).to_string(index=False))


if __name__ == "__main__":
    main()
