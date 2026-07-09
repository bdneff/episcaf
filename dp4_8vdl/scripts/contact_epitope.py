#!/usr/bin/env python3
"""
contact_epitope.py -- the AbDb/IEDB "contact" epitope definition: every antigen residue with a
heavy (non-H) atom within CUTOFF (default 4.0 A) of any antibody heavy atom.

For 8VDL this flags the chain-C residues that actually touch the C7 Fab (chains H/L) -- a third
epitope definition to scaffold, alongside the contiguous 651-670 window (epitope20) and the
paper's F655/F656/E666 hotspots (hotspots). Run this first just to SEE which residues it flags.

Usage:
  python dp4_8vdl/scripts/contact_epitope.py                       # 8VDL C vs H,L at 4.0 A
  python dp4_8vdl/scripts/contact_epitope.py --cutoff 5 --pdb ... --antigen-chain A --antibody-chains B,C
"""
from __future__ import annotations
import argparse
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


AA3 = {"ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS",
       "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"}


def is_heavy(atom_name: str, element: str) -> bool:
    el = element.strip()
    if el:
        return el not in ("H", "D")
    return not atom_name.strip().lstrip("0123456789").startswith("H")


def parse_pdb(pdb: Path):
    """-> dict chain -> list of (resseq:int, icode, resname, atom_name, xyz). Heavy atoms only."""
    atoms: dict[str, list] = {}
    for line in pdb.read_text().splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        name = line[12:16]
        element = line[76:78] if len(line) >= 78 else ""
        if not is_heavy(name, element):
            continue
        chain = line[21]
        resseq = int(line[22:26])
        icode = line[26]
        resname = line[17:20].strip()
        xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        atoms.setdefault(chain, []).append((resseq, icode, resname, name.strip(), xyz))
    return atoms


def contiguous_runs(resids: list[int]) -> list[tuple[int, int]]:
    runs, s = [], None
    for i, r in enumerate(resids):
        if s is None:
            s = prev = r
        elif r == prev + 1:
            prev = r
        else:
            runs.append((s, prev)); s = prev = r
    if s is not None:
        runs.append((s, prev))
    return runs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pdb", type=Path, default=Path("dp4_8vdl/data/8VDL.pdb"))
    ap.add_argument("--antigen-chain", default="C")
    ap.add_argument("--antibody-chains", default="H,L")
    ap.add_argument("--cutoff", type=float, default=4.0)
    args = ap.parse_args()

    atoms = parse_pdb(args.pdb)
    ab_chains = [c.strip() for c in args.antibody_chains.split(",")]
    # amino-acid heavy atoms only (drop waters / hetero groups) on both sides
    ab_xyz = np.array([a[4] for c in ab_chains for a in atoms.get(c, []) if a[2] in AA3], float)
    if not len(ab_xyz):
        raise SystemExit(f"no antibody heavy atoms in chains {ab_chains}")
    tree = cKDTree(ab_xyz)

    # min antibody distance per antigen residue (amino acids only)
    per_res: dict[tuple, list] = {}
    for a in atoms.get(args.antigen_chain, []):
        if a[2] not in AA3:
            continue
        per_res.setdefault((a[0], a[1], a[2]), []).append(a[4])
    flagged = []
    for (resseq, icode, resname), coords in sorted(per_res.items()):
        dmin, _ = tree.query(np.asarray(coords, float), k=1)
        if dmin.min() <= args.cutoff:
            flagged.append((resseq, resname, float(dmin.min())))

    resids = [r for r, _, _ in flagged]
    print(f"8VDL contact epitope: chain {args.antigen_chain} residues with a heavy atom "
          f"<= {args.cutoff} A of chains {ab_chains}")
    print(f"  {len(flagged)} residues flagged\n")
    seq1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G",
            "HIS":"H","ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S",
            "THR":"T","TRP":"W","TYR":"Y","VAL":"V"}
    for resseq, resname, d in flagged:
        aa = seq1.get(resname, "?")
        marks = []
        if 651 <= resseq <= 670: marks.append("in 651-670")
        if resseq in (655, 656, 666): marks.append("HOTSPOT")
        print(f"  {resname}{resseq} ({aa})  min_dist={d:.2f}  {'  '.join(marks)}")

    print(f"\n  contiguous runs (islands): "
          + ", ".join(f"{a}-{b}" for a, b in contiguous_runs(resids)))
    inwin = [r for r in resids if 651 <= r <= 670]
    print(f"  overlap with 651-670 window: {len(inwin)}/{len(resids)} flagged residues; "
          f"window residues NOT flagged: {sorted(set(range(651,671)) - set(resids))}")
    print(f"  hotspots F655/F656/E666 flagged: "
          f"{[r for r in (655,656,666) if r in resids]}")


if __name__ == "__main__":
    main()
