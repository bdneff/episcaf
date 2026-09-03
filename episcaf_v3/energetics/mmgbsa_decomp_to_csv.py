#!/usr/bin/env python3
"""Per-residue MM-GBSA binding contribution for the antigen, parsed straight from sander TDC output.

gmx_MMPBSA's own decomp writer crashes on our -merge all topology (H/L chains share residue numbers),
but the energies it computed are fine -- in the sander mdout `TDC` lines. This parses those directly.

TDC line: `TDC <resnum> <internal> <vdW> <electrostatic> <polar-solv(EGB)> <nonpolar-solv(surf)>`,
one per decomposed residue per frame. Per-residue binding contribution (single-trajectory MM-GBSA):
    dG_i = <TDC_total_i(complex)>_frames - <TDC_total_i(ligand)>_frames
The internal term is identical in complex and ligand (same coordinates), so it cancels; what remains
is residue i's van der Waals + electrostatic interaction with the ANTIBODY plus its desolvation
(the EGB change) -- i.e. the desolvation-inclusive per-residue binding energy. dG < 0 is favorable.

Antigen residues: in the complex (merged 1-558) they are 430-558 -> lysozyme = resnum-429; in the
ligand-alone file they are already lysozyme-numbered (1-129). We match on lysozyme number.

Run:
    python mmgbsa_decomp_to_csv.py --complex md/3hfm/complex_tdc.txt --ligand md/3hfm/ligand_tdc.txt \
        --resnames md/3hfm/holo_ie_mean.csv --out md/3hfm/mmgbsa_perres.csv
"""
import argparse
from collections import defaultdict
import csv
from pathlib import Path


def read_tdc(path, antigen_offset):
    """resnum->list of per-frame (int, vdw, eel, egb, surf). If antigen_offset given, keep only
    residues > that offset (the antigen in the complex) and shift to lysozyme numbering."""
    per = defaultdict(list)
    for line in Path(path).read_text().splitlines():
        p = line.split()
        if not p or p[0] != "TDC":
            continue
        rn = int(p[1])
        vals = [float(x) for x in p[2:7]]  # int, vdw, eel, egb(polar), surf(nonpolar)
        if antigen_offset is not None:
            if rn <= antigen_offset:      # antibody residue in the complex; skip
                continue
            rn -= antigen_offset          # -> lysozyme numbering
        per[rn].append(vals)
    return per


def mean_terms(per):
    """resnum -> mean [int, vdw, eel, egb, surf] over frames."""
    out = {}
    for rn, frames in per.items():
        n = len(frames)
        out[rn] = [sum(f[k] for f in frames) / n for k in range(5)]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--complex", required=True, help="complex_tdc.txt (grep '^TDC' of the complex mdout)")
    ap.add_argument("--ligand", required=True, help="ligand_tdc.txt (grep '^TDC' of the ligand mdout)")
    ap.add_argument("--offset", type=int, default=429,
                    help="antigen residues in the complex are numbered > offset; lysozyme = resnum-offset")
    ap.add_argument("--resnames", default=None, help="a csv with ag_res_idx,resname to attach names")
    ap.add_argument("--out", default="mmgbsa_perres.csv")
    args = ap.parse_args()

    comp = mean_terms(read_tdc(args.complex, antigen_offset=args.offset))
    lig = mean_terms(read_tdc(args.ligand, antigen_offset=None))
    print(f"complex antigen residues: {len(comp)}   ligand residues: {len(lig)}")
    shared = sorted(set(comp) & set(lig))
    missing = sorted((set(comp) | set(lig)) - set(shared))
    if missing:
        print(f"  WARNING: residues not in both files: {missing}")
    print(f"  matched: {len(shared)}")

    names = {}
    if args.resnames:
        for r in csv.DictReader(open(args.resnames)):
            names[int(r["ag_res_idx"])] = r["resname"]

    # The 'sas' column is SASA-like (Å^2), not a kcal/mol energy (per-residue values reach ~100);
    # the nonpolar solvation energy is surften * SASA. gmx_MMPBSA's default surften is 0.0072.
    GAMMA = 0.0072
    rows = []
    for rn in shared:
        d = [comp[rn][k] - lig[rn][k] for k in range(5)]   # complex - ligand, per term
        dvdw, deel, dpol = d[1], d[2], d[3]                 # d[0] (internal) cancels ~0
        dnonpol = GAMMA * d[4]                               # scale the SASA-like column to energy
        mmgbsa_dg = dvdw + deel + dpol + dnonpol             # per-residue binding contribution
        rows.append((rn, names.get(rn, "?"), round(d[0], 3), round(dvdw, 3), round(deel, 3),
                     round(dpol, 3), round(dnonpol, 3), round(mmgbsa_dg, 3)))

    with open(args.out, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["ag_res_idx", "resname", "d_int", "d_vdw", "d_eel", "d_polar", "d_nonpolar", "mmgbsa_dg"])
        w.writerows(rows)
    print(f"wrote {args.out}  ({len(rows)} antigen residues)")

    # quick look: most favorable (most negative dG)
    print("\nmost favorable antigen residues by MM-GBSA dG (kcal/mol):")
    for rn, nm, di, dv, de, dp, dn, dg in sorted(rows, key=lambda r: r[-1])[:10]:
        print(f"  {nm}{rn:<4d} dG={dg:8.2f}   (vdw {dv:7.2f}  eel {de:8.2f}  polar {dp:8.2f}  nonpol {dn:6.2f})")


if __name__ == "__main__":
    main()
