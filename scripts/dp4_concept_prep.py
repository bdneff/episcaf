#!/usr/bin/env python3
"""Prepare real-structure PDBs for the two "why it matters" concept renders (diagnostics, vaccine).

Source of truth: 8PWW = PfRH5 (chain A) bound to the scFv antibody MAD8-151 (chain B). The DP4
designs graft the RH5 epitope (two motifs, AVDAFIKKINEA + ICMDMKNYGTNLFE) onto a scaffold. We
superpose each design onto native chain A by the epitope CA atoms, so the design's scaffold lands
in the native frame -- the same frame the real scFv already occupies. That lets a single molecular
render show the real antibody sitting on the scaffolded epitope (diagnostics), and the conserved
epitope as a small patch on the whole antigen vs. isolated on the scaffold (vaccine).

Writes PDBs into data/dp4_binding/concept/ (gitignored, under data/). Rendered by
scripts/render_concept_figs.sh.

Run:  /usr/bin/python3 scripts/dp4_concept_prep.py
"""
from pathlib import Path
import gemmi

ROOT = Path(__file__).resolve().parents[1]
NATIVE = ROOT/"data/dp4_binding/structures/8PWW.cif"
DES = {"low": ROOT/"data/dp4_binding/designs/low_rmsd0.22.cif",
       "high": ROOT/"data/dp4_binding/designs/high_rmsd3.17.cif"}
OUT = ROOT/"data/dp4_binding/concept"; OUT.mkdir(parents=True, exist_ok=True)
MOTIFS = ["AVDAFIKKINEA", "ICMDMKNYGTNLFE"]      # the two grafted epitope segments
EPI_NATIVE = [(205, 216), (279, 292)]            # their native resid ranges on chain A

def one(r): return gemmi.find_tabulated_residue(r.name).one_letter_code.upper()

def chain(path, name):
    st = gemmi.read_structure(str(path)); m = st[0]
    ch = [c for c in m if c.name == name][0]
    return st, ch

def ca_in_ranges(ch, ranges):
    out = []
    for r in ch:
        if any(a <= r.seqid.num <= b for a, b in ranges):
            ca = r.find_atom("CA", "*")
            if ca: out.append(ca.pos)
    return out

def write_subset(st, ch, ranges, path):
    """Write only residues whose seqid falls in `ranges` (one chain), as PDB."""
    ns = gemmi.Structure(); ns.name = st.name; nm = gemmi.Model("1")
    nc = gemmi.Chain(ch.name)
    for r in ch:
        if any(a <= r.seqid.num <= b for a, b in ranges):
            nc.add_residue(r)
    nm.add_chain(nc); ns.add_model(nm); ns.setup_entities()
    path.write_text(ns.make_pdb_string())

def apply_transform(st, tr):
    for m in st:
        for ch in m:
            for r in ch:
                for a in r:
                    a.pos = gemmi.Position(tr.apply(a.pos))

def find_motif_ranges(ch):
    seq = "".join(one(r) for r in ch)
    nums = [r.seqid.num for r in ch]
    ranges = []
    for mot in MOTIFS:
        i = seq.find(mot)
        if i < 0: raise SystemExit(f"motif {mot} not found in chain {ch.name}")
        ranges.append((nums[i], nums[i + len(mot) - 1]))
    return ranges

def main():
    nst, nA = chain(NATIVE, "A")                       # antigen
    natca = ca_in_ranges(nA, EPI_NATIVE)               # 26 epitope CA (native frame)
    write_subset(nst, nA, [(-10**9, 10**9)], OUT/"native_antigen.pdb")   # full antigen
    write_subset(nst, nA, EPI_NATIVE, OUT/"native_epi.pdb")              # epitope only
    _, nB = chain(NATIVE, "B")                          # scFv antibody
    write_subset(nst, nB, [(-10**9, 10**9)], OUT/"scfv.pdb")
    print(f"native: antigen {sum(1 for _ in nA)} res, epitope CA {len(natca)}, scFv {sum(1 for _ in nB)} res")

    for tag, path in DES.items():
        st, ch = chain(path, list(gemmi.read_structure(str(path))[0])[0].name)
        rng = find_motif_ranges(ch)
        dca = ca_in_ranges(ch, rng)
        sup = gemmi.superpose_positions(natca, dca)     # design epitope -> native epitope
        apply_transform(st, sup.transform)
        st.setup_entities()
        (OUT/f"aligned_{tag}.pdb").write_text(st.make_pdb_string())
        write_subset(st, [c for c in st[0]][0], rng, OUT/f"{tag}_epi.pdb")
        print(f"{tag}: epitope resids {rng}  superpose RMSD {sup.rmsd:.2f} A")

if __name__ == "__main__":
    main()
