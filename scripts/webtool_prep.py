#!/usr/bin/env python3
"""Prep one antibody-antigen complex for the Epitope Scaffolder web tool render.

Keeps only the antigen + antibody chains, computes the contact epitope (antigen residues with any
heavy atom within CUTOFF of an antibody heavy atom), writes a trimmed PDB, and prints the epitope
resid list (space-separated) for the VMD selection.

Usage: webtool_prep.py <in.pdb|.cif> <antigen_chain> <antibody_chains_csv> <out.pdb> [cutoff]
"""
import sys, gemmi, numpy as np

pin, ag, abch, out = sys.argv[1], sys.argv[2], sys.argv[3].split(","), sys.argv[4]
cut = float(sys.argv[5]) if len(sys.argv) > 5 else 4.5

st = gemmi.read_structure(pin); st.setup_entities(); m = st[0]
abset, keep = set(abch), set([ag]) | set(abch)

def heavy(res):
    return [[a.pos.x, a.pos.y, a.pos.z] for a in res if a.element.name != "H"]

ab = []
for ch in m:
    if ch.name in abset:
        for r in ch:
            if r.name != "HOH": ab += heavy(r)
ab = np.array(ab)

epi = []
agch = [c for c in m if c.name == ag][0]
for r in agch:
    if r.name == "HOH": continue
    pts = heavy(r)
    if not pts: continue
    pts = np.array(pts)
    d = np.sqrt(((pts[:, None, :] - ab[None, :, :]) ** 2).sum(-1))
    if d.min() <= cut: epi.append(r.seqid.num)

ns = gemmi.Structure(); ns.name = st.name; nm = gemmi.Model("1")
for ch in m:
    if ch.name in keep:
        nc = gemmi.Chain(ch.name)
        for r in ch:
            if r.name != "HOH": nc.add_residue(r)
        nm.add_chain(nc)
ns.add_model(nm); ns.setup_entities()
open(out, "w").write(ns.make_pdb_string())
print(" ".join(str(x) for x in epi))
