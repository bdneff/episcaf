#!/usr/bin/env python3
"""Blind, structure-based epitope propensity for the web tool's predict mode (PROVISIONAL).

Per antigen residue: a contact-number (protrusion) surrogate — surface/protruding residues
(few CA neighbours) score high. This is a classic first-order epitope-likelihood signal and a
stand-in until real DiscoTope-3.0 per-residue scores are wired in (same 0..1 per-residue format,
so swapping is a drop-in). Uses NO antibody info (it is a prediction, not the known contact set).

Writes an antigen-only PDB with the propensity (0..100) in the B-factor column:
  docs/pdb/<id>_scored.pdb
"""
import sys, gemmi, numpy as np
TARGETS = [("8PWW","A"),("7OX3","C"),("5FHX","A"),("8VDL","C")]
R = 10.0  # neighbour radius (Å)
for pid, ag in TARGETS:
    st = gemmi.read_structure(f"docs/pdb/{pid}.pdb"); m = st[0]
    ch = [c for c in m if c.name == ag][0]
    res = [r for r in ch if r.find_atom("CA","*")]
    ca = np.array([[r.find_atom("CA","*").pos.x, r.find_atom("CA","*").pos.y, r.find_atom("CA","*").pos.z] for r in res])
    d = np.sqrt(((ca[:,None,:]-ca[None,:,:])**2).sum(-1))
    cn = (d < R).sum(1) - 1                       # neighbours within R (excl self)
    lo, hi = np.percentile(cn,5), np.percentile(cn,95)
    prop = np.clip((hi - cn)/(hi - lo + 1e-9), 0, 1)   # low CN (protruding) -> high propensity
    # write antigen-only PDB with B = prop*100
    ns = gemmi.Structure(); ns.name = st.name; nm = gemmi.Model("1"); nc = gemmi.Chain(ag)
    for r, p in zip(res, prop):
        for a in r: a.b_iso = float(p*100)
        nc.add_residue(r)
    nm.add_chain(nc); ns.add_model(nm); ns.setup_entities()
    open(f"docs/pdb/{pid}_scored.pdb","w").write(ns.make_pdb_string())
    print(f"  {pid}: {len(res)} res, propensity range {prop.min():.2f}-{prop.max():.2f} -> docs/pdb/{pid}_scored.pdb")
