#!/usr/bin/env python3
"""Emit a growPatch viewer data file (residues + relSASA + candidate patches) for docs/growpatch.html.

Reads an antigen CIF, runs the growPatch enumerator (neighbor-based growth), and writes a small JS
file setting `window.GROWPATCH_DATA = {...}` that the viewer page loads. Regenerable:

  /Users/bneff/micromamba/envs/bcell-repair/bin/python scripts/growpatch_viz_data.py \
      --cif data/dp4_binding/structures/8PWW.cif --chain A \
      --scored-pdb docs/pdb/8PWW_scored.pdb \
      --out docs/data/growpatch_8pww.js --pdb-id 8pww

The residue `p` is the side-chain-centroid representative point (what the enumerator grows on), used
for the flat-surface projection; `resi` is the PDB residue number, so the 3D viewer can select it.
Each candidate carries `resis` (the full gap-filled residue set sent to scaffold design) and
`resis_rank` (the growth-selected surface residues used for DiscoTope ranking).

DiscoTope ranking (optional, --scored-pdb): if a DiscoTope-3.0-scored PDB is given, each candidate
gets a within-antigen epitope-likelihood summary. The score is read from the CA B-factor of the
scored PDB (baked by scripts/webtool_discotope_bake.py as a per-residue DiscoTope-3.0 likelihood,
min-max normalized across the chain and written on a 0-100 scale). We divide by --disco-scale
(default 100) back to a 0-1 normalized score and attach, per candidate:
    disco_mean   mean normalized DiscoTope-3.0 likelihood over the candidate's scored ranking residues
    disco_n      how many of the ranking residues had a score
    disco_total  how many ranking residues the candidate has (keys_raw)
The mean is taken over the growth-selected surface residues (resis_rank), excluding gap-fill residues.
Residues with no score are skipped, never treated as 0; a candidate with no scored ranking residue
gets disco_mean null ("No score"). This is a PREDICTION of epitope likelihood, distinct from the
geometric growPatch sampling, and is not antibody-binding truth. Candidates are written in descending
disco_mean order (No-score candidates last). Targets with no scored PDB omit the disco fields and the
viewer shows the full candidate set unranked.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import warnings

warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from episcaf_analysis.growpatch import read as R   # noqa: E402


def _resi(label):
    m = re.match(r"-?\d+", str(label))
    return int(m.group()) if m else None


def load_disco_scores(scored_pdb, chain, scale):
    """Map PDB residue number -> normalized DiscoTope-3.0 likelihood (0-1) from CA B-factors.

    The scored PDB is produced by scripts/webtool_discotope_bake.py, which writes the per-residue
    DiscoTope-3.0 likelihood (min-max normalized across the chain) into the B-factor column on a
    0-100 scale. Dividing by `scale` (100) recovers the 0-1 normalized score. We read the CA atom's
    B-factor as the residue score.
    """
    from Bio.PDB import PDBParser

    structure = PDBParser(QUIET=True).get_structure("scored", scored_pdb)
    model = next(iter(structure))
    if chain not in model:
        sys.exit(f"  scored PDB {scored_pdb} has no chain {chain}; chains={[c.id for c in model]}")
    scores = {}
    for res in model[chain]:
        if not res.has_id("CA"):
            continue
        scores[res.id[1]] = res["CA"].get_bfactor() / scale
    if not scores:
        sys.exit(f"  scored PDB {scored_pdb} chain {chain} has no CA atoms to read scores from")
    return scores


def attach_disco(patches, scores):
    """Attach disco_mean / disco_n / disco_total to each candidate, in place.

    The mean is taken over the candidate's growth-selected ranking residues (resis_rank), which are
    the exposed surface residues growPatch selected, excluding gap-fill residues. Residues with no
    score are skipped, never treated as 0. A candidate with no scored ranking residue gets
    disco_mean = None so the viewer can show "No score". disco_total is the number of ranking residues
    so the viewer can show "N of M residues scored".
    """
    total_missing = 0
    for p in patches:
        rank = p["resis_rank"]
        vals = [scores[r] for r in rank if r in scores]
        total_missing += len(rank) - len(vals)
        p["disco_mean"] = round(sum(vals) / len(vals), 4) if vals else None
        p["disco_n"] = len(vals)
        p["disco_total"] = len(rank)
    return total_missing


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cif", required=True)
    ap.add_argument("--chain", default="A")
    ap.add_argument("--pdb-id", default="")
    ap.add_argument("--sasa-cut", type=float, default=0.2)
    ap.add_argument("--scored-pdb", default="",
                    help="DiscoTope-3.0-scored PDB (CA B-factor = per-residue likelihood); enables ranking")
    ap.add_argument("--scored-chain", default="",
                    help="chain in the scored PDB (defaults to --chain)")
    ap.add_argument("--disco-scale", type=float, default=100.0,
                    help="divisor to map the scored B-factor back to a 0-1 score (bake writes 0-100)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    ag = R.read_antigen(a.cif, a.chain, relsasa=True, verbose=True)
    cands, _ = R.enumerate_candidates(a.cif, a.chain, antigen=ag, sasa_cut=a.sasa_cut)

    residues = []
    for r in ag.residues:
        c, slot = r["key"]
        residues.append(dict(resi=_resi(ag.label(c, slot)), resn=r["resname"],
                             rel=round(float(ag.rs(r["key"])), 3),
                             surf=bool(ag.rs(r["key"]) >= a.sasa_cut),
                             p=[round(float(x), 2) for x in r["rep"]]))

    patches = []
    for p in cands:
        resis = sorted({_resi(ag.label(c, s)) for (c, s) in p["keys"]})
        resis_rank = sorted({_resi(ag.label(c, s)) for (c, s) in p["keys_raw"]})
        patches.append(dict(i=p["patch"], epitope=p["epitope"],
                            n_res=p["n_res"], n_isl=p["n_islands"],
                            resis=resis, resis_rank=resis_rank))

    has_disco = bool(a.scored_pdb)
    if has_disco:
        scores = load_disco_scores(a.scored_pdb, a.scored_chain or a.chain, a.disco_scale)
        total_missing = attach_disco(patches, scores)
        # rank within the antigen by descending mean DiscoTope score; No-score candidates last;
        # stable id `i` as the deterministic tie-breaker. Candidate `i` stays its enumeration id.
        patches.sort(key=lambda p: (p["disco_mean"] is None, -(p["disco_mean"] or 0.0), p["i"]))
        n_scored = len(scores)
        print(f"  DiscoTope: {n_scored} residues scored from {a.scored_pdb}; "
              f"{total_missing} ranking-residue lookups had no score (skipped, not zeroed)")

    pdb_id = a.pdb_id or os.path.basename(a.cif).split(".")[0]
    data = dict(pdb=pdb_id, chain=a.chain, sasa_cut=a.sasa_cut,
                n_residues=len(residues), n_surface=sum(r["surf"] for r in residues),
                n_patches=len(patches), has_disco=has_disco,
                residues=residues, patches=patches)

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("// auto-generated by scripts/growpatch_viz_data.py -- do not edit by hand\n")
        fh.write("window.GROWPATCH_DATA = ")
        json.dump(data, fh, separators=(",", ":"))
        fh.write(";\n")
    print(f"wrote {a.out}: {len(residues)} residues ({data['n_surface']} surface), "
          f"{len(patches)} growPatch candidates (has_disco={has_disco})")


if __name__ == "__main__":
    main()
