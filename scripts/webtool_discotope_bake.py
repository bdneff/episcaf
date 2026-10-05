#!/usr/bin/env python3
"""Bake real DiscoTope-3.0 (or BepiPred-3.0) per-residue epitope scores into the web-tool viewer.

The predict-mode viewer colors the antigen by the B-factor of docs/pdb/<id>_scored.pdb. This script
replaces the provisional protrusion surrogate (scripts/webtool_epitope_propensity.py) with real
predictor scores: it reads a per-residue CSV, maps each score onto the matching antigen residue, and
rewrites docs/pdb/<id>_scored.pdb with score*100 in the B-factor column. Nothing else in the viewer
changes, so this is a drop-in once you have the CSVs.

HOW TO GET THE CSVs (the one step that can't run on a laptop here):
  DiscoTope-3.0 (structure-based; recommended) — https://services.healthtech.dtu.dk/services/DiscoTope-3.0/
  or the Colab linked there. Upload each ANTIGEN-ONLY PDB (already served, antigen chain only):
      docs/pdb/8PWW_scored.pdb  docs/pdb/7OX3_scored.pdb  docs/pdb/5FHX_scored.pdb  docs/pdb/8VDL_scored.pdb
  Download each per-residue result CSV, then:
      python scripts/webtool_discotope_bake.py 8PWW=path/to/8pww.csv 7OX3=... 5FHX=... 8VDL=...

The CSV needs a residue-number column and a score column; this script auto-detects common DiscoTope-3.0
/ BepiPred-3.0 headers and prints what it used. Override with --resid-col / --score-col if needed.
"""
import sys, csv, argparse, gemmi

ANTIGEN = {"8PWW": "A", "7OX3": "C", "5FHX": "A", "8VDL": "C"}
SCORE_HINTS = ["calibrated", "discotope", "bepipred", "epitope_score", "score", "prob"]
RESID_HINTS = ["res_id", "resid", "residue_number", "position", "resi", "pos", "num"]

def pick(cols, hints):
    low = {c.lower(): c for c in cols}
    for h in hints:
        for lc, orig in low.items():
            if h in lc:
                return orig
    return None

def load_scores(path, resid_col, score_col):
    with open(path, newline="") as f:
        # sniff delimiter (DiscoTope uses comma; some exports use tab/space)
        sample = f.read(4096); f.seek(0)
        delim = "\t" if sample.count("\t") > sample.count(",") else ","
        rows = list(csv.DictReader(f, delimiter=delim))
    cols = rows[0].keys()
    rc = resid_col or pick(cols, RESID_HINTS)
    sc = score_col or pick(cols, SCORE_HINTS)
    if not rc or not sc:
        sys.exit(f"  could not find columns in {path}: resid={rc} score={sc}; columns={list(cols)}")
    out = {}
    for r in rows:
        try:
            out[int(float(r[rc]))] = float(r[sc])
        except (ValueError, TypeError):
            continue
    return out, rc, sc

def bake(pid, csv_path, resid_col, score_col):
    ag = ANTIGEN[pid]
    scores, rc, sc = load_scores(csv_path, resid_col, score_col)
    vals = list(scores.values()); lo, hi = min(vals), max(vals)
    rng = (hi - lo) or 1.0
    st = gemmi.read_structure(f"docs/pdb/{pid}.pdb"); m = st[0]
    ch = [c for c in m if c.name == ag][0]
    res = [r for r in ch if r.find_atom("CA", "*")]
    ns = gemmi.Structure(); ns.name = st.name; nm = gemmi.Model("1"); nc = gemmi.Chain(ag)
    matched = 0
    for r in res:
        s = scores.get(r.seqid.num)
        b = ((s - lo) / rng * 100.0) if s is not None else 0.0
        if s is not None: matched += 1
        for a in r: a.b_iso = float(b)
        nc.add_residue(r)
    nm.add_chain(nc); ns.add_model(nm); ns.setup_entities()
    open(f"docs/pdb/{pid}_scored.pdb", "w").write(ns.make_pdb_string())
    print(f"  {pid}: cols[res={rc}, score={sc}] matched {matched}/{len(res)} residues, "
          f"score {lo:.3f}-{hi:.3f} -> docs/pdb/{pid}_scored.pdb")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("pairs", nargs="+", help="ID=csv_path (e.g. 8PWW=8pww_discotope.csv)")
    ap.add_argument("--resid-col"); ap.add_argument("--score-col")
    a = ap.parse_args()
    for p in a.pairs:
        pid, path = p.split("=", 1)
        if pid not in ANTIGEN: sys.exit(f"unknown target {pid}; expected one of {list(ANTIGEN)}")
        bake(pid, path, a.resid_col, a.score_col)
    print("done. rebuild not needed; the viewer reads docs/pdb/<id>_scored.pdb directly.")
