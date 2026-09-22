#!/usr/bin/env python3
"""Native-context structural features for each epitope, from its RCSB antibody:antigen complex.

Hypothesis (Brandon): how well an epitope can be scaffolded is contingent on its native context, so
add features of the epitope's own structure. John already computed each epitope's helix/strand/loop;
we use those as an answer key to find the antigen chain + epitope inside each complex (the chain whose
antibody-interface secondary structure best matches John's numbers), which both avoids guessing chain
roles and self-validates. Then we read new features off the located epitope:

  epi_size            number of epitope (interface) residues
  epi_rel_sasa        mean relative SASA of the epitope in the FREE antigen (native exposure)
  epi_frac_hydrophobic / _charged / _aromatic / _gly / _pro   composition
  epi_rg              radius of gyration of the epitope Cα (compact vs extended)
  epi_helix/strand/loop_recomputed   (validation vs John)

Run (needs biotite):  <structvenv>/bin/python scripts/dp4_epitope_features.py
"""
from pathlib import Path
import numpy as np, pandas as pd
import biotite.structure as struc
import biotite.structure.io.pdbx as pdbx

ROOT = Path(__file__).resolve().parents[1]
STRUCT = ROOT / "data/dp4_binding/structures"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
OUT = ROOT / "data/dp4_binding/epitope_struct_features.csv"

TIEN = {"ALA":129,"ARG":274,"ASN":195,"ASP":193,"CYS":167,"GLU":223,"GLN":225,"GLY":104,"HIS":224,
        "ILE":197,"LEU":201,"LYS":236,"MET":224,"PHE":240,"PRO":159,"SER":155,"THR":172,"TRP":285,
        "TYR":263,"VAL":174}
HYD = set("AVILMFW"); CHG = set("DEKR"); ARO = set("FWY")
THREE1 = {"ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLU":"E","GLN":"Q","GLY":"G","HIS":"H",
          "ILE":"I","LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W",
          "TYR":"Y","VAL":"V"}

def ss_fracs(sse):
    n = len(sse)
    if n == 0: return (0.0, 0.0, 0.0)
    return (np.mean(sse == "a"), np.mean(sse == "b"), np.mean(sse == "c"))

def analyze(cif, john_hsl):
    f = pdbx.CIFFile.read(str(cif))
    arr = pdbx.get_structure(f, model=1)
    arr = arr[struc.filter_amino_acids(arr) & (arr.element != "H")]
    if arr.array_length() == 0: return None
    # interface residues: any atom within 4 A of an atom in a different chain
    cl = struc.CellList(arr, cell_size=4.0)
    iface = np.zeros(arr.array_length(), dtype=bool)
    for i in range(arr.array_length()):
        nb = cl.get_atoms(arr.coord[i], radius=4.0)
        nb = nb[nb >= 0]
        if np.any(arr.chain_id[nb] != arr.chain_id[i]):
            iface[i] = True

    best = None
    for ch in np.unique(arr.chain_id):
        cmask = arr.chain_id == ch
        chain = arr[cmask]
        try:
            sse = struc.annotate_sse(chain)                 # per-residue, chain order
        except Exception:
            continue
        res_ids = struc.get_residues(chain)[0]
        sse_by_res = dict(zip(res_ids, sse))
        # epitope = this chain's interface residues
        epi_res = np.unique(arr.res_id[cmask & iface])
        if len(epi_res) < 3: continue
        epi_sse = np.array([sse_by_res.get(r, "c") for r in epi_res])
        hsl = ss_fracs(epi_sse)
        d = np.linalg.norm(np.array(hsl) - np.array(john_hsl))
        if best is None or d < best[0]:
            best = (d, ch, epi_res, hsl)
    if best is None: return None
    d, ag_chain, epi_res, hsl = best

    # features on the located epitope, computed on the FREE antigen chain (native exposure)
    ag = arr[arr.chain_id == ag_chain]
    sasa_atom = struc.sasa(ag, vdw_radii="Single")
    # per-residue SASA sum
    res_ids_ag, res_names_ag = struc.get_residues(ag)
    sasa_res = struc.apply_residue_wise(ag, np.nan_to_num(sasa_atom), np.sum)
    sasa_by_res = dict(zip(res_ids_ag, sasa_res))
    rname_by_res = dict(zip(res_ids_ag, res_names_ag))

    rels, comp = [], []
    for r in epi_res:
        nm = rname_by_res.get(r);
        if nm is None: continue
        mx = TIEN.get(nm)
        if mx: rels.append(min(sasa_by_res.get(r, 0.0) / mx, 1.5))
        comp.append(THREE1.get(nm, "X"))
    comp = [c for c in comp if c != "X"]
    ca = ag[(np.isin(ag.res_id, epi_res)) & (ag.atom_name == "CA")]
    rg = float(np.sqrt(np.mean(np.sum((ca.coord - ca.coord.mean(0))**2, axis=1)))) if ca.array_length() > 2 else np.nan

    frac = lambda S: (np.mean([c in S for c in comp]) if comp else np.nan)
    return dict(ag_chain=ag_chain, epi_size=len(epi_res),
                epi_rel_sasa=float(np.mean(rels)) if rels else np.nan, epi_rg=rg,
                epi_frac_hydrophobic=frac(HYD), epi_frac_charged=frac(CHG), epi_frac_aromatic=frac(ARO),
                epi_frac_gly=frac({"G"}), epi_frac_pro=frac({"P"}),
                helix_recomp=hsl[0], strand_recomp=hsl[1], loop_recomp=hsl[2], ss_match_dist=float(d))

def main():
    summ = pd.read_csv(SUMM, encoding="latin-1"); summ.columns = [c.strip() for c in summ.columns]
    rows = []
    for _, r in summ.iterrows():
        epi = str(r["epitope"]).upper(); cif = STRUCT / f"{epi}.cif"
        if not cif.exists(): continue
        try:
            res = analyze(cif, (r["helix"], r["strand"], r["loop"]))
        except Exception as e:
            print(f"  {epi}: FAILED ({type(e).__name__}: {e})"); continue
        if res is None: print(f"  {epi}: no epitope located"); continue
        res["epitope"] = epi.lower(); res["helix_john"] = r["helix"]; rows.append(res)
    df = pd.DataFrame(rows)
    print(f"\nepitopes with features: {len(df)}/{len(summ)}")
    if len(df):
        mae = (df["helix_recomp"] - df["helix_john"]).abs().mean()
        rho = df["helix_recomp"].corr(df["helix_john"])
        print(f"validation vs John's helix:  MAE {mae:.3f}   Pearson r {rho:.3f}   "
              f"(median SS-match dist {df['ss_match_dist'].median():.3f})")
        print("\nnew feature ranges (median [min, max]):")
        for c in ["epi_size","epi_rel_sasa","epi_rg","epi_frac_hydrophobic","epi_frac_charged","epi_frac_aromatic"]:
            print(f"  {c:22s} {df[c].median():7.3f}  [{df[c].min():.3f}, {df[c].max():.3f}]")
        cols = ["epitope","ag_chain","epi_size","epi_rel_sasa","epi_rg","epi_frac_hydrophobic",
                "epi_frac_charged","epi_frac_aromatic","epi_frac_gly","epi_frac_pro",
                "helix_recomp","strand_recomp","loop_recomp","helix_john","ss_match_dist"]
        df[cols].to_csv(OUT, index=False)
        print(f"\nwrote {OUT.relative_to(ROOT)}")

if __name__ == "__main__":
    main()
