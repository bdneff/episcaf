"""Tier-1 interpretable native-context features for the binding model (per epitope).

Computed from the epitope's native antibody:antigen complex (``data/dp4_binding/structures/<ID>.cif``),
reusing the epitope-location rule validated in ``scripts/dp4_epitope_features.py``: the epitope is the
set of antigen-chain residues with any heavy atom within 4 A of another chain (the antibody). The
antigen chain is taken from the self-validated ``epitope_struct_features.csv`` (``ag_chain``), and we
reproduce its ``epi_size`` exactly as a check (see ``extract_all``).

These are all static, apo/native-context quantities (causal before an antibody exists), the direction
the DP4 result said carried the signal. Each is a single interpretable number per epitope:

  epi_depth        mean residue depth of the epitope (A): distance from epitope atoms to the nearest
                   solvent-exposed ("surface") atom. SASA saturates once a residue is exposed; depth
                   keeps ranking how buried / pocketed the patch is.
  epi_cn           mean Ca contact-number of the epitope (Thornton protrusion surrogate, reused from
                   scripts/webtool_epitope_propensity.py): Ca neighbours within 10 A. LOW => the patch
                   sticks out => more antibody-accessible.
  ag_out25         native antigen mass the paratope must clear: antigen atoms on the antibody-approach
  ag_out25_hi      side within 25 A of the epitope frame origin (ag_out25), the subset >5 A out along
  ag_near15        the approach normal (ag_out25_hi), and antigen atoms within 15 A (ag_near15).
                   Ported verbatim from episcaf_analysis.aomclash.features.target_features; dropping
                   this block was AOMclash's single biggest ablation hit.
  frame_ratio      epitope-plane anisotropy from the same frame (free by-product; elongated vs round).
  epi_glyco_dist   min Ca distance from the epitope to the nearest N-linked sequon (N-X-S/T, X!=P) Asn
                   on the antigen chain (A); SMALL => glycan-shielded => poorly accessible. 30 A cap
                   when no sequon exists.
  epi_net_charge   sum of formal charges over the epitope (K,R +1; D,E -1; H +0.5).
  epi_abs_charge   sum of |formal charge| (total charged content, sign-blind).
  epi_coulomb      cheap intra-epitope Coulomb sum over Ca: sum_{i<j} q_i q_j / r_ij (kept
                   dimensionless; a static electrostatic-organisation scalar, not a PB energy).

Run (regenerates results/feature_model/native_tier1_features.csv):
  <env>/bin/python -m episcaf_analysis.binding_model.native_features
"""
from __future__ import annotations

import re
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[2]
STRUCT = ROOT / "data/dp4_binding/structures"
STRUCTF = ROOT / "data/dp4_binding/epitope_struct_features.csv"
SUMM = ROOT / "data/dp4_binding/john/scaffoldedEpitopeSummary.csv"
OUT = ROOT / "results/feature_model/native_tier1_features.csv"

# element VdW radii (A) for the SASA calc; default 1.7 for anything unlisted
VDW = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80, "H": 1.20}
THREE1 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLU": "E", "GLN": "Q",
          "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
          "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
FORMAL = {"K": 1.0, "R": 1.0, "D": -1.0, "E": -1.0, "H": 0.5}
SEQUON = re.compile(r"N[^P][ST]")

NEW_FEATURES = ["epi_depth", "epi_cn", "ag_out25", "ag_out25_hi", "ag_near15",
                "frame_ratio", "epi_glyco_dist", "epi_net_charge", "epi_abs_charge", "epi_coulomb"]


# ---------------------------------------------------------------------------
# structure reading + epitope location (Bio.PDB; no biotite/DSSP needed)
# ---------------------------------------------------------------------------
def _load_chains(cif_path):
    """First model as {chain_id: [(resid, resname, {atom_name: coord}, element_list)]}, ATOM only."""
    from Bio.PDB import MMCIFParser
    model = next(iter(MMCIFParser(QUIET=True).get_structure("x", str(cif_path))))
    chains = {}
    for ch in model:
        residues = []
        for r in ch:
            if r.id[0] != " ":          # skip hetero / water
                continue
            atoms = {}
            for a in r:
                if a.element == "H":
                    continue
                atoms[a.get_name()] = np.asarray(a.coord, float)
            if atoms:
                residues.append((r.id[1], r.resname, atoms))
        if residues:
            chains[ch.id] = residues
    return chains


def locate_epitope(chains, ag_chain, cutoff=4.0):
    """Epitope = antigen-chain residues with any heavy atom within `cutoff` of another chain.

    Returns (epi_resids, epitope_resname_by_resid) for the antigen chain. Matches the interface rule
    in scripts/dp4_epitope_features.py (there the "other chain" is the antibody).
    """
    from scipy.spatial import cKDTree
    other = [xyz for cid, res in chains.items() if cid != ag_chain
             for _, _, atoms in res for xyz in atoms.values()]
    other = np.asarray(other, float)
    tree = cKDTree(other)
    epi, names = [], {}
    for resid, resname, atoms in chains[ag_chain]:
        coords = np.asarray(list(atoms.values()), float)
        if tree.query(coords)[0].min() < cutoff:
            epi.append(resid)
            names[resid] = resname
    return epi, names


# ---------------------------------------------------------------------------
# individual feature extractors (each takes already-parsed arrays; unit-testable)
# ---------------------------------------------------------------------------
def residue_depth(ag_atom_xyz, ag_atom_elem, epi_atom_xyz, surf_cut=1.0):
    """Mean epitope residue depth (A): distance from each epitope heavy atom to the nearest
    solvent-exposed antigen atom. Surface atoms are those with per-atom SASA > surf_cut (A^2),
    computed on the free antigen chain with freesasa.calcCoord.
    """
    import freesasa
    from scipy.spatial import cKDTree
    radii = [VDW.get(e, 1.70) for e in ag_atom_elem]
    res = freesasa.calcCoord(ag_atom_xyz.ravel().tolist(), radii)
    area = np.array([res.atomArea(i) for i in range(res.nAtoms())])
    surf = ag_atom_xyz[area > surf_cut]
    if len(surf) == 0:
        return float("nan")
    d, _ = cKDTree(surf).query(epi_atom_xyz)
    return float(d.mean())


def contact_number(ag_ca, epi_ca_mask, radius=10.0):
    """Mean Ca contact-number of the epitope (protrusion surrogate, reused from
    webtool_epitope_propensity.py): neighbours within `radius` A, excluding self. Low = protruding.
    """
    D = np.sqrt(((ag_ca[:, None, :] - ag_ca[None, :, :]) ** 2).sum(-1))
    cn = (D < radius).sum(1) - 1
    return float(cn[epi_ca_mask].mean())


def glyco_distance(chains, ag_chain, epi_resids, cap=30.0):
    """Min Ca distance (A) from the epitope to the nearest N-X-S/T sequon Asn on the antigen chain.
    Sequence is read in residue order from the antigen chain; `cap` returned when no sequon exists.
    """
    res = chains[ag_chain]
    seq = "".join(THREE1.get(rn, "X") for _, rn, _ in res)
    ca = {resid: atoms["CA"] for resid, _, atoms in res if "CA" in atoms}
    sequon_ca = [ca[res[m.start()][0]] for m in SEQUON.finditer(seq)
                 if res[m.start()][0] in ca]
    epi_ca = np.array([ca[r] for r in epi_resids if r in ca])
    if not sequon_ca or len(epi_ca) == 0:
        return float(cap)
    from scipy.spatial import cKDTree
    d, _ = cKDTree(np.array(sequon_ca)).query(epi_ca)
    return float(min(d.min(), cap))


def charge_features(epi_resnames, epi_ca):
    """(net_charge, abs_charge, coulomb) from formal charges over the epitope. Coulomb is the cheap
    dimensionless pairwise sum sum_{i<j} q_i q_j / r_ij over charged-residue Ca.
    """
    q1 = {THREE1.get(rn, "X"): FORMAL.get(THREE1.get(rn, "X"), 0.0) for rn in set(epi_resnames)}
    q = np.array([FORMAL.get(THREE1.get(rn, "X"), 0.0) for rn in epi_resnames])
    net = float(q.sum())
    absq = float(np.abs(q).sum())
    chg = np.where(q != 0.0)[0]
    coulomb = 0.0
    for a in range(len(chg)):
        for b in range(a + 1, len(chg)):
            i, j = chg[a], chg[b]
            r = np.linalg.norm(epi_ca[i] - epi_ca[j]) + 1e-6
            coulomb += q[i] * q[j] / r
    return net, absq, float(coulomb)


# ---------------------------------------------------------------------------
# per-epitope driver
# ---------------------------------------------------------------------------
def extract_one(cif_path, ag_chain, n_isl=1):
    """All Tier-1 native features for one epitope. Returns (feature_dict, epi_size)."""
    from episcaf_analysis.aomclash.features import target_features

    chains = _load_chains(cif_path)
    if ag_chain not in chains:
        raise ValueError(f"ag_chain {ag_chain} not in {sorted(chains)}")
    epi_resids, names = locate_epitope(chains, ag_chain)
    if len(epi_resids) < 3:
        raise ValueError("fewer than 3 epitope residues located")

    ag_res = chains[ag_chain]
    ag_atom_xyz = np.array([xyz for _, _, atoms in ag_res for xyz in atoms.values()])
    ag_atom_elem = [nm[0] for _, _, atoms in ag_res for nm in atoms.keys()]
    ag_ca = np.array([atoms["CA"] for _, _, atoms in ag_res if "CA" in atoms])
    ca_resids = [resid for resid, _, atoms in ag_res if "CA" in atoms]
    epi_set = set(epi_resids)
    epi_ca_mask = np.array([r in epi_set for r in ca_resids])
    epi_ca = ag_ca[epi_ca_mask]
    epi_atom_xyz = np.array([xyz for resid, _, atoms in ag_res if resid in epi_set
                             for xyz in atoms.values()])
    epi_resnames = [rn for resid, rn, _ in ag_res if resid in epi_set]

    tf = target_features(epi_ca, ag_ca, ag_atom_xyz, n_isl)
    if tf is None:
        crowd = dict(ag_out25=np.nan, ag_out25_hi=np.nan, ag_near15=np.nan, frame_ratio=np.nan)
    else:
        _, _, _, _, t = tf
        crowd = dict(ag_out25=t["ag_out25"], ag_out25_hi=t["ag_out25_hi"],
                     ag_near15=t["ag_near15"], frame_ratio=t["frame_ratio"])

    net, absq, coulomb = charge_features(epi_resnames, epi_ca)
    feat = dict(
        epi_depth=residue_depth(ag_atom_xyz, ag_atom_elem, epi_atom_xyz),
        epi_cn=contact_number(ag_ca, epi_ca_mask),
        epi_glyco_dist=glyco_distance(chains, ag_chain, epi_resids),
        epi_net_charge=net, epi_abs_charge=absq, epi_coulomb=coulomb,
        **crowd,
    )
    return feat, len(epi_resids)


def extract_all(struct_dir=STRUCT, structf=STRUCTF, summ=SUMM, out=OUT):
    """Compute Tier-1 native features for every epitope with a local CIF; validate epi_size."""
    st = pd.read_csv(structf)
    isl = pd.read_csv(summ, encoding="latin-1")
    isl.columns = [c.strip() for c in isl.columns]
    isl = isl.rename(columns={"#islands": "n_islands"})
    isl["epitope"] = isl["epitope"].astype(str).str.lower()
    nmap = dict(zip(isl["epitope"], isl["n_islands"].fillna(1).astype(int)))

    rows, mismatches = [], []
    for epi, agc, size in st[["epitope", "ag_chain", "epi_size"]].values:
        cif = struct_dir / f"{str(epi).upper()}.cif"
        if not cif.exists():
            continue
        try:
            feat, my_size = extract_one(cif, agc, n_isl=nmap.get(str(epi).lower(), 1))
        except Exception as e:  # noqa: BLE001
            print(f"  {epi}: FAILED ({type(e).__name__}: {e})")
            continue
        if my_size != int(size):
            mismatches.append((epi, int(size), my_size))
        feat["epitope"] = str(epi).lower()
        feat["epi_size_check"] = my_size
        rows.append(feat)

    df = pd.DataFrame(rows)[["epitope"] + NEW_FEATURES + ["epi_size_check"]]
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"epitopes with Tier-1 features: {len(df)}/{len(st)}")
    print(f"epi_size mismatches vs epitope_struct_features.csv: {len(mismatches)}", mismatches or "")
    print("\nfeature ranges (median [min, max]):")
    for c in NEW_FEATURES:
        print(f"  {c:16s} {df[c].median():9.3f}  [{df[c].min():.3f}, {df[c].max():.3f}]")
    print(f"\nwrote {out.relative_to(ROOT)}")
    return df


if __name__ == "__main__":
    extract_all()
