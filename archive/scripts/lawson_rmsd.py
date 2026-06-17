from __future__ import annotations
import gzip
from pathlib import Path
from typing import Tuple, List, Optional

import numpy as np
import gemmi
import MDAnalysis as mda
from MDAnalysis.analysis import rms

AA3TO1 = {
 "ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G","HIS":"H","ILE":"I","LYS":"K","LEU":"L",
 "MET":"M","ASN":"N","PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V","TRP":"W","TYR":"Y"
}

def read_gemmi_structure(p: Path) -> gemmi.Structure:
    s = str(p)
    if s.endswith(".gz"):
        with gzip.open(p, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
        return gemmi.make_structure_from_block(doc.sole_block())
    if p.suffix.lower() == ".cif":
        doc = gemmi.cif.read(str(p))
        return gemmi.make_structure_from_block(doc.sole_block())
    return gemmi.read_structure(str(p))

def get_chainA(st: gemmi.Structure) -> gemmi.Chain:
    m = st[0]
    for ch in m:
        if ch.name == "A":
            return ch
    return m[0]

def chain_seq_1letter(chain: gemmi.Chain) -> str:
    return "".join(AA3TO1.get(r.name.upper(), "X") for r in chain)

def find_subseq_allowX(hay: str, needle: str) -> int:
    # exact fast path
    i = hay.find(needle)
    if i >= 0:
        return i
    # allow X in hay to match anything
    H = hay
    N = needle
    Lh, Ln = len(H), len(N)
    for i in range(Lh - Ln + 1):
        ok = True
        for a, b in zip(H[i:i+Ln], N):
            if a != "X" and a != b:
                ok = False
                break
        if ok:
            return i
    return -1

def gemmi_bb_positions(chain: gemmi.Chain, start: int, end: int) -> np.ndarray:
    res = [r for r in chain][start:end]
    coords = []
    for r in res:
        for name in ("N","CA","C","O"):
            a = r.find_atom(name, altloc="*")
            if a:
                p = a.pos
                coords.append([p.x,p.y,p.z])
    return np.array(coords, float)

def gemmi_bb_positions_for_resindices(chain: gemmi.Chain, ris: List[int]) -> np.ndarray:
    res = [r for r in chain]
    coords = []
    for i in ris:
        r = res[i]
        for name in ("N","CA","C","O"):
            a = r.find_atom(name, altloc="*")
            if a:
                p = a.pos
                coords.append([p.x,p.y,p.z])
    return np.array(coords, float)

def mda_bb_positions_pdb(pdb: Path) -> np.ndarray:
    u = mda.Universe(str(pdb))
    ag = u.select_atoms("segid A and backbone")
    if len(ag) == 0:
        ag = u.select_atoms("chainid A and backbone")
    return ag.positions.copy()

def mda_bb_positions_pdb_resindices(pdb: Path, ris0: List[int]) -> np.ndarray:
    u = mda.Universe(str(pdb))
    selA = u.select_atoms("segid A")
    if len(selA) == 0:
        selA = u.select_atoms("chainid A")
    pos = selA.residues[ris0].atoms.select_atoms("backbone").positions.copy()
    return pos

def rmsd_superpose(P: np.ndarray, Q: np.ndarray) -> float:
    n = min(len(P), len(Q))
    if n < 3:
        return float("nan")
    return float(rms.rmsd(P[:n], Q[:n], superposition=True))

def overall_rmsd_mpnn_vs_af3_window(mpnn_pdb: Path, af3_cif: Path) -> Tuple[float, int, int]:
    # mpnn seq from mda
    u = mda.Universe(str(mpnn_pdb))
    selA = u.select_atoms("segid A")
    if len(selA) == 0:
        selA = u.select_atoms("chainid A")
    mpnn_seq = selA.residues.sequence(format="string")

    st = read_gemmi_structure(af3_cif)
    chA = get_chainA(st)
    af3_seq = chain_seq_1letter(chA)

    start = find_subseq_allowX(af3_seq, mpnn_seq)
    if start < 0:
        raise ValueError("Could not locate mpnn antigen sequence inside AF3 chain A sequence.")
    end = start + len(mpnn_seq)

    P = mda_bb_positions_pdb(mpnn_pdb)
    Q = gemmi_bb_positions(chA, start, end)
    return rmsd_superpose(P, Q), start, end

def epitope_chunk_rmsd(mpnn_pdb: Path, af3_cif: Path, mpnn_ris0: List[int], af3_ris0: List[int]) -> float:
    st = read_gemmi_structure(af3_cif)
    chA = get_chainA(st)
    P = mda_bb_positions_pdb_resindices(mpnn_pdb, mpnn_ris0)
    Q = gemmi_bb_positions_for_resindices(chA, af3_ris0)
    return rmsd_superpose(P, Q)
