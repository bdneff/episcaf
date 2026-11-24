"""Structure readers (stdlib only — no gemmi/Bio dependency, runs identically here and on Gemini).

Two readers:
  read_rfd4_cif(path)   the RFD4-Proteina output CIF. Beyond coordinates it carries per-atom
                        annotation columns that mark the conditioning, which is how we recover
                        where the epitope landed (replacing John's motif_residue_map):
                          condition_sequence_1_residue                     "<M>" = designed scaffold;
                                                                           a residue code = fixed epitope
                          annotation_src-coordinate-conditioned-region..   "." scaffold; "A1","A2",... =
                                                                           native source residue of an
                                                                           epitope position (the mapping)
  read_pdb(path)        a standard PDB (the AF2 model, the native antigen/antibody) → atoms by chain.

Both return atom arrays the clash/score code consumes. Tuned to these two formats, not a general
mmCIF/PDB parser.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

_3TO1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E",
    "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F",
    "PRO": "P", "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
_BACKBONE = ("N", "CA", "C", "O")
_MASK_TOKENS = {"<M>", ".", "?", ""}
_SRC_RE = re.compile(r"^([A-Za-z]+)(-?\d+)$")


@dataclass
class Atom:
    seq_id: int           # residue index within the chain
    name: str             # atom name (CA, N, ...)
    element: str          # type symbol (C, N, O, ...)
    comp: str             # residue 3-letter code
    xyz: tuple            # (x, y, z)
    chain: str = "A"


@dataclass
class RFD4Design:
    atoms: list = field(default_factory=list)
    epitope_seq_ids: tuple = ()              # design residue indices that are the fixed epitope
    native_of: dict = field(default_factory=dict)   # design seq_id -> (native_chain, native_resnum)

    @property
    def seq_ids(self) -> list:
        seen, out = set(), []
        for a in self.atoms:
            if a.seq_id not in seen:
                seen.add(a.seq_id)
                out.append(a.seq_id)
        return out

    @property
    def sequence(self) -> str:
        comp = {}
        for a in self.atoms:
            comp.setdefault(a.seq_id, a.comp)
        return "".join(_3TO1.get(comp[s], "X") for s in self.seq_ids)

    @property
    def scaffold_seq_ids(self) -> tuple:
        epi = set(self.epitope_seq_ids)
        return tuple(s for s in self.seq_ids if s not in epi)

    def _xyz(self, atoms) -> np.ndarray:
        return np.array([a.xyz for a in atoms], float) if atoms else np.empty((0, 3))

    def ca_xyz(self, seq_ids=None) -> np.ndarray:
        want = set(self.seq_ids if seq_ids is None else seq_ids)
        return self._xyz([a for a in self.atoms if a.name == "CA" and a.seq_id in want])

    def backbone_xyz(self, seq_ids) -> np.ndarray:
        want = set(seq_ids)
        # ordered by (seq_id, N/CA/C/O) so two structures' epitope backbones pair up atom-for-atom
        order = {n: i for i, n in enumerate(_BACKBONE)}
        atoms = [a for a in self.atoms if a.seq_id in want and a.name in order]
        atoms.sort(key=lambda a: (a.seq_id, order[a.name]))
        return self._xyz(atoms)

    def heavy_xyz_with_resid(self, seq_ids):
        """Heavy (non-H) atom coords + their residue ids, for the scaffold clash set."""
        want = set(seq_ids)
        atoms = [a for a in self.atoms if a.seq_id in want and a.element.upper() != "H"]
        return self._xyz(atoms), np.array([a.seq_id for a in atoms])


def read_rfd4_cif(path) -> RFD4Design:
    columns, rows = [], []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("_atom_site."):
                columns.append(s.split(".", 1)[1].split()[0])
            elif s.startswith("ATOM ") or s.startswith("HETATM "):
                rows.append(s.split())
    if not columns or not rows:
        raise ValueError(f"no _atom_site loop/atoms found in {path}")
    ix = {c: i for i, c in enumerate(columns)}
    need = ("label_atom_id", "type_symbol", "label_comp_id", "label_seq_id",
            "Cartn_x", "Cartn_y", "Cartn_z", "condition_sequence_1_residue",
            "annotation_src-coordinate-conditioned-region_1_atom")
    for c in need:
        if c not in ix:
            raise ValueError(f"RFD4 CIF missing column {c!r}")

    atoms, epitope, native = [], set(), {}
    for r in rows:
        seq = int(r[ix["label_seq_id"]])
        atoms.append(Atom(seq_id=seq, name=r[ix["label_atom_id"]], element=r[ix["type_symbol"]],
                          comp=r[ix["label_comp_id"]],
                          xyz=(float(r[ix["Cartn_x"]]), float(r[ix["Cartn_y"]]), float(r[ix["Cartn_z"]]))))
        if r[ix["condition_sequence_1_residue"]] not in _MASK_TOKENS:
            epitope.add(seq)
        src = r[ix["annotation_src-coordinate-conditioned-region_1_atom"]]
        if seq not in native:
            m = _SRC_RE.match(src)
            if m:
                native[seq] = (m.group(1), int(m.group(2)))
                epitope.add(seq)
    return RFD4Design(atoms=atoms, epitope_seq_ids=tuple(sorted(epitope)), native_of=native)


# --- standard PDB (AF2 model, native antigen/antibody) ---

@dataclass
class Structure:
    atoms: list = field(default_factory=list)   # Atom with .chain set

    def chain_atoms(self, chains):
        chains = set(chains)
        return [a for a in self.atoms if a.chain in chains]

    def heavy_xyz(self, chains=None):
        a = self.atoms if chains is None else self.chain_atoms(chains)
        a = [x for x in a if x.element.upper() != "H"]
        return np.array([x.xyz for x in a], float) if a else np.empty((0, 3))

    def backbone_xyz(self, chain, resnums):
        want = set(resnums)
        order = {n: i for i, n in enumerate(_BACKBONE)}
        a = [x for x in self.atoms if x.chain == chain and x.seq_id in want and x.name in order]
        a.sort(key=lambda x: (x.seq_id, order[x.name]))
        return np.array([x.xyz for x in a], float) if a else np.empty((0, 3))

    def ca_xyz(self, chains=None):
        a = self.atoms if chains is None else self.chain_atoms(chains)
        a = [x for x in a if x.name == "CA"]
        return np.array([x.xyz for x in a], float) if a else np.empty((0, 3))


def read_pdb(path) -> Structure:
    atoms = []
    with open(path) as fh:
        for line in fh:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            name = line[12:16].strip()
            comp = line[17:20].strip()
            chain = line[21].strip() or "A"
            resseq = int(line[22:26])
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            element = line[76:78].strip() or name[:1]
            atoms.append(Atom(seq_id=resseq, name=name, element=element, comp=comp,
                              xyz=(x, y, z), chain=chain))
    return Structure(atoms=atoms)
