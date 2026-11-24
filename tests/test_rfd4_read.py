#!/usr/bin/env python3
"""The structure readers: the RFD4 annotated-CIF ingest (epitope vs scaffold via the conditioning
columns + the native-residue mapping) and the standard PDB reader. Stdlib + numpy, no gemmi.

The CIF rows mirror the real RFD4 output format (confirmed on Gemini 2026-11-24): a scaffold residue
has condition_sequence '<M>' and src '.', an epitope residue has its code + a src label like 'A1'
giving the native residue. Run:  /Users/bneff/micromamba/envs/bcell-repair/bin/python tests/test_rfd4_read.py
"""
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from episcaf_service.rfd4 import read as R  # noqa: E402

# minimal CIF in the real column shape (incl. the hyphenated conditioning tag + '<M>'/'.' placeholders)
CIF = """data_x
loop_
_atom_site.group_PDB
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_comp_id
_atom_site.label_seq_id
_atom_site.condition_sequence_1_residue
_atom_site.annotation_src-coordinate-conditioned-region_1_atom
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
ATOM C CA TYR 0 <M> . 3.5 10.6 -7.1
ATOM N N ALA 5 ALA A1 7.0 6.0 -0.5
ATOM C CA ALA 5 ALA A1 8.4 7.0 -0.7
ATOM C C ALA 5 ALA A1 9.0 6.0 0.0
ATOM O O ALA 5 ALA A1 9.5 5.0 0.2
ATOM C CA LEU 6 LEU A2 5.5 5.8 1.4
ATOM C CA ILE 13 <M> . 8.8 -2.9 7.4
"""


def _cif_file():
    p = tempfile.NamedTemporaryFile("w", suffix=".cif", delete=False)
    p.write(CIF)
    p.close()
    return p.name


def test_rfd4_epitope_vs_scaffold_and_native_map():
    d = R.read_rfd4_cif(_cif_file())
    assert d.seq_ids == [0, 5, 6, 13], d.seq_ids
    assert d.sequence == "YALI", d.sequence               # TYR,ALA,LEU,ILE
    assert d.epitope_seq_ids == (5, 6), d.epitope_seq_ids
    assert d.scaffold_seq_ids == (0, 13), d.scaffold_seq_ids
    assert d.native_of == {5: ("A", 1), 6: ("A", 2)}, d.native_of


def test_rfd4_coord_extraction():
    d = R.read_rfd4_cif(_cif_file())
    assert d.ca_xyz().shape == (4, 3)                      # one CA per residue
    assert d.ca_xyz([5])[0].tolist() == [8.4, 7.0, -0.7]
    # backbone of residue 5 ordered N, CA, C, O
    bb = d.backbone_xyz([5])
    assert bb.shape == (4, 3)
    assert bb[0].tolist() == [7.0, 6.0, -0.5] and bb[1].tolist() == [8.4, 7.0, -0.7]
    xyz, resid = d.heavy_xyz_with_resid([5])
    assert xyz.shape == (4, 3) and set(resid.tolist()) == {5}   # N,CA,C,O all heavy


PDB = """\
ATOM      1  N   ALA A   1      11.000  12.000  13.000  1.00  0.00           N
ATOM      2  CA  ALA A   1      12.000  12.000  13.000  1.00  0.00           C
ATOM      3  C   ALA A   1      13.000  12.000  13.000  1.00  0.00           C
ATOM      4  O   ALA A   1      14.000  12.000  13.000  1.00  0.00           O
ATOM      5  CA  GLY B   7      20.000  20.000  20.000  1.00  0.00           C
HETATM    6  O   HOH B   8      30.000  30.000  30.000  1.00  0.00           O
"""


def test_read_pdb():
    p = tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False)
    p.write(PDB)
    p.close()
    s = R.read_pdb(p.name)
    assert len(s.atoms) == 6
    assert s.ca_xyz(["A"]).tolist() == [[12.0, 12.0, 13.0]]
    assert s.backbone_xyz("A", [1]).shape == (4, 3)        # N,CA,C,O for chain A res 1
    # chain B heavy atoms = the GLY CA + the HOH O
    assert s.heavy_xyz(["B"]).shape == (2, 3)


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} passed")


if __name__ == "__main__":
    main()
