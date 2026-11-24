#!/usr/bin/env python3
"""campaign.evaluate_design on synthetic structures: a design whose AF2 model exactly matches both
its own backbone and the native epitope must score ~0 RMSDs, 0 clashes, and pass. Also the ranked
CSV/zip writer. The subprocess orchestration (rfd4-generate, colabfold) is Gemini-only and not tested
here. Run:  /Users/bneff/micromamba/envs/bcell-repair/bin/python tests/test_rfd4_campaign.py
"""
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from episcaf_service.rfd4 import campaign as C  # noqa: E402
from episcaf_service.rfd4 import read as R  # noqa: E402


def residue(seq_id, comp, ca, chain="A"):
    x, y, z = ca
    return [R.Atom(seq_id, "N", "N", comp, (x - 1, y, z), chain),
            R.Atom(seq_id, "CA", "C", comp, (x, y, z), chain),
            R.Atom(seq_id, "C", "C", comp, (x + 1, y, z), chain),
            R.Atom(seq_id, "O", "O", comp, (x + 1, y + 1, z), chain)]


def _structures():
    cas = {0: (0, 0, 0), 1: (10, 0, 0), 2: (20, 0, 0), 3: (30, 0, 0)}
    comps = {0: "GLY", 1: "ALA", 2: "LEU", 3: "SER"}
    design_atoms = []
    for s, ca in cas.items():
        design_atoms += residue(s, comps[s], ca)
    design = R.RFD4Design(atoms=design_atoms, epitope_seq_ids=(1, 2),
                          native_of={1: ("A", 10), 2: ("A", 11)})
    # native: the two epitope residues at the SAME coords (perfect epitope match), no antibody
    native = R.Structure(atoms=residue(10, "ALA", cas[1]) + residue(11, "LEU", cas[2]))
    # AF2 model: design seq_id + 1 numbering, identical coords (perfect fold)
    af2 = R.Structure(atoms=sum((residue(s + 1, comps[s], cas[s]) for s in cas), []))
    return design, native, af2


EPI_NATIVE = [("A", 10), ("A", 11)]   # the epitope string's ordered native residues


def test_evaluate_perfect_design():
    design, native, af2 = _structures()
    m = C.evaluate_design(design, native, af2, mean_pae=3.0, antibody_chains=(), epitope_native=EPI_NATIVE)
    assert m["epitope_rmsd"] < 1e-6, m["epitope_rmsd"]
    assert m["overall_rmsd"] < 1e-6, m["overall_rmsd"]
    assert m["af_clashes"] == 0, m["af_clashes"]
    assert m["composite"] is not None and m["composite"] > 1.8, m["composite"]   # a clean design passes


def test_backbone_prefilter_no_antibody():
    design, native, _ = _structures()
    assert C.backbone_prefilter_clashes(design, native, antibody_chains=(), epitope_native=EPI_NATIVE) == 0


def test_write_results_ranks_by_composite():
    rows = [
        C.DesignResult("d_lo", "AAAA", 0.9, 1.0, 3.5, 0, 0, 0.7, 0, ""),
        C.DesignResult("d_hi", "LLLL", 0.3, 0.8, 3.0, 0, 0, 2.9, 0, ""),
        C.DesignResult("d_na", "GGGG", float("nan"), 1.0, 3.0, 0, 0, None, 0, ""),
    ]
    with tempfile.TemporaryDirectory() as d:
        csv_path, ranked = C.write_results(Path(d), rows)
        assert ranked[0].design_id == "d_hi"              # highest composite first
        assert ranked[-1].design_id == "d_na"             # None composite sorts last
        text = Path(csv_path).read_text().splitlines()
        assert text[0].startswith("design_id,sequence,epitope_rmsd")
        assert text[1].startswith("d_hi,")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} passed")


if __name__ == "__main__":
    main()
