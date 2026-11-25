"""Contact-based antibody detection (John §1.4 fallback) used by both the viewer and the scoring path:
a non-antigen chain with >= min_contacts heavy-atom pairs within cutoff of the epitope is the antibody.

Run: python tests/test_rfd4_viewer.py
"""
import numpy as np

from episcaf_service.rfd4 import viewer


def _atom(chain, resi, xyz, name="CA", elem="C"):
    return {"name": name, "resn": "ALA", "chain": chain, "resi": resi,
            "x": xyz[0], "y": xyz[1], "z": xyz[2], "elem": elem}


def test_detects_contacting_chain_not_distant_one():
    epi = np.array([[0.0, 0.0, float(i)] for i in range(6)])          # the epitope heavy atoms
    atoms = []
    for i in range(12):                                               # chain B: Fab, hugs the epitope
        atoms.append(_atom("B", i, (1.0, 0.0, float(i % 6))))
    for i in range(12):                                               # chain C: far away, not antibody
        atoms.append(_atom("C", i, (100.0, 100.0, float(i))))
    ab = viewer.detect_antibody_chains(atoms, ["A"], epi)
    assert ab == ["B"], ab


def test_ignores_hydrogens_and_antigen_chain():
    epi = np.array([[0.0, 0.0, float(i)] for i in range(6)])
    atoms = [_atom("A", i, (1.0, 0.0, float(i))) for i in range(12)]  # antigen chain -> never the antibody
    atoms += [_atom("B", i, (1.0, 0.0, float(i % 6)), elem="H") for i in range(12)]  # H-only -> no contacts
    assert viewer.detect_antibody_chains(atoms, ["A"], epi) == []


def test_below_threshold_not_detected():
    epi = np.array([[0.0, 0.0, 0.0]])
    atoms = [_atom("B", i, (1.0, 0.0, 0.0)) for i in range(3)]        # only 3 contacts < min 10
    assert viewer.detect_antibody_chains(atoms, ["A"], epi) == []


if __name__ == "__main__":
    import sys
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    bad = 0
    for fn in fns:
        try:
            fn(); print(f"ok   {fn.__name__}")
        except Exception as e:  # noqa: BLE001
            bad += 1; print(f"FAIL {fn.__name__}: {e}")
    print(f"\n{len(fns)-bad}/{len(fns)} passed")
    sys.exit(1 if bad else 0)
