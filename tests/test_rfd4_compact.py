"""Compact input: single-chain, renumbered motif for multi-chain / high-resnum epitopes (John §1.1).
The renumbering and antigen rewrite are pure and tested here; the full RFD4 run needs the cluster.

Run: python tests/test_rfd4_compact.py
"""
import tempfile
from pathlib import Path

from episcaf_service.rfd4 import generate, topology
from episcaf_service.rfd4.viewer import parse_pdb


def _atom(serial, name, resn, chain, resi, x=0.0, elem="C"):
    nm = name if len(name) == 4 else " " + name
    return ("ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s"
            % (serial, nm.ljust(4), resn, chain, resi, x, 0.0, 0.0, elem.rjust(2)))


def test_compact_islands_renumber_with_gap():
    isl = topology.parse_islands("A:205-216;F:623-629")        # len 12 + len 7, two chains
    c = topology.compact_islands(isl, gap=9)
    assert [(x.chain, x.start, x.end) for x in c] == [("A", 1, 12), ("A", 22, 28)]
    # gap of 9 (residues 13-21) separates the islands so they aren't read as one.
    assert c[1].start - c[0].end - 1 == 9


def test_needs_compact_input_rules():
    assert topology.needs_compact_input(topology.parse_islands("A:205-216;F:623-629"))   # multi-chain
    assert topology.needs_compact_input(topology.parse_islands("A:10200-10210"))         # high resnum
    assert not topology.needs_compact_input(topology.parse_islands("A:205-216;A:278-292"))  # single chain, ok


def test_write_compact_antigen_round_trip():
    lines = [
        _atom(1, "N", "ALA", "A", 205), _atom(2, "CA", "ALA", "A", 205),
        _atom(3, "C", "ALA", "A", 205), _atom(4, "O", "ALA", "A", 205),
        _atom(5, "CB", "ALA", "A", 205),                       # a side-chain atom (must survive)
        _atom(6, "N", "GLY", "A", 206), _atom(7, "CA", "GLY", "A", 206),
        _atom(8, "N", "SER", "F", 623), _atom(9, "CA", "SER", "F", 623),
        _atom(10, "N", "LEU", "A", 300),                       # NOT in the epitope -> dropped
    ]
    with tempfile.TemporaryDirectory() as d:
        native = Path(d) / "native.pdb"; native.write_text("\n".join(lines) + "\n")
        out = Path(d) / "compact.pdb"
        compact = generate.write_compact_antigen(native, "A:205-206;F:623", out, gap=9)
        atoms = parse_pdb(out)

    assert [(x.chain, x.start, x.end) for x in compact] == [("A", 1, 2), ("A", 12, 12)]
    assert {a["chain"] for a in atoms} == {"A"}                # single chain
    assert sorted({a["resi"] for a in atoms}) == [1, 2, 12]    # renumbered, gap, epitope only
    assert any(a["name"] == "CB" and a["resi"] == 1 for a in atoms)   # side chain kept on A:205 -> 1
    assert all(a["resn"] != "LEU" for a in atoms)              # non-epitope A:300 dropped


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
