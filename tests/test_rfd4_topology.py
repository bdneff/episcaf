#!/usr/bin/env python3
"""Topology math for the RFD4 engine: island parsing, the even scaffold split, and parsing
RFD4's motif_residue_map back into a per-design layout. Pure logic, no cluster/GPU.

Parity note: the even-split and the map parse are ported from tamarind.R; the (26,25,25) and
map cases below are the reference values to hold against John's R.

Run:  python3 tests/test_rfd4_topology.py
"""
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from episcaf_service.rfd4 import topology as T  # noqa: E402


def test_parse_islands():
    isl = T.parse_islands("A:205-216;A:278-292")
    assert [(i.chain, i.start, i.end, i.length) for i in isl] == [
        ("A", 205, 216, 12), ("A", 278, 292, 15)], isl
    one = T.parse_islands("A:252")
    assert len(one) == 1 and one[0].length == 1, one
    # whitespace + multi-char chain id tolerated
    assert T.parse_islands(" H : 10 - 12 ")[0].length == 3
    for bad in ("", "205-216", "A:300-299"):
        try:
            T.parse_islands(bad)
        except ValueError:
            pass
        else:
            raise AssertionError(f"expected ValueError for {bad!r}")


def test_even_split():
    assert T.even_split(76, 3) == (26, 25, 25)
    assert T.even_split(76, 1) == (76,)
    assert T.even_split(10, 4) == (3, 3, 2, 2)
    assert sum(T.even_split(97, 5)) == 97
    assert T.even_split(0, 2) == (0, 0)


def test_plan_topology_two_islands():
    top = T.plan_topology("A:205-216;A:278-292")        # lengths 12, 15 -> epitope 27, scaffold 76
    assert top.scaffold_segments == (26, 25, 25), top.scaffold_segments
    assert top.total_epitope_length == 27
    # layout: 26 gap | 12 epi (27..38) | 25 gap | 15 epi (64..78) | 25 gap
    assert top.epitope_positions == tuple(range(27, 39)) + tuple(range(64, 79)), top.epitope_positions
    assert len(top.scaffold_positions) == 76
    assert max(top.epitope_positions) <= top.protein_length
    # epitope and scaffold positions partition 1..103
    assert set(top.epitope_positions) | set(top.scaffold_positions) == set(range(1, 104))


def test_plan_topology_single_island():
    top = T.plan_topology("A:1-7")                       # len 7 -> scaffold 96 over 2 gaps
    assert top.scaffold_segments == (48, 48)
    assert top.epitope_positions == tuple(range(49, 56))


def test_plan_topology_rejects_oversized_epitope():
    try:
        T.plan_topology("A:1-103")                       # no room for a scaffold
    except ValueError:
        pass
    else:
        raise AssertionError("expected ValueError for an epitope filling the protein")


def test_topology_from_map():
    # native 312-328 (len17) -> design 28-44; native 362-364 (len3) -> design 83-85
    out = T.topology_from_map("R312-328>R28-44;R362-364>R83-85", total_epitope_length=20)
    assert out["epitope_positions"] == tuple(range(28, 45)) + tuple(range(83, 86))
    # segments: lead 27, middle 83-44-1=38, trail 103-85=18 ; sums with epitope to 103
    assert out["scaffold_segments"] == (27, 38, 18)
    assert sum(out["scaffold_segments"]) + 20 == 103


def test_topology_from_map_reordered():
    # unindexed run can place the second native island first in the design
    out = T.topology_from_map("A50-52>A10-12;A10-14>A40-44", total_epitope_length=8)
    assert out["epitope_positions"] == tuple(range(10, 13)) + tuple(range(40, 45))


def test_topology_from_map_rejects_bad():
    cases = [
        ("A1-5>A1-3", 5),                 # source/dest length mismatch
        ("A1-5>A1-5", 10),                # wrong total epitope length
        ("A1-5>A10-14;A3-4>A12-13", 7),   # overlapping dest ranges
        ("garbage", 5),
    ]
    for m, n in cases:
        try:
            T.topology_from_map(m, total_epitope_length=n)
        except T.BadResidueMap:
            pass
        else:
            raise AssertionError(f"expected BadResidueMap for {m!r}")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} passed")


if __name__ == "__main__":
    main()
