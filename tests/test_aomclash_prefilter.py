#!/usr/bin/env python3
"""Unit tests for the AOMclash backbone pre-filter's contig parser + keep rule (no structures)."""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import aomclash_prefilter as PF   # noqa: E402


def test_design_epi_positions_single_island():
    # 55 scaffold, then native A377-419 (43 res) at output 56..98, then 5 scaffold
    pos = PF.design_epi_positions("55-55/A377-419/5-5")
    assert pos == list(range(56, 99)), pos
    assert len(pos) == 43


def test_design_epi_positions_multi_island_multichain():
    # 10 scaffold / A5-9 (5) at 11..15 / 6 linker / B20-25 (6) at 22..27 / 8 scaffold
    pos = PF.design_epi_positions("10-10/A5-9/6-6/B20-25/8-8")
    assert pos == list(range(11, 16)) + list(range(22, 28)), pos


def test_design_epi_positions_offset_and_single_residue():
    assert PF.design_epi_positions("4-4/C100/3-3") == [5]            # single-residue motif 'C100'
    assert PF.design_epi_positions("55-55/A377-419/5-5", resid_offset=0) == list(range(55, 98))


def test_keep_mask_fraction_per_group():
    rows = [dict(id="e1", s_p1_free=v) for v in (1.0, 2.0, 3.0, 4.0)] + \
           [dict(id="e2", s_p1_free=v) for v in (5.0, 6.0)] + \
           [dict(id="e1", s_p1_free=None)]                            # unscored -> never kept
    PF.keep_mask(rows, "id", keep_frac=0.5, cutoff=None)
    kept = {(r["id"], r["s_p1_free"]) for r in rows if r["kept"]}
    assert kept == {("e1", 1.0), ("e1", 2.0), ("e2", 5.0)}, kept     # lowest half per epitope
    assert all(not r["kept"] for r in rows if r["s_p1_free"] is None)


def test_keep_mask_cutoff():
    rows = [dict(id="e", s_p1_free=v) for v in (1.0, 5.0, 9.0)]
    PF.keep_mask(rows, "id", keep_frac=0.5, cutoff=7.0)
    assert [r["kept"] for r in rows] == [True, True, False]


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"  ok: {name}")
    print("all aomclash prefilter tests passed")
