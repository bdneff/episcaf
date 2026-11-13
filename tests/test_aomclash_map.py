#!/usr/bin/env python3
"""Sanity tests for the vendored antibody-occupancy map (episcaf_analysis/aomclash/aom_map.pkl).

Run: /Users/bneff/micromamba/envs/bcell-repair/bin/python tests/test_aomclash_map.py
Needs numpy + scipy; skips cleanly if the map artifact is absent.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from episcaf_analysis.aomclash import aom_map as M   # noqa: E402

MAP = os.path.join(os.path.dirname(__file__), "..", "episcaf_analysis", "aomclash", "aom_map.pkl")


def test_map_loads_and_is_physical():
    if not os.path.exists(MAP):
        print("  skip: aom_map.pkl not present")
        return
    art = M.load_map(MAP)
    assert art["field_side"] == 33 and art["field_length"] == 33 ** 3
    assert len(art["keys"]) == art["fields_packed"].shape[0] == len(art["pdbs"])

    field, n = M.mean_field(art)
    assert field.shape == (33, 33, 33)
    assert n == art["fields_packed"].shape[0]
    assert 0.0 <= field.min() and field.max() <= 1.0            # it is a mean of booleans

    P1, P2 = M.smoothed_maps(field)
    assert P1.shape == P2.shape == (33, 33, 33)
    assert 0.0 <= P1.min() and P1.max() <= 1.0
    # physical: occupancy concentrates on the antibody side (z>0; axis 2, center index 16)
    zprof = [P1[:, :, k].max() for k in range(33)]
    assert int(np.argmax(zprof)) > 16, "occupancy should peak on the z>0 antibody side"
    # P2 is smoother than P1 (lower peak)
    assert P2.max() <= P1.max() + 1e-6
    print(f"  ok: test_map_loads_and_is_physical ({n} complexes, peak z-index {int(np.argmax(zprof))})")


def test_leave_one_out_excludes():
    if not os.path.exists(MAP):
        print("  skip: aom_map.pkl not present")
        return
    art = M.load_map(MAP)
    full, nf = M.mean_field(art)
    pdb0 = str(art["pdbs"][0])
    excl, ne = M.mean_field(art, exclude_pdb=pdb0)
    assert ne < nf
    assert not np.allclose(full, excl)
    print(f"  ok: test_leave_one_out_excludes (exclude {pdb0}: {nf} -> {ne})")


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    print("all aomclash map tests passed")
