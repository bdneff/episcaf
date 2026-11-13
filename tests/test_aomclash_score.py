#!/usr/bin/env python3
"""Sanity test for the AOMclash scorer glue (episcaf_analysis/aomclash/score.py).

Run: /Users/bneff/micromamba/envs/bcell-repair/bin/python tests/test_aomclash_score.py
Needs numpy/scipy/biopython + the vendored map + docs/pdb/8PWW.pdb; skips cleanly otherwise.
"""
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

HERE = os.path.dirname(__file__)
PDB = os.path.join(HERE, "..", "docs", "pdb", "8PWW.pdb")
MAP = os.path.join(HERE, "..", "episcaf_analysis", "aomclash", "aom_map.pkl")


def test_score_native_as_design():
    if not (os.path.exists(PDB) and os.path.exists(MAP)):
        print("  skip: need docs/pdb/8PWW.pdb + aom_map.pkl")
        return
    from episcaf_analysis.aomclash import score as S

    isl = [("A", 205, 205), ("A", 209, 209), ("A", 212, 213), ("A", 216, 216), ("A", 278, 278),
           ("A", 281, 282), ("A", 284, 286), ("A", 288, 289), ("A", 292, 292)]
    tgt = S.build_target(PDB, isl, exclude_pdb="8pww")
    assert tgt is not None and tgt.n_complexes >= 1000
    assert tgt.P1.shape == (33, 33, 33) and len(tgt.ref) == 14 * 4   # 14 epitope residues x N/CA/C/O

    ep = S.native_epi_pos(isl)
    feat = S.score_design(tgt, PDB, ep)
    assert feat is not None
    # the native antigen's own atoms are inside antigen space (dag<=4 => free=False), so the carved
    # occupancy sum is ~0, while the uncarved s_p1 is positive -- the free mask is doing its job
    assert feat["s_p1_free"] < 1.0, feat["s_p1_free"]
    assert feat["s_p1"] > feat["s_p1_free"]
    assert S.aomclash_score(tgt, PDB, ep) == feat["s_p1_free"]

    # epitope string parsing round-trips to the same spans
    assert S.parse_epitope("A:205-216;A:278-292") == [("A", 205, 216), ("A", 278, 292)]
    print(f"  ok: test_score_native_as_design (s_p1_free={feat['s_p1_free']:.3f}, s_p1={feat['s_p1']:.2f})")


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
    print("all aomclash score tests passed")
