"""growPatch core tests.

Two layers:
  1. Self-contained synthetic tests of each deterministic unit (islands, gap-fill, trim, jaccard,
     NMS dedupe, set-cover, grow invariants) -- no external data, CI-safe.
  2. A guarded byte-parity cross-check against John Altin's reference `growPatch.py` on a real
     antigen chain. Runs only when GROWPATCH_REF (path to his growPatch.py) and GROWPATCH_CIF
     (path to a .cif) are set; otherwise skipped. This is the growPatch analogue of the AOMclash
     byte-for-byte feature check.

Run: python tests/test_growpatch_core.py
"""
import importlib.util
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from episcaf_analysis.growpatch import core  # noqa: E402


# ------------------------------------------------------------------ helpers

def _line_residues(xs):
    """Residues at given x positions on a line; two atoms each, rep = their mean."""
    res = []
    for i, x in enumerate(xs):
        atoms = np.array([[x, 0.0, 0.0], [x, 0.8, 0.0]], dtype=float)
        res.append(dict(chain="A", resno=i + 1, key=("A", i + 1),
                        resname="ALA", atoms=atoms, rep=atoms.mean(0)))
    return res


def _ok(name):
    print(f"  ok: {name}")


# ------------------------------------------------------------------ 1. synthetic units

def test_islands():
    segs = core.islands({("A", 1), ("A", 2), ("A", 3), ("A", 7), ("A", 8)})
    assert segs == [("A", 1, 3), ("A", 7, 8)], segs
    # never bridges chains, sorted by chain
    segs2 = core.islands({("B", 5), ("A", 1), ("A", 2)})
    assert segs2 == [("A", 1, 2), ("B", 5, 5)], segs2
    _ok("test_islands")


def test_gap_fill_merges_and_respects_present():
    keys = {("A", 1), ("A", 2), ("A", 3), ("A", 7), ("A", 8)}
    present_all = {"A": frozenset(range(1, 9))}
    # gap 3->7 is 4 <= gap+1 (11), so it merges and span-fills 1..8
    assert core.gap_fill(keys, 10, present_all) == {("A", n) for n in range(1, 9)}
    # a residue missing from `present` is not filled, so the island breaks there
    present_hole = {"A": frozenset({1, 2, 3, 4, 6, 7, 8})}
    gf = core.gap_fill(keys, 10, present_hole)
    assert gf == {("A", n) for n in (1, 2, 3, 4, 6, 7, 8)}, gf
    assert core.islands(gf) == [("A", 1, 4), ("A", 6, 8)]
    # a gap larger than gap+1 is NOT bridged
    far = {("A", 1), ("A", 2), ("A", 20)}
    assert core.gap_fill(far, 10, {"A": frozenset(range(1, 21))}) == {("A", 1), ("A", 2), ("A", 20)}
    _ok("test_gap_fill_merges_and_respects_present")


def test_trim_short_islands():
    keys = {("A", n) for n in range(1, 6)} | {("A", 20), ("A", 21)}   # len-5 island + len-2 island
    # min_len 3, nothing protected -> the len-2 island is dropped
    assert core.trim_short_islands(keys, 3, frozenset()) == {("A", n) for n in range(1, 6)}
    # protecting a residue in the short island keeps it
    kept = core.trim_short_islands(keys, 3, frozenset({("A", 20)}))
    assert ("A", 20) in kept and ("A", 21) in kept
    _ok("test_trim_short_islands")


def test_jaccard():
    assert core.jaccard({1, 2, 3}, {2, 3, 4}) == 2 / 4
    assert core.jaccard(set(), set()) == 0.0
    assert core.jaccard({1}, set()) == 0.0
    _ok("test_jaccard")


def test_grow_patch_contiguous_interval():
    # strictly irregular spacing (all <= contact apart, all pairwise distances distinct) so there
    # are no ties: a grown patch is then a deterministic contiguous index interval containing the seed
    xs = np.cumsum([0, 3.6, 4.0, 3.8, 4.4, 3.7, 4.1, 3.9, 4.3, 3.5])
    res = _line_residues(xs)
    adj = core.build_adjacency(res, contact=4.5)
    # consecutive residues are adjacent; non-consecutive are not
    assert adj[3] == {2, 4}, adj[3]
    surf = frozenset(range(len(res)))
    for seed in (0, 4, len(res) - 1):
        p = sorted(core.grow_patch(seed, surf, res, adj, size=5, link=10.0))
        assert len(p) == 5
        assert p == list(range(p[0], p[0] + 5)), p          # contiguous interval
        assert seed in p
    _ok("test_grow_patch_contiguous_interval")


def test_nms_dedupe_order_and_suppression():
    big = frozenset(("A", n) for n in range(1, 13))             # 12 res, 1 island
    near = frozenset(("A", n) for n in range(1, 12))            # 11 res, Jaccard 11/12 ~ .92 vs big
    far = frozenset(("A", n) for n in range(40, 52))            # 12 res, disjoint
    two_isl = frozenset(list(("A", n) for n in range(1, 7)) + list(("A", n) for n in range(60, 66)))
    cands = [(set(), near, 1), (set(), big, 1), (set(), far, 1), (set(), two_isl, 2)]
    kept = core.nms_dedupe(cands, nms_jaccard=0.8)
    kept_sets = [k[1] for k in kept]
    assert big in kept_sets and far in kept_sets               # big beats its near-duplicate
    assert near not in kept_sets                               # suppressed (Jaccard >= 0.8)
    assert two_isl in kept_sets                                # different footprint, kept
    # fewer islands sort first; among equal islands, larger first
    assert kept[-1][2] == 2, "two-island candidate should sort last"
    _ok("test_nms_dedupe_order_and_suppression")


def test_greedy_set_cover():
    target = {("A", n) for n in range(1, 11)}
    patch_sets = [
        ("p1", {("A", n) for n in range(1, 7)}),               # covers 1-6
        ("p2", {("A", n) for n in range(5, 11)}),              # covers 5-10
        ("p3", {("A", 1), ("A", 2)}),                          # redundant
    ]
    chosen, uncovered = core.greedy_set_cover(target, patch_sets)
    assert not uncovered
    assert [c[0] for c in chosen] == ["p1", "p2"] or [c[0] for c in chosen] == ["p2", "p1"]
    # an uncoverable residue is reported
    chosen2, unc2 = core.greedy_set_cover(target | {("A", 99)}, patch_sets)
    assert unc2 == {("A", 99)}
    _ok("test_greedy_set_cover")


def test_candidate_patches_end_to_end_small():
    # a compact blob of residues on a line: one 16-seed grow, gap-fill, keep if >= MIN_RES
    xs = np.cumsum([0.0] + [4.0] * 19)                          # 20 residues, all adjacent
    res = _line_residues(xs)
    adj = core.build_adjacency(res, contact=4.5)
    idx_key = {i: r["key"] for i, r in enumerate(res)}
    present = {"A": frozenset(n for (_, n) in idx_key.values())}
    surf = frozenset(range(len(res)))
    kept = core.candidate_patches(res, adj, surf, idx_key, present,
                                  size=16, min_res=12, max_islands=3)
    assert kept, "expected at least one candidate"
    for raw_r, gf_r, nisl in kept:
        assert len(gf_r) >= 12 and nisl <= 3
        assert gf_r == core.gap_fill(raw_r, core.GAP, present) or gf_r <= core.gap_fill(raw_r, core.GAP, present)
    _ok("test_candidate_patches_end_to_end_small")


# ------------------------------------------------------------------ 2. guarded cross-check vs John's code

def test_parity_against_reference():
    ref_path, cif = os.environ.get("GROWPATCH_REF"), os.environ.get("GROWPATCH_CIF")
    if not (ref_path and cif and os.path.exists(ref_path) and os.path.exists(cif)):
        print("  skip: test_parity_against_reference (set GROWPATCH_REF + GROWPATCH_CIF)")
        return
    chain = os.environ.get("GROWPATCH_CHAIN", "A")
    spec = importlib.util.spec_from_file_location("growPatch_ref", ref_path)
    gp = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(gp)

    res = gp.parse_chain(cif, chain)                            # his reader -> residue dicts
    idx_key = {i: r["key"] for i, r in enumerate(res)}
    present = {chain: frozenset(n for (_, n) in idx_key.values())}

    # adjacency identical
    adj_ours = core.build_adjacency(res, 4.5)
    adj_ref = gp.build_adjacency(res, 4.5)
    assert adj_ours == adj_ref, "adjacency differs from reference"

    surf = frozenset(range(len(res)))
    # grow identical from a spread of seeds
    for seed in range(0, len(res), max(1, len(res) // 25)):
        assert core.grow_patch(seed, surf, res, adj_ours, 16, 5.0) == \
               gp.grow_patch(seed, surf, res, adj_ref, 16, 5.0), f"grow differs at seed {seed}"

    # gap_fill / islands / trim identical on the grown patches
    for seed in range(0, len(res), max(1, len(res) // 15)):
        raw = {idx_key[i] for i in core.grow_patch(seed, surf, res, adj_ours, 16, 5.0)}
        assert core.gap_fill(raw, 10, present) == gp.gap_fill(raw, 10, present)
        g = core.gap_fill(raw, 10, present)
        assert core.islands(g) == gp.islands(g)
        assert core.trim_short_islands(g, 3, frozenset()) == gp.trim_short_islands(g, 3, frozenset())

    # full candidate enumeration identical (our composed function vs his main()-loop primitives)
    ours = core.candidate_patches(res, adj_ours, surf, idx_key, present,
                                  size=16, min_res=12, max_islands=3, min_island_len=3)
    raw_ref = {p for s in surf for p in [gp.grow_patch(s, surf, res, adj_ref, 16, 5.0)] if len(p) >= 3}
    cands_ref = []
    for p in raw_ref:
        raw_r = {idx_key[i] for i in p}
        gf = gp.trim_short_islands(gp.gap_fill(raw_r, 10, present), 3, frozenset())
        if gf and len(gf) >= 12 and len(gp.islands(gf)) <= 3:
            cands_ref.append((raw_r, gf, len(gp.islands(gf))))
    cands_ref.sort(key=lambda t: (t[2], -len(t[1]), sorted(t[1])))
    kept_ref = []
    for raw_r, gf, nisl in cands_ref:
        if all(gp.jaccard(gf, k[1]) < 0.8 for k in kept_ref):
            kept_ref.append((raw_r, gf, nisl))
    assert [k[1] for k in ours] == [k[1] for k in kept_ref], "candidate sets differ from reference"
    print(f"  ok: test_parity_against_reference ({len(ours)} candidates on {os.path.basename(cif)}:{chain})")


if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for fn in fns:
        fn()
    print("all growpatch core tests passed")
