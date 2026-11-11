#!/usr/bin/env python3
"""Unit tests for the AOMclash feature extractor (no data files required).

Run:  python tests/test_aomclash_features.py     (or: pytest tests/)

The exact faithfulness of the port is gated separately against John's own `aom_features.design_feats`
(a one-time cross-check on a synthetic design; see episcaf_analysis/aomclash/features.py). These tests
are the durable, dependency-free guards: the embedded cylinder must reproduce our committed
`native_cylinder_core`, and the feature set must stay frozen.
"""
import os
import sys

import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.spatial import cKDTree

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from episcaf_analysis import native_cylinder_core as ncc          # noqa: E402
from episcaf_analysis.aomclash import features as feats           # noqa: E402

NAMED = [
    "cyl_aware", "cyl_plain", "apr", "n_scaf", "n_seg", "seg_max", "seg_sum",
    "s_p1", "s_p2", "s_p1_free", "s_p1_near",
    "n_p05", "n_p15", "n_p30", "n_p50", "n_p15_free", "n_p30_free",
    "r_sum", "r_sum_free", "r_n30", "r_n50", "r_n70",
    "b_near20", "b_out25", "b_out15", "b_out35_hi", "ag3", "ag4", "ag6", "ag_free_out25",
    "z_mean", "z_p90", "z_max", "r_mean", "ca_rg",
]
EXPECTED_KEYS = set(NAMED) | {f"h{j}" for j in range(35)} | {f"hf{j}" for j in range(35)}
TARGET_KEYS = {"frame_ratio", "n_epi", "n_isl", "ag_out25", "ag_out25_hi", "ag_near15"}


def _synthetic_design(seed=0):
    """A 2-island design + its native antigen, as coordinate arrays (no PDB I/O)."""
    rng = np.random.default_rng(seed)
    segs = [5, 4, 6]
    isl = [("A", 100, 103), ("A", 200, 202)]       # island lengths 4, 3
    pos, ep_pos = 1, []
    for k, (_, a, b) in enumerate(isl):
        pos += segs[k]
        ep_pos.append(list(range(pos, pos + b - a + 1)))
        pos += b - a + 1
    epi_pos = [p for rg in ep_pos for p in rg]
    epi_set = set(epi_pos)
    n_res = pos - 1 + segs[-1]

    atoms = []
    for r in range(1, n_res + 1):
        names = ("N", "CA", "C", "O") if r in epi_set else ("N", "CA", "C", "O", "CB")
        for nm in names:
            x, y, z = rng.normal(0, 12, 3)
            atoms.append((nm, "A", r, x, y, z))

    n_epi_res = len(epi_pos)
    ref = rng.normal(0, 10, (n_epi_res * 4, 3))
    epi_ca = rng.normal(0, 6, (n_epi_res, 3))
    ag_ca = np.vstack([epi_ca, rng.normal(0, 10, (30, 3))])
    o, Rf, _ = feats.frame_from(epi_ca, ag_ca)
    ag_L = (ag_ca - o) @ Rf.T
    ag_tree = cKDTree(ag_L)
    epi_L = (epi_ca - o) @ Rf.T
    field = rng.random((33, 33, 33)).astype(np.float32)
    P1, P2 = gaussian_filter(field, 1.0), gaussian_filter(field, 2.5)
    return dict(atoms=atoms, epi_pos=epi_pos, segs=segs, ref=ref, o=o, Rf=Rf,
                ag_tree=ag_tree, epi_L=epi_L, P1=P1, P2=P2, epi_ca=epi_ca, ag_ca=ag_ca, n_isl=len(isl))


def test_cylinder_counts_matches_native_core():
    """cylinder_counts must reproduce native_cylinder_core.native_aware_scaffold_count exactly."""
    for seed in range(6):
        rng = np.random.default_rng(100 + seed)
        epi_ca = np.array([[0, 0, 0], [6, 0, 0], [0, 6, 0], [6, 6, 0]], float)
        scaf_ca = rng.normal([3, 3, 10], [10, 10, 14], size=(40, 3))
        native_heavy = rng.normal([3, 3, 2], [8, 8, 6], size=(60, 3))
        # match the normal-orientation reference cylinder_counts uses internally
        base, normal = ncc.cylinder_frame(epi_ca, np.vstack([epi_ca, scaf_ca]))
        aware_ref, plain_ref = ncc.native_aware_scaffold_count(
            scaf_ca, base, normal, native_heavy, exclude_dist=1.0)
        aware, plain = feats.cylinder_counts(scaf_ca, epi_ca, cKDTree(native_heavy))
        assert (aware, plain) == (aware_ref, plain_ref), (seed, (aware, plain), (aware_ref, plain_ref))


def test_design_feats_key_set_frozen():
    d = _synthetic_design(0)
    f = feats.design_feats(d["atoms"], d["epi_pos"], d["segs"], d["ref"], d["o"], d["Rf"],
                           d["ag_tree"], d["epi_L"], d["P1"], d["P2"])
    assert f is not None
    assert set(f) == EXPECTED_KEYS, set(f) ^ EXPECTED_KEYS
    assert len(f) == 105
    assert all(np.isfinite(float(v)) for v in f.values())


def test_design_feats_deterministic():
    d = _synthetic_design(1)
    a = feats.design_feats(d["atoms"], d["epi_pos"], d["segs"], d["ref"], d["o"], d["Rf"],
                           d["ag_tree"], d["epi_L"], d["P1"], d["P2"])
    b = feats.design_feats(d["atoms"], d["epi_pos"], d["segs"], d["ref"], d["o"], d["Rf"],
                           d["ag_tree"], d["epi_L"], d["P1"], d["P2"])
    assert a.keys() == b.keys()
    assert all(a[k] == b[k] for k in a)


def test_design_feats_scaffold_structure():
    """n_scaf / n_seg / seg_max / seg_sum must reflect the passed layout, epitope excluded."""
    d = _synthetic_design(2)
    f = feats.design_feats(d["atoms"], d["epi_pos"], d["segs"], d["ref"], d["o"], d["Rf"],
                           d["ag_tree"], d["epi_L"], d["P1"], d["P2"])
    assert f["n_scaf"] == 15 and f["n_seg"] == 3 and f["seg_max"] == 6 and f["seg_sum"] == 15
    assert abs(f["apr"] - 5.0) < 1e-9          # 5 heavy atoms per scaffold residue


def test_design_feats_missing_epitope_backbone_returns_none():
    d = _synthetic_design(3)
    atoms = [a for a in d["atoms"] if not (a[2] in set(d["epi_pos"]) and a[0] == "CA")]  # drop epitope CA
    f = feats.design_feats(atoms, d["epi_pos"], d["segs"], d["ref"], d["o"], d["Rf"],
                           d["ag_tree"], d["epi_L"], d["P1"], d["P2"])
    assert f is None


def test_target_features_keys_and_undefined_frame():
    d = _synthetic_design(4)
    t = feats.target_features(d["epi_ca"], d["ag_ca"], d["ag_ca"], d["n_isl"])
    assert t is not None
    o, Rf, ag_tree, epi_L, tf = t
    assert set(tf) == TARGET_KEYS, set(tf) ^ TARGET_KEYS
    assert tf["n_epi"] == len(d["epi_ca"]) and tf["n_isl"] == 2
    # fewer than 20 local antigen CA -> frame undefined
    assert feats.frame_from(d["epi_ca"], d["epi_ca"]) is None


def test_read_pdb_atoms_roundtrip(tmp_path=None):
    """read_pdb_atoms parses the ATOM columns John's extractor relies on."""
    import tempfile
    line = "ATOM      2  CA  ALA A  14      12.345 -6.789   3.210  1.00  0.00"
    with tempfile.NamedTemporaryFile("w", suffix=".pdb", delete=False) as fh:
        fh.write(line + "\n")
        path = fh.name
    try:
        atoms = feats.read_pdb_atoms(path)
        assert len(atoms) == 1
        name, chain, resid, x, y, z = atoms[0]
        assert name == "CA" and chain == "A" and resid == 14
        assert abs(x - 12.345) < 1e-6 and abs(y + 6.789) < 1e-6 and abs(z - 3.210) < 1e-6
    finally:
        os.unlink(path)


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print(f"  ok: {name}")
    print("all aomclash feature tests passed")
