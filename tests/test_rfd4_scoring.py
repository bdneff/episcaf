#!/usr/bin/env python3
"""Scoring + clash math for the RFD4 engine.

PARITY GATE: the composite is checked against the five worked examples in
RFD4_pipeline_description.pdf section 4.2 (target 0.5 A) — our numbers must equal John's to the
doc's rounding. The clash tests use synthetic geometry to exercise the strict count and the §3.3
shadow rule. (Full clash parity vs John's R still wants a fixture of real designs from him.)

Run:  /Users/bneff/micromamba/envs/bcell-repair/bin/python tests/test_rfd4_scoring.py
"""
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from episcaf_service.rfd4 import clash as C  # noqa: E402
from episcaf_service.rfd4 import geom as G  # noqa: E402
from episcaf_service.rfd4 import score as S  # noqa: E402


def test_sigmoid():
    assert abs(S.sigmoid(5, 5, 1.0) - 0.5) < 1e-12          # at the midpoint
    assert S.sigmoid(0, 5, 1.0) > S.sigmoid(10, 5, 1.0)     # lower-is-better: decreasing
    assert 0.0 <= S.sigmoid(1e6, 5, 12) <= 1e-9             # clipped, no overflow


def test_composite_worked_examples():
    # section 4.2 table: (AFclashes, epitopeRMSD, overallRMSD, meanPAE) -> (base, gate, composite)
    cases = [
        ((0, 0.3, 0.8, 3.0), (0.918, 0.989, 2.897)),
        ((0, 0.45, 1.5, 4.0), (0.823, 0.765, 2.353)),
        ((0, 0.9, 1.0, 3.5), (0.712, 0.000, 0.712)),
        ((3, 0.4, 1.0, 3.5), (0.814, 0.000, 0.814)),
        ((12, 3.0, 4.0, 12.0), (0.022, 0.000, 0.022)),
    ]
    for (clashes, epi, overall, pae), (want_base, want_gate, want_comp) in cases:
        got = S.composite_breakdown(clashes, epi, overall, pae)
        assert abs(got["base"] - want_base) < 2e-3, (clashes, epi, overall, pae, got["base"], want_base)
        assert abs(got["gate"] - want_gate) < 2e-3, (got["gate"], want_gate)
        assert abs(got["composite"] - want_comp) < 2e-3, (got["composite"], want_comp)
    # the passing designs clear 1.8, the failing ones do not
    assert S.passes(2.897) and S.passes(2.353)
    assert not S.passes(0.712) and not S.passes(0.814)


def test_composite_na_safe():
    assert S.composite(None, 0.3, 0.8, 3.0) is None           # missing clash -> NA
    assert S.composite(0, float("nan"), 0.8, 3.0) is None
    assert not S.passes(None)


def test_kabsch_recovers_transform():
    rng = np.random.default_rng(0)
    ref = rng.normal(size=(12, 3))
    # a known rotation + translation applied to ref
    theta = 0.7
    rot = np.array([[np.cos(theta), -np.sin(theta), 0],
                    [np.sin(theta), np.cos(theta), 0], [0, 0, 1]])
    mob = ref @ rot.T + np.array([5.0, -2.0, 1.0])
    assert G.kabsch_rmsd(ref, mob) < 1e-9                     # perfectly superposable
    assert G.kabsch_rmsd(ref, ref) < 1e-12


def _identity_epitope():
    # 4 backbone atoms; using the same coords for native and design gives an identity alignment
    epi = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]], float)
    return epi, epi


def test_clash_strict_and_shadow():
    native_epi, design_epi = _identity_epitope()
    # two scaffold atoms (residues 1 and 2); atom 0 sits 3.0 A from an antibody atom (within 4.0)
    scaffold_xyz = np.array([[10.0, 0, 0], [100.0, 0, 0]], float)
    scaffold_resid = np.array([1, 2])
    antibody_xyz = np.array([[13.0, 0, 0]], float)            # 3.0 A from scaffold atom 0

    # no antigen -> strict only: residue 1 clashes
    shadow, strict = C.count_clashing_residues(design_epi, native_epi, scaffold_xyz, scaffold_resid,
                                               antibody_xyz, antigen_xyz=None)
    assert (shadow, strict) == (1, 1), (shadow, strict)

    # a native antigen atom coincident with scaffold atom 0 (d_ag=0<=1.5) and d_ab=3.0>=2.5 guard
    # -> the shadow rule excuses it: strict still 1, shadow 0
    antigen_xyz = np.array([[10.0, 0, 0]], float)
    shadow, strict = C.count_clashing_residues(design_epi, native_epi, scaffold_xyz, scaffold_resid,
                                               antibody_xyz, antigen_xyz=antigen_xyz)
    assert (shadow, strict) == (0, 1), (shadow, strict)


def test_clash_no_antibody_is_zero():
    native_epi, design_epi = _identity_epitope()
    shadow, strict = C.count_clashing_residues(design_epi, native_epi,
                                               np.array([[1.0, 2, 3]]), np.array([1]),
                                               antibody_xyz=np.empty((0, 3)))
    assert (shadow, strict) == (0, 0)


def test_backbone_prefilter():
    assert C.passes_backbone_prefilter(3) and not C.passes_backbone_prefilter(4)
    assert C.passes_backbone_prefilter(7, method="aom") and not C.passes_backbone_prefilter(8, method="aom")


def main():
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]
    for t in tests:
        t()
        print(f"ok  {t.__name__}")
    print(f"\n{len(tests)} passed")


if __name__ == "__main__":
    main()
