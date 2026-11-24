"""Antibody-clash filter — the strict count and John's antigen-occupancy "shadow" rule (§3.3).

Ported from the §3.6 Python reference and clash_method.R. A design's epitope backbone is Kabsch-
aligned onto the native epitope; the SAME transform is applied to the scaffold heavy atoms, which are
then tested against the antibody. A scaffold residue clashes if any of its heavy atoms lies within
`cutoff` (4.0 Å) of an antibody heavy atom — EXCEPT, under the shadow rule (John's recent §3.3 fix
that removes over-called clashes), an atom that sits within `shadow` (1.5 Å) of a native antigen atom
AND is at least `guard` (2.5 Å) from the antibody is not counted: the native complex already places an
atom there next to the antibody, so the position is compatible.

Used twice: on RFD4 backbones as a pre-filter (drop >threshold before AF2) and on the AF2 model as the
composite's `AFclashes` term. `AFclashes_strict` is the same count with the shadow rule off.
"""
from __future__ import annotations

import numpy as np
from scipy.spatial.distance import cdist

from episcaf_service.rfd4.config import DEFAULT, ClashConfig
from episcaf_service.rfd4.geom import apply_transform, kabsch


def _min_dist(points: np.ndarray, others: np.ndarray) -> np.ndarray:
    """Per-point minimum distance to any row of `others`. inf where `others` is empty."""
    points = np.asarray(points, float)
    if others is None or len(others) == 0:
        return np.full(len(points), np.inf)
    return cdist(points, np.asarray(others, float)).min(axis=1)


def count_clashing_residues(design_epitope_bb: np.ndarray, native_epitope_bb: np.ndarray,
                            scaffold_xyz: np.ndarray, scaffold_resid, antibody_xyz: np.ndarray,
                            antigen_xyz: np.ndarray | None = None,
                            cfg: ClashConfig = DEFAULT.clash) -> tuple[int, int]:
    """Return (shadow_count, strict_count): clashing scaffold *residues* against the antibody.

    design_epitope_bb / native_epitope_bb: matched (m,3) epitope backbone atoms (N,CA,C,O), paired
      in the same order, defining the alignment. scaffold_xyz: (a,3) scaffold heavy atoms in the
      design frame; scaffold_resid: length-a residue id per atom. antigen_xyz: native non-antibody
      heavy atoms for the shadow rule (None/empty -> strict only).
    """
    scaffold_xyz = np.asarray(scaffold_xyz, float)
    scaffold_resid = np.asarray(scaffold_resid)
    if len(scaffold_xyz) == 0:
        return 0, 0

    rot, rc, mc = kabsch(np.asarray(native_epitope_bb, float), np.asarray(design_epitope_bb, float))
    xyz = apply_transform(scaffold_xyz, rot, rc, mc)

    d_ab = _min_dist(xyz, antibody_xyz)
    near = d_ab <= cfg.af_clash_cutoff

    strict_resids = set(scaffold_resid[near].tolist())

    if antigen_xyz is not None and len(antigen_xyz) > 0:
        d_ag = _min_dist(xyz, antigen_xyz)
        compatible = (d_ag <= cfg.shadow) & (d_ab >= cfg.guard)
        shadow_near = near & ~compatible
    else:
        shadow_near = near
    shadow_resids = set(scaffold_resid[shadow_near].tolist())

    return len(shadow_resids), len(strict_resids)


def passes_backbone_prefilter(clash_count: int, cfg: ClashConfig = DEFAULT.clash,
                              method: str = "antibody") -> bool:
    """Keep an RFD4 backbone for AF2 iff its clash count is within the pre-filter threshold."""
    threshold = cfg.backbone_threshold_aom if method == "aom" else cfg.backbone_threshold_antibody
    return clash_count <= threshold
