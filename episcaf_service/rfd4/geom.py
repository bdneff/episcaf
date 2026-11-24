"""Kabsch superposition primitives, shared by the clash filter and the RMSD metrics.

Ported from the Python reference in RFD4_pipeline_description.pdf section 3.6 (kabsch) and the
structure-metric functions in tamarind.R (compute/apply_kabsch_transform, calculate_rmsd).

Convention: kabsch(ref, mob) returns (R, ref_center, mob_center) such that a mobile-frame point
set X maps into the reference frame by  apply_transform(X, R, ref_center, mob_center).
"""
from __future__ import annotations

import numpy as np


def kabsch(ref: np.ndarray, mob: np.ndarray):
    """Optimal rotation aligning `mob` onto `ref` (both (n,3)), with reflection correction (det)."""
    ref = np.asarray(ref, float)
    mob = np.asarray(mob, float)
    if ref.shape != mob.shape or ref.ndim != 2 or ref.shape[1] != 3:
        raise ValueError(f"kabsch needs matching (n,3) arrays, got {ref.shape} and {mob.shape}")
    rc = ref.mean(axis=0)
    mc = mob.mean(axis=0)
    h = (mob - mc).T @ (ref - rc)
    u, _s, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    rot = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    return rot, rc, mc


def apply_transform(points: np.ndarray, rot: np.ndarray, ref_center: np.ndarray,
                    mob_center: np.ndarray) -> np.ndarray:
    """Map mobile-frame `points` into the reference frame: (points - mc) @ R.T + rc."""
    points = np.asarray(points, float)
    return (points - mob_center) @ rot.T + ref_center


def kabsch_rmsd(ref: np.ndarray, mob: np.ndarray) -> float:
    """RMSD of `mob` superposed onto `ref` by Kabsch (both (n,3), paired row by row)."""
    rot, rc, mc = kabsch(ref, mob)
    aligned = apply_transform(mob, rot, rc, mc)
    diff = aligned - np.asarray(ref, float)
    return float(np.sqrt((diff * diff).sum(axis=1).mean()))
