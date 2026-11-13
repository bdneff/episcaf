"""Antibody's-eye-view epitope enumerator -- the direction-sampling alternative to seed-and-grow.

An antibody approaches an antigen from *outside*, along some direction, and contacts the cap of the
surface that faces it. So instead of growing a patch along a residue-adjacency graph (growPatch) or
enclosing residues in a ball (surfBall), this enumerates over approach *directions*: sample directions
roughly uniformly over the sphere, and for each one take the surface residues on the facing cap. Each
direction yields one candidate conformational epitope.

Why this shape:
  * physically motivated -- a candidate is literally "what an antibody coming from here would touch";
  * uniform coverage -- directions are evenly spaced (Fibonacci sphere), so flat regions are sampled
    as well as residue-dense ones (residue-seeded growth over-samples dense regions);
  * caps are roughly planar contact patches, the shape a paratope presents.

This is an EVALUATION alternative, benchmarked head-to-head against growPatch's seed-and-grow on the
known-epitope sample (`scripts/growpatch_benchmark.py`). It reuses the exact same post-processing as
growPatch -- gap-10 fill, short-island trim, >=min_res / <=max_islands footprint filter, Jaccard NMS
(`core`) -- so the ONLY thing that differs is how a raw patch is formed. That keeps the comparison
about the enumeration geometry, not the cleanup.

The method, given a read `Antigen` and its surface set:
  1. interior reference point c0 = centroid of all residue representative points.
  2. for each of `n_dirs` directions d (Fibonacci sphere): project every surface residue onto d; the
     apex is the residue of greatest height (the first the antibody reaches). The raw patch is the
     surface residues within `footprint_radius` of the apex laterally (perpendicular to d) and no
     more than `depth` below the apex height (so a residue on the far side, laterally close but deep,
     is not pulled in). That is the facing cap.
  3. the same core cleanup + NMS as growPatch.
"""
from __future__ import annotations

import numpy as np

from . import core


def fibonacci_sphere(n):
    """`n` roughly-uniform unit vectors on the sphere (Fibonacci spiral)."""
    i = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)                 # polar angle, evenly in cos
    golden = np.pi * (1.0 + 5.0 ** 0.5)                # golden angle
    theta = golden * i
    return np.column_stack([np.sin(phi) * np.cos(theta),
                            np.sin(phi) * np.sin(theta),
                            np.cos(phi)])


def enumerate_eyeview(antigen, sasa_cut=core.SASA_CUT, n_dirs=150, footprint_radius=12.0, depth=12.0,
                      gap=core.GAP, min_res=core.MIN_RES, max_islands=core.MAX_ISLANDS,
                      min_island_len=core.MIN_ISLAND_LEN, nms_jaccard=core.NMS_JACCARD,
                      protect=frozenset()):
    """Enumerate candidate epitopes by direction-sampled facing caps on a read `Antigen`.

    Returns the kept candidates as the same dicts `read.enumerate_patches` returns
    ({patch, n_raw, n_res, n_islands, epitope, islands, keys, vmd_selection}), NMS-ordered. `protect`
    (a set of (chain, slot) keys) shields a known epitope from short-island trimming (benchmark use).
    """
    res = antigen.residues
    reps = np.array([r["rep"] for r in res], dtype=float)
    if len(reps) == 0:
        return []
    c0 = reps.mean(axis=0)

    surf = sorted(antigen.surface_indices(sasa_cut))
    if len(surf) < 3:
        return []
    surf = np.array(surf)
    V = reps[surf] - c0                                 # surface reps relative to the interior point

    raw_patches = set()
    for d in fibonacci_sphere(n_dirs):
        h = V @ d                                       # height of each surface residue along d
        a = int(h.argmax())                             # apex = first residue the antibody reaches
        lat = V - np.outer(h, d)                        # lateral component (perp to d)
        dl = np.linalg.norm(lat - lat[a], axis=1)       # lateral distance from the apex
        mask = (dl <= footprint_radius) & (h >= h[a] - depth)
        p = frozenset(int(surf[j]) for j in np.nonzero(mask)[0])
        if len(p) >= 3:
            raw_patches.add(p)

    cands = []
    for p in raw_patches:
        raw_r = {antigen.idx_key[i] for i in p}
        gf = core.trim_short_islands(core.gap_fill(raw_r, gap, antigen.present), min_island_len, protect)
        if not gf:
            continue
        nisl = len(core.islands(gf))
        if len(gf) >= min_res and nisl <= max_islands:
            cands.append((raw_r, gf, nisl))
    kept = core.nms_dedupe(cands, nms_jaccard)

    out = []
    for k, (raw_r, gf_r, nisl) in enumerate(kept):
        out.append(dict(
            patch=k, n_raw=len(raw_r), n_res=len(gf_r), n_islands=nisl,
            epitope=antigen.fmt_epitope(gf_r), islands=antigen.fmt_islands(gf_r),
            keys=sorted(gf_r), vmd_selection=antigen.vmd_selection(gf_r)))
    return out
