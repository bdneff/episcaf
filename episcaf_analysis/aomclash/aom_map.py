"""Load the antibody-occupancy map and build the smoothed P1/P2 fields the scorer reads.

The map artifact (`aom_map.pkl`, extracted by `scripts/aom_extract_map.py` from John Altin's committed
AOMclash model, or rebuilt locally from SAbDab) holds one boolean 33^3 occupancy grid per SAbDab
complex, bit-packed. The deployable map is their mean = P(antibody occupies this epitope-frame point),
which `features.design_feats` reads at each scaffold atom to get `s_p1_free`, our shipped clash score.

Exactly matches John's `aom_score.mean_map` + `aom_features.target_context`:
  mean over complexes -> P1 = gaussian_filter(mean, 1.0), P2 = gaussian_filter(mean, 2.5).
No x/y flip is applied (his aom_build docstring mentions one, but the committed code does not).
"""
from __future__ import annotations

import pickle

import numpy as np

P1_SIGMA = 1.0    # ~2 A smoothing (STEP = 2 A/voxel)
P2_SIGMA = 2.5    # ~5 A smoothing


def load_map(path):
    """Load the occupancy-map artifact (dict with fields_packed, field_length, field_side, keys, pdbs)."""
    with open(path, "rb") as fh:
        return pickle.load(fh)


def mean_field(art, exclude_key=None, exclude_pdb=None):
    """Mean occupancy grid (field_side^3), i.e. P(antibody) per epitope-frame point (port of mean_map).

    `exclude_key` / `exclude_pdb` drop a target's own SAbDab complexes for a leave-one-antigen-out
    read (benchmark use); omit both for a genuinely new antigen. Returns (field, n_complexes_used).
    """
    grids = np.unpackbits(art["fields_packed"], axis=1)[:, :art["field_length"]]
    keep = np.ones(len(grids), bool)
    if exclude_key:
        keep &= art["keys"] != exclude_key
    if exclude_pdb:
        keep &= art["pdbs"] != str(exclude_pdb).lower()
    if keep.sum() == 0:
        raise ValueError("every complex in the map was excluded")
    n = art["field_side"]
    return grids[keep].mean(0).reshape(n, n, n).astype(np.float32), int(keep.sum())


def smoothed_maps(field):
    """(P1, P2) = the mean field smoothed at P1_SIGMA and P2_SIGMA voxels (what design_feats reads)."""
    from scipy.ndimage import gaussian_filter
    return gaussian_filter(field, P1_SIGMA), gaussian_filter(field, P2_SIGMA)


def occupancy_maps(path, exclude_key=None, exclude_pdb=None):
    """Convenience: load the artifact and return (P1, P2, n_complexes_used)."""
    art = load_map(path)
    field, n = mean_field(art, exclude_key, exclude_pdb)
    P1, P2 = smoothed_maps(field)
    return P1, P2, n
