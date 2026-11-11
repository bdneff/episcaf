"""growPatch: no-antibody epitope-candidate enumerator (port of John Altin's segmentation/growPatch).

`core` is the pure seed-and-grow algorithm (numpy only); `read` (added next) is the CIF + freesasa
adapter that builds the residue list and the relSASA surface set the core consumes. See
PARAMS.md for the fitted defaults and the method write-up.
"""
from . import core
from .core import (
    SASA_CUT, PATCH_SIZE, GAP, MIN_RES, MAX_ISLANDS, MIN_ISLAND_LEN,
    CONTACT, GROW_LINK, NMS_JACCARD,
    build_adjacency, grow_patch, islands, gap_fill, trim_short_islands, jaccard,
    enumerate_raw_patches, nms_dedupe, candidate_patches, greedy_set_cover,
)

__all__ = [
    "core",
    "SASA_CUT", "PATCH_SIZE", "GAP", "MIN_RES", "MAX_ISLANDS", "MIN_ISLAND_LEN",
    "CONTACT", "GROW_LINK", "NMS_JACCARD",
    "build_adjacency", "grow_patch", "islands", "gap_fill", "trim_short_islands", "jaccard",
    "enumerate_raw_patches", "nms_dedupe", "candidate_patches", "greedy_set_cover",
]
