"""growPatch: no-antibody epitope-candidate enumerator (port of John Altin's segmentation/growPatch).

`core` is the pure seed-and-grow algorithm (numpy only); `read` is the structure + freesasa adapter
that builds the residue list and the relSASA surface set the core consumes, and exposes the
production entry `enumerate_patches(cif, chains, ...) -> candidate epitopes`. See PARAMS.md for the
fitted defaults and the method write-up.
"""
from . import core
from .core import (
    SASA_CUT, PATCH_SIZE, GAP, MIN_RES, MAX_ISLANDS, MIN_ISLAND_LEN,
    CONTACT, GROW_LINK, NMS_JACCARD,
    build_adjacency, grow_patch, islands, gap_fill, trim_short_islands, jaccard,
    enumerate_raw_patches, nms_dedupe, candidate_patches, greedy_set_cover,
)

# read.py pulls in Bio.PDB + freesasa, so import it lazily: structure-free users of `core` (and the
# synthetic tests) still work without those deps installed.
try:
    from . import read, eyeview
    from .read import (Antigen, read_antigen, enumerate_patches, enumerate_candidates,
                       compute_relsasa, MAXASA)
    from .eyeview import enumerate_eyeview
    _HAS_READ = True
except ImportError:  # Bio.PDB / freesasa not available
    _HAS_READ = False

__all__ = [
    "core",
    "SASA_CUT", "PATCH_SIZE", "GAP", "MIN_RES", "MAX_ISLANDS", "MIN_ISLAND_LEN",
    "CONTACT", "GROW_LINK", "NMS_JACCARD",
    "build_adjacency", "grow_patch", "islands", "gap_fill", "trim_short_islands", "jaccard",
    "enumerate_raw_patches", "nms_dedupe", "candidate_patches", "greedy_set_cover",
]
if _HAS_READ:
    __all__ += ["read", "eyeview", "Antigen", "read_antigen", "enumerate_patches",
                "enumerate_candidates", "enumerate_eyeview", "compute_relsasa", "MAXASA"]
