"""AOMclash: antibody-free clash prediction via an antibody occupancy map.

Port of John Altin's `aomclash/` (AltinLab/epitope-universe). The occupancy map is antibody geometry
only (built from ~1,029 SAbDab complexes) and is AF-version agnostic.

What we ship (decided 2026-11-13, see results/aomclash/feature_set_decision.md): a SINGLE interpretable
feature, `s_p1_free` -- the occupancy map summed over the design's free scaffold atoms -- not John's
109-feature gradient-boosted model, which we keep only as a demoted head-to-head baseline.

Modules:
    aom_map   -- load the vendored occupancy map (aom_map.pkl) and build the smoothed P1/P2 fields the
                 scorer reads (mean over complexes + gaussian 1.0 / 2.5; leave-one-antigen-out aware).
    features  -- per-design geometry features in the epitope frame (port of aom_features.design_feats),
                 incl. the carved occupancy sum `s_p1_free`, plus the cylinder cross-check.
    model     -- the scorer: `fit_prior`/`score_prior` (isotonic calibration of s_p1_free onto the AF3
                 clash count); FEATURES_39 + the GBM are the kept-but-rejected baseline.

The map artifact is extracted from John's committed model by scripts/aom_extract_map.py (and is
rebuildable from SAbDab as a cross-check). Source reports: ~/Downloads/{AOMclash,growPatch,epiScaf}.pdf.
"""
from . import aom_map  # noqa: F401  (numpy/pickle only; scipy is imported lazily inside)
