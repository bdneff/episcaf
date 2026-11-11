"""AOMclash: antibody-free clash prediction via an antibody occupancy map.

Port of John Altin's `aomclash/` (AltinLab/epitope-universe @ feat/aomclash, 27ed71d) into
episcaf_v2, retrained on AlphaFold3 clash labels (John's reference trains on AlphaFold2).

The occupancy map is antibody geometry only (built from SAbDab complexes) and is AF-version
agnostic; only the regressor's labels and the design structures differ between AF2 and AF3.

Modules:
    features  -- per-design geometry features in the epitope frame (port of aom_features.design_feats),
                 plus the cylinder cross-check (port of aom_features.cylinder_counts).

See docs: the three source reports live in ~/Downloads/{AOMclash,growPatch,epiScaf}.pdf; the map
builder is John's aom_build.py (reused unchanged apart from a local SAbDab I/O adapter).
"""
