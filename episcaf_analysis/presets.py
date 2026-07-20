"""Scorer dials — one preset per dataset. Edit these to tune; no CLI needed.

A preset is a dict:
  gate        : (column, threshold)      # keep rows where column <= threshold (lower-better gate)
  scope       : "pooled" | "per_antigen" # population used for the per-metric transforms
  antigen_col : column name used when scope == "per_antigen"
  select      : dict(group=<col>, topk=<int>)   # keep top-k by composite within each group
  metrics     : { column: dict(weight, better="low"|"high", transform, **params) }

transform in {percentile, minmax, zscore, sigmoid, identity}
  - sigmoid takes extra params: midpoint, k   ->  1 / (1 + exp(s*k*(x-midpoint))),
    with s = +1 when better="low" else -1 (a saturating activation around `midpoint`).

This is a single-layer perceptron with a per-input nonlinearity and hand-set
weights — "feed-forward, no backprop". To LEARN the weights instead, fit a logistic
regression of the transformed features against the DP3 `is_pass` label
(see docs/REORG.md, the v2 hook).
"""

# Weights are a PRIOR set from the DP3 binding data (manuscript sec:whatpredicts): correlating
# each metric with experimental cognate enrichment WITHIN antibody (the deconfounded view),
# ACCESSIBILITY is the strongest predictor and epitope RMSD next (-0.14); overall RMSD (-0.08)
# and PAE (~0, mean_pae tested) carry little binding signal. So accessibility + epitope RMSD get
# the weight (0.35 each); the other two are kept small. A hand-set prior from an all-passing set,
# to be re-fit once DP4 spans the full metric space.
#
# NO GATE: we do not hard-filter before ranking. The composite already penalizes bad metrics
# softly, so a hard threshold would only discard designs the ranking buries anyway (and risk
# dropping a design that is weak on one axis but excellent overall). Selection = rank all designs
# on the composite, keep the top-k per group.
#
# The accessibility term depends on what we know: with a known antibody (ANTIBODY) we use the
# REAL clash af3_n_clash_res directly; with no antibody (TWELVEMER) we use the cylinder surrogate
# the clash is approximated by (the cylinder exists precisely to stand in where there is no
# antibody). Same accessibility weight either way.
TWELVEMER = dict(   # no-antibody set — accessibility via the native-aware cylinder surrogate
    gate=None,
    scope="per_antigen",                 # rank within each antigen, not pooled across all
    antigen_col="antigen",
    select=dict(group="id", topk=5),     # top-5 per epitope  (verify `id` is the epitope key)
    metrics={
        "cylinder_native_aware":  dict(weight=0.35, better="low", transform="percentile"),
        "epitope_chunk_rmsd":     dict(weight=0.35, better="low", transform="percentile"),
        "overall_rmsd":           dict(weight=0.15, better="low", transform="percentile"),
        "epitope_pae":            dict(weight=0.15, better="low", transform="percentile"),
    },
)

ANTIBODY = dict(   # DP3 / mAb set — known antibody, so accessibility via the REAL clash
    gate=None,
    scope="pooled",
    antigen_col="antigen",
    select=dict(group="id", topk=5),     # per-island top-5 deliverable (group per island in the run)
    metrics={
        "af3_n_clash_res":      dict(weight=0.35, better="low", transform="percentile"),
        "epitope_chunk_rmsd":   dict(weight=0.35, better="low", transform="percentile"),
        "overall_rmsd":         dict(weight=0.15, better="low", transform="percentile"),
        "epitope_pae":          dict(weight=0.15, better="low", transform="percentile"),
    },
)

# ---- ANTIBODY_SOFTGATE: the adopted known-antibody scorer (soft gate + global-pass promotion) -----
# Percentile (ANTIBODY, above) is scale-blind: it over-resolves the dense sub-1A RMSD pile and lets
# meaningless RMSD gains outrank real clash reductions, even though accessibility is the one metric that
# tracks experimental binding (sec:whatpredicts). Soft-gate fixes this with saturating activations and a
# steep-but-finite gate on fold quality:
#   - af3_n_clash_res: sigmoid, midpoint 6, k=0.5, weight .45 -- BROAD, so it RANKS accessibility across
#     the 1..~15 band (the mass) and heavily penalizes >12. Weighted heavily (it predicts binding).
#   - epitope_chunk_rmsd / overall_rmsd: STEEP sigmoids (k=4) at the four-filter thresholds (1, 2).
#     Steep k acts as a SOFT GATE on fold quality -- misfolds are crushed toward 0 but never TO 0, so no
#     island is ever dropped (unlike a hard fold floor, which empties 3/87 C2 islands). As k -> inf this
#     is Lawson's hard filter; finite k keeps it soft and coverage-safe.
#   - epitope_pae: gentler sigmoid (k=1.2), midpoint 2.5 -- a rank nudge toward a rigid epitope, not a
#     gate. Midpoint set from the data (see the provenance note below), NOT the global mean_pae<5.
# On C2 (single-island, cluster) this cuts clash 6->2 pooled (6cyf 14.5->3) while keeping the fold and
# all 87 islands. On C1 it helps mildly and never hurts (that arm is generation-limited). Chosen over
# both percentile (leaves clashes on the table) and hard-weight-only mixes (backfire / misfold blowup).
#
# pass_bonus = John's rule (2026-07-16): every design clearing ALL four Lawson filters ranks above every
# design that does not. Done as a soft AND -- P = product of steep sigmoids at the thresholds, ~1 iff all
# pass -- scaled by a gain > the composite range (weights sum to 1, scores in [0,1], so gain 2 > 1). See
# score.py. Uses the four-filter's OWN metrics (global mean_pae<5, clash==0). Promotes, never excludes:
# an epitope with no passer just ranks on the composite. On C1: 725 soft-passers vs 727 hard (0.52%).
#
# Midpoints: the fold gates (epitope_chunk_rmsd 1, overall_rmsd 2) and the clash centre (6) sit at the
# DP3 thresholds / the clash mass. epitope_pae is the EXCEPTION -- its midpoint is 2.5, set from the
# data, NOT the global mean_pae<5 threshold. Provenance (measured 2026-07-20 on metrics_whole_epitope_103
# .csv, 140,716 designs): epitope_pae is the intra-epitope PAE BLOCK, made of short-range pairs, so it
# runs far below the whole-matrix mean_pae -- four-filter passers median 1.85 A, pool median 9.30. It
# correlates with epitope_chunk_rmsd (r=0.80), so that 1.85 is partly a conditioning artifact; de-
# confounded (pass the other three filters, NO epi-rmsd cut) it is still 1.98, and designs failing only
# epi-rmsd sit at 3.57. So 2.5 is the half-credit point between a good epitope (~2) and a marginal one
# (~3.6); the old 5.0 (borrowed from the global threshold) sat in the tail and barely discriminated.
# k / weights / gain remain provisional dials to be re-fit on DP4 binding data (that is what C5 is for).
ANTIBODY_SOFTGATE = dict(
    gate=None,
    scope="pooled",
    antigen_col="antigen",
    select=dict(group="id", topk=5),
    metrics={
        "af3_n_clash_res":    dict(weight=0.45, better="low", transform="sigmoid", midpoint=6.0, k=0.5),
        "epitope_chunk_rmsd": dict(weight=0.25, better="low", transform="sigmoid", midpoint=1.0, k=4.0),
        "overall_rmsd":       dict(weight=0.20, better="low", transform="sigmoid", midpoint=2.0, k=4.0),
        "epitope_pae":        dict(weight=0.10, better="low", transform="sigmoid", midpoint=2.5, k=1.2),
    },
    pass_bonus=dict(
        gain=2.0, k=12.0,
        criteria={"epitope_chunk_rmsd": 1.0, "overall_rmsd": 2.0, "mean_pae": 5.0, "af3_n_clash_res": 0.5},
    ),
)

PRESETS = {"twelvemer": TWELVEMER, "antibody": ANTIBODY, "antibody_softgate": ANTIBODY_SOFTGATE}
