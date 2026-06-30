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

# ---- example of the sigmoid squish, if you want absolute (population-independent) scores ----
# "overall_rmsd": dict(weight=0.15, better="low", transform="sigmoid", midpoint=3.0, k=1.5),

PRESETS = {"twelvemer": TWELVEMER, "antibody": ANTIBODY}
