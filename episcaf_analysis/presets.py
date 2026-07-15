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

# ---- CANDIDATE (John 2026-07-15): saturate RMSD/PAE, penalize clash absolutely ----------------
# John flagged mAb designs that are excellent on RMSD/PAE but poor on AF3 clashes, and asked whether
# the scorer over-optimizes RMSD at the expense of clashes that likely kill binding. Percentile is the
# culprit: it is scale-blind, so it over-resolves the dense sub-1A RMSD pile and lets meaningless RMSD
# gains outrank real clash reductions. Fix = saturating transforms:
#   - RMSD/PAE -> sigmoid with midpoint at the four-filter threshold (epitope RMSD 1, overall 2, PAE 5).
#     Below the threshold the score is ~flat, so 0.3A and 0.9A are treated the same (John's "cap"); the
#     term stops paying for over-optimization and leaves room for accessibility to break ties.
#   - af3_n_clash_res -> sigmoid ABSOLUTE (not percentile), midpoint 6 (the clash distribution's mass
#     sits 0-10, median 4), k=0.5 so it discriminates across 1..~15 and heavily penalizes >12. Absolute
#     so "this whole cohort is too clashy" is expressible. Accessibility is the one metric that tracks
#     experimental binding (sec:whatpredicts), so penalizing it hard is the point.
#
# REALITY CHECK on the actual 103 C1 pool (metrics_whole_epitope_103.csv, 140,716; scripts/scoring_worlds.py):
# the effect is SMALL. Pooled over C1 targets, mixed moves clash median only 2 -> 1 vs percentile, and for
# a hard target (6cyf) ALL scorings give the same top-10 (clash ~7) -- its whole-epitope pool simply has no
# lower-clash designs, so that is a GENERATION limit, not a scoring one. Earlier alarming results (a naive-
# sigmoid "backfire" to clash 11, and a misfold blowup to overall_rmsd 9.4) were ARTIFACTS of the stale
# 104-mer pool and DO NOT reproduce on real 103 data. So clash 0.50 here is a mild, safe tilt toward
# accessibility, not a dramatic fix. The real test is C2 (single-island, cluster metrics_dual_island.parquet)
# -- that is where the low-clash 6cyf designs John cited live. ALL PROVISIONAL: re-fit weights + k on DP4
# binding data (that is what C5 is for).
ANTIBODY_SIGMOID = dict(
    gate=None,
    scope="pooled",
    antigen_col="antigen",
    select=dict(group="id", topk=5),
    metrics={
        "af3_n_clash_res":    dict(weight=0.50, better="low", transform="sigmoid", midpoint=6.0, k=0.5),
        "epitope_chunk_rmsd": dict(weight=0.25, better="low", transform="sigmoid", midpoint=1.0, k=4.0),
        "overall_rmsd":       dict(weight=0.15, better="low", transform="sigmoid", midpoint=2.0, k=2.0),
        "epitope_pae":        dict(weight=0.10, better="low", transform="sigmoid", midpoint=5.0, k=1.0),
    },
)

PRESETS = {"twelvemer": TWELVEMER, "antibody": ANTIBODY, "antibody_sigmoid": ANTIBODY_SIGMOID}
