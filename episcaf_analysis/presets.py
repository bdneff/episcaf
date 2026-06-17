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

TWELVEMER = dict(
    gate=("epitope_chunk_rmsd", 2.5),
    scope="per_antigen",                 # rank within each antigen, not pooled across all
    antigen_col="antigen",
    select=dict(group="id", topk=5),     # top-5 per epitope  (verify `id` is the epitope key)
    metrics={
        "epitope_chunk_rmsd":     dict(weight=0.35, better="low", transform="percentile"),
        "epitope_pae":            dict(weight=0.25, better="low", transform="percentile"),
        "overall_rmsd":           dict(weight=0.15, better="low", transform="percentile"),
        "cylinder_native_aware":  dict(weight=0.25, better="low", transform="percentile"),
    },
)

ANTIBODY = dict(   # DP3 / mAb set — has real af3_n_clash_res ground truth + is_pass
    gate=("epitope_chunk_rmsd", 2.5),
    scope="pooled",
    antigen_col="antigen",
    select=dict(group="id", topk=15),
    metrics={
        "epitope_chunk_rmsd":   dict(weight=0.35, better="low", transform="percentile"),
        "epitope_pae":          dict(weight=0.25, better="low", transform="percentile"),
        "overall_rmsd":         dict(weight=0.15, better="low", transform="percentile"),
        "cylinder_ca_clashes":  dict(weight=0.25, better="low", transform="percentile"),
    },
)

# ---- example of the sigmoid squish, if you want absolute (population-independent) scores ----
# "overall_rmsd": dict(weight=0.15, better="low", transform="sigmoid", midpoint=3.0, k=1.5),

PRESETS = {"twelvemer": TWELVEMER, "antibody": ANTIBODY}
