# DP4 IM0276 — binder-prediction ML write-up

`dp4_ml_findings.tex` (build with `tectonic dp4_ml_findings.tex`) is the findings report for the DP4
IM0276 assay: which features predict binding, which model form fits, and how a learned re-ranker compares
to the composite and Lawson's four-filter.

The figures in `figs/` are a committed snapshot. They are produced by the named scripts in `scripts/`
(see the reproducibility table at the end of the report), which write into `data/dp4_binding/figs/` from
the gitignored assay data. To refresh:

```
# rerun the scripts you changed (each writes into data/dp4_binding/figs/), then:
cp data/dp4_binding/figs/{binder_ml_v2,additive_shapes,headtohead,target_compare,context_rules,model_vs_composite,roc_comparison}.png results/dp4_ml/figs/
cd results/dp4_ml && tectonic dp4_ml_findings.tex
```

The assay data itself lives under `data/dp4_binding/` and is gitignored (see the repo `.gitignore`).
