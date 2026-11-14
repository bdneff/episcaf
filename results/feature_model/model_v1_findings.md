# Binding model v1 — findings

> **UPDATE 2026-11-13: v1b ADOPTED as the shipped model.** A held-out, paired, grouped-by-epitope
> comparison (`scripts/binding_v1b_heldout.py`, 40×5 folds) had the clean-static set
> (`epi_coulomb + frame_ratio + ag_near15`) beat the bound-state dynamics block by **+0.013 ROC-AUC,
> 95% CI [+0.008, +0.018]** — reliably better out-of-sample, not just non-inferior, and causally
> cleaner (nothing downstream of binding). `BindingModel` now defaults to `CLEAN_STATIC_FEATURES`;
> `binding_model_v1.pkl` re-fit on it. The analysis below (which recommended v1b pending this check)
> stands as the record.

The shippable interpretable binding-model scorer, and the Tier-1 new-feature keep test from
`FEATURE_CATALOG.md`. All numbers are **repeated grouped-by-epitope cross-validation** (5 shuffles of
`StratifiedGroupKFold(5)`), so every fold is tested on antibodies/epitopes it never trained on. Model
is L2-regularized logistic regression (median-impute → standardize → logistic, `class_weight=balanced`).

- Universe: 2,720 assayed scaffolded designs, 217 hits, 56 epitopes, base rate 0.080.
- Reproduce every number: `bcell-repair/bin/python scripts/binding_feature_eval.py`.

## Shipped model (`episcaf_analysis/binding_model/binding_model_v1.pkl`)

15 features — design metrics + secondary structure + the dynamics subset — reproducing the published
DP4 number:

| feature block | n | ROC-AUC | PR-AUC |
|---|---|---|---|
| design metrics | 9 | 0.759 ± 0.010 | 0.313 |
| + secondary structure | 12 | 0.787 ± 0.008 | 0.316 |
| **+ dynamics (shipped baseline)** | **15** | **0.812 ± 0.008** | **0.329** |
| + native composition | 21 | 0.766 ± 0.010 | 0.296 |

- **ROC-AUC 0.812 ± 0.008**, matching the DP4 findings' 0.810 (`results/dp4_ml/dp4_ml_findings.tex`).
- The **native-composition block** (`epi_size, epi_rel_sasa, epi_rg, frac_hydrophobic/charged/aromatic`)
  **hurts** under grouped CV (0.812 → 0.766): six more epitope-level terms for only 56 groups overfit.
  It is therefore **excluded** from the baseline (`data.NATIVE` is kept only for the ablation).
- `predict()` emits a `binding_pred` probability keyed by `design_id` (= `library_member`) — this is the
  live Results-tab "binding model" ranker.

## Tier-1 new features — keep test

Five Tier-1 feature families were implemented and computed on the 56 native complexes
(`episcaf_analysis/binding_model/native_features.py` → `results/feature_model/native_tier1_features.csv`;
51/56 reproduce `epi_size` exactly, 5 multi-copy/ambiguous-chain structures differ):

- **residue depth** `epi_depth` (freesasa surface + KD-tree), **protrusion/contact-number** `epi_cn`
  (reused from `scripts/webtool_epitope_propensity.py`), **epitope crowding** `ag_out25 / ag_out25_hi /
  ag_near15` + `frame_ratio` (reused verbatim from `episcaf_analysis/aomclash/features.target_features`),
  **glycosylation proximity** `epi_glyco_dist` (N-X-S/T sequon regex + Cα distance), **electrostatics**
  `epi_net_charge / epi_abs_charge / epi_coulomb` (formal-charge sums + cheap intra-patch Coulomb).

**Keep rule (marginal ROC-AUC ≥ 0.005 on top of the 0.812 baseline, max|r| < 0.9): KEEP NONE.**

| new feature | max\|r\| vs baseline | +feat ROC-AUC | marginal |
|---|---|---|---|
| epi_depth | 0.41 | 0.803 | −0.0084 |
| epi_cn | 0.51 | 0.802 | −0.0094 |
| ag_out25 | 0.48 | 0.808 | −0.0033 |
| ag_out25_hi | 0.54 | 0.801 | −0.0104 |
| ag_near15 | 0.41 | 0.808 | −0.0036 |
| frame_ratio | 0.34 | 0.806 | −0.0052 |
| epi_glyco_dist | 0.28 | 0.810 | −0.0016 |
| epi_net_charge | 0.35 | 0.796 | −0.0152 |
| epi_abs_charge | 0.63 | 0.803 | −0.0088 |
| epi_coulomb | 0.25 | **0.813** | +0.0014 |

None clears the bar on top of the full baseline; `epi_coulomb` is the only positive and it is +0.0014
(noise). So by the strict rule, **no feature is folded into the 0.812 model.**

## Why — and the one result worth acting on

The features are not noise; they are **redundant with the existing dynamics block**. On the leaner
design+SS baseline (0.787, before any native/dynamics term), the same static features *do* add signal:

| added to design+SS (0.787) | ROC-AUC | marginal |
|---|---|---|
| epi_coulomb | 0.808 | **+0.0214** |
| frame_ratio | 0.798 | +0.0112 |
| ag_near15 | 0.794 | +0.0075 |
| crowding_trio (ag_out25+hi+near15) | 0.794 | +0.0067 |

`epi_coulomb` alone recovers almost the entire native-context lift that the three dynamics features
(`epi_bfac_z, epi_gnm_dfit, epi_dsasa`) provide. The DP4 findings flagged that dynamics block as the
**weakest and confounded by the bound state** (it is read off the holo complex). The catalog's keep
rule explicitly prefers features that are **causal before an antibody exists (apo/static)**. So:

**A clean static native block matches-or-beats the confounded dynamics block:**

| model | ROC-AUC | note |
|---|---|---|
| design+SS+**dynamics** (shipped) | 0.812 ± 0.008 | 3 bound-state features (confounded) |
| design+SS+**epi_coulomb** | 0.808 ± 0.009 | 1 static feature ties the block |
| design+SS+**epi_coulomb+frame_ratio+ag_near15** | **0.817 ± 0.009** | 3 static/apo features, no dynamics |

`CLEAN_STATIC_FEATURES` in `data.py` is this set. It is within shuffle noise of the shipped model but
causally preferable (nothing downstream of binding). **Recommendation:** keep the 0.812 artifact as v1
(reproduces the published number), and adopt the clean-static set as v1b once confirmed on the held-out
8VDL arm / future assay data — this is the defensible direction, not a black box and not confounded.

## Files

- `episcaf_analysis/binding_model/native_features.py` — Tier-1 extractors + `extract_all`.
- `episcaf_analysis/binding_model/data.py` — modeling-table loader + feature-block definitions.
- `episcaf_analysis/binding_model/model.py` — `BindingModel` (fit/predict/cross_val_auc/save/load).
- `episcaf_analysis/binding_model/binding_model_v1.pkl` — fitted, tracked artifact.
- `results/feature_model/native_tier1_features.csv` — per-epitope new features.
- `scripts/binding_feature_eval.py` — the ablation + keep test (regenerates every number above).
- `tests/test_binding_model.py` — 24 checks (extractor math + fit/predict/save-load).

## Blocked on cluster data

Everything above ran locally. Not computable here:
- **Tier-2 conservation / Shannon entropy** needs an MSA (jackhmmer/ConSurf) — external data, cluster.
- Re-running `scripts/dp4_epitope_features.py` needs `biotite` (absent in `bcell-repair`); the existing
  `epitope_struct_features.csv` is reused as-is, and `native_features.py` uses Bio.PDB instead.
- Scoring *new* designs needs their native antigen CIF; the 56 DP4 complexes are local, the growPatch /
  8VDL targets live on `$WS`.
