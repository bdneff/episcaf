# AOMclash: 39-feature model + branch reconciliation

Staged 2026-11-10. Code: `episcaf_analysis/aomclash/model.py`, tested by
`tests/test_aomclash_model.py`. Nothing trained yet — the feature matrix is cluster-only.

## 1. Branch reconciliation (which version of John's code is current)

Checked all remote branches of `~/Desktop/projects/epitope-universe` (read-only, already fetched):

| branch | tip date | has `aomclash/` |
|---|---|---|
| `upstream/feat/aomclash` | 2026-10-25 | yes |
| `upstream/feat/clash-method` | **2026-11-04** (newest) | yes |
| `upstream/main` | 2026-09-14 | no (only `segmentation/`) |

- `aom_model.py` is **byte-identical** between `feat/aomclash` and `feat/clash-method` (empty diff).
- Only `aom_features.py` changed between them (+67/−25): the **island-order fix** (`target_context`,
  `design_islands`). The feature **keys are unchanged** (same `feat` dict, same `h*/hf*`, same
  `target_feats`).

**So the newest committed code (Nov 4) still:**
- folds **leave-one-ANTIGEN-out** (`groups = biological_antigen_key`), and
- **excludes** the cylinder counts from the model (`NOT_FEATURES`) → **109** model inputs
  (111 computed − 2 cylinder; the 70 histograms are kept).

The 8-Nov report's "leave-one-epitope-out / 111 inputs incl. cylinder" is **not reflected in any
committed branch.** The report is ahead of (or diverges from) the code.

## 2. Our decision

- **Features: 39** (drop the 70 position histograms *and* the 2 cylinder counts). John's own ablation
  (`results/ablation.csv`) shows dropping the histograms costs ~0.007 post-AF within-campaign rho
  (0.919→0.912) and cylinder adds ~0.003 — a safe plateau. [UNVERIFIED on AF3: the ablation is John's
  AF2 data; re-confirm the plateau holds once we have the AF3 matrix.]
- **Folds: leave-one-EPITOPE-out** (default), per the report and our project decision — *not* the
  antigen folds the committed code uses. Flagged in `model.py`; switching to antigen folds is a
  one-arg change (pass an antigen-keyed `groups`), and the full 111 bank is `feats=None`.

## 3. The 39 features (names match `features.design_feats` / `target_features`)

- **Occupancy (15):** `s_p1 s_p2 s_p1_free s_p1_near n_p05 n_p15 n_p30 n_p50 n_p15_free n_p30_free
  r_sum r_sum_free r_n30 r_n50 r_n70`
- **Outward-side counts (4):** `b_near20 b_out25 b_out15 b_out35_hi`
- **Antigen contact (4):** `ag3 ag4 ag6 ag_free_out25`
- **Scaffold shape (10):** `z_mean z_p90 z_max r_mean ca_rg n_scaf n_seg seg_max seg_sum apr`
- **Epitope descriptors (6):** `frame_ratio n_epi n_isl ag_out25 ag_out25_hi ag_near15`

Dropped from the model: 70 histograms (`h0..h34`, `hf0..hf34`) + 2 cylinder (`cyl_aware`,
`cyl_plain`). `tests/test_aomclash_model.py` asserts the 39 equal the extractor's full 111 key set
minus exactly those 72, so the list can't silently drift from `features.py`.

## 4. What's staged

- `episcaf_analysis/aomclash/model.py` — `FEATURES_39`, `gbm()` (John's exact hyperparameters),
  `stage_weights` (equal per-stage weight), `leave_one_group_out` (+ isotonic `aom_prior_sum`
  baseline), `fit_final`, `predict`, `per_group_spearman`, `epitope_groups`, and a `--frame` CLI.
  Operates on a prepared feature DataFrame; no structures, no map build.
- `tests/test_aomclash_model.py` — 8 tests, all pass locally (`/usr/bin/python3`). Existing
  `tests/test_aomclash_features.py` still passes (7/7).

## 5. Open questions (confirm cluster-side before training)

1. **Build the AF3 feature matrix** on `$WS` (needs the design structures) — this module consumes it;
   it does not create it. [staged job still to be written]
2. **Label column** is `af3_n_clash_res` (stage05) / `af3_clashes` (shipped libraries); `resolve_label`
   accepts either. Confirm which the assembled matrix carries.
3. **Epitope group key**: `epitope_groups()` derives it as the prefix of a `target` column
   (`target.split("_")[0]`, matching `dp4_binder_ml.py`). Confirm the matrix has that column, or pass
   an explicit epitope/patch id as `groups`.
4. **Stage column**: if the matrix pools post-RFD3 (stage 0) and post-AF (stage 1), name it `stage`
   for the equal-weight scheme; otherwise weighting is uniform.
5. Re-confirm the 39-feature plateau and the fold choice (epitope vs antigen) on the AF3 data, since
   both were decided from John's AF2 ablation and report.
