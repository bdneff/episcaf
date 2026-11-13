# AOMclash explorations — what we tried before shipping one feature

Depth record for the AOMclash scorer. The **shipped** decision (single feature `s_p1_free`) lives one
level up in `results/aomclash/feature_set_decision.md`; this folder holds the alternatives we weighed
and the evidence, so the choice is auditable without being front and center.

All margins below are from **John Altin's AF2 campaign** (`aom_ablation.py`, `results/ablation.csv` in
his repo) and the 8-Nov AOMclash report, unless marked otherwise. Nothing here has been re-run on our
AF3 labels yet — the AF3 feature matrix is cluster-only (`$WS`) and not built.

## 1. The ladder of models we compared

| model | inputs | within-campaign ρ (post-AF) | verdict |
|---|---|---|---|
| map-only prior | **1** — `s_p1_free`, isotonic-calibrated | ~0.87–0.89 | **shipped** |
| 39-feature GBM | 39 (drop 70 histograms + 2 cylinder) | ~0.912 | rejected (+~0.02, overfit risk) |
| full bank GBM | 109–111 (John's committed model) | ~0.919 | rejected (no gain over 39) |
| native-aware cylinder | geometric surrogate (Brandon's prior method) | ~0.63 (report) | the thing AOMclash replaces |

The jump that matters is cylinder → occupancy map. Everything past the single carved-occupancy feature
is a fraction of a point for a large, correlated feature bank.

## 2. Ablation (why 39, and then why 1)

- Dropping the **70 position-histogram bins** (`h0..h34`, `hf0..hf34`) costs ~0.007 within-campaign ρ
  (0.919 → 0.912). They are ≈ redundant with the occupancy sums.
- The **2 cylinder counts** (`cyl_aware`, `cyl_plain`) add ~0.003.
- Dropping the **epitope-descriptor group** (6) was the single biggest hit in the 39-set
  (0.919 → 0.891) — but it is per-target context, identical for every design of one epitope, so it
  helps *across* targets, not for ranking designs *within* one target (the deployment metric).
- Collapsing all the way to the **single `s_p1_free`** feature lands at ~0.87–0.89. The ~0.02–0.05 the
  multi-feature models add is not worth the overfitting exposure on a set this size, and costs the
  interpretability and map-only reproducibility of the single feature.

## 3. The 39-feature set (the rejected alternative, kept as baseline)

Names match `features.design_feats` / `target_features`; defined in full in
`episcaf_analysis/aomclash/FEATURES.md`.

- **Occupancy (15):** `s_p1 s_p2 s_p1_free s_p1_near n_p05 n_p15 n_p30 n_p50 n_p15_free n_p30_free
  r_sum r_sum_free r_n30 r_n50 r_n70`
- **Outward-side counts (4):** `b_near20 b_out25 b_out15 b_out35_hi`
- **Antigen contact (4):** `ag3 ag4 ag6 ag_free_out25`
- **Scaffold shape (10):** `z_mean z_p90 z_max r_mean ca_rg n_scaf n_seg seg_max seg_sum apr`
- **Epitope descriptors (6):** `frame_ratio n_epi n_isl ag_out25 ag_out25_hi ag_near15`

Dropped from that model: 70 histograms + 2 cylinder counts. `tests/test_aomclash_model.py` asserts the
39 equal the extractor's full 111 key set minus exactly those 72, so the list can't silently drift
from `features.py`. The GBM uses John's exact hyperparameters (`model.py:GBM_KW`). `feats=None`
restores the full 111 bank.

## 4. Branch reconciliation (which version of John's code is current)

Checked all remote branches of `~/Desktop/projects/epitope-universe` (read-only, already fetched):

| branch | tip date | has `aomclash/` |
|---|---|---|
| `upstream/feat/aomclash` | 2026-10-25 | yes |
| `upstream/feat/clash-method` | **2026-11-04** (newest) | yes |
| `upstream/main` | 2026-09-14 | no (only `segmentation/`) |

- `aom_model.py` is **byte-identical** between `feat/aomclash` and `feat/clash-method`.
- Only `aom_features.py` changed between them (+67/−25): the **island-order fix** (`target_context`,
  `design_islands`). Feature keys are unchanged.
- The newest committed code (Nov 4) folds **leave-one-ANTIGEN-out** (`biological_antigen_key`) and
  excludes the cylinder counts → **109** model inputs (keeps the 70 histograms).
- The 8-Nov report's "leave-one-epitope-out / 111 incl. cylinder" is **not reflected in any committed
  branch** — the report is ahead of (or diverges from) the code. We followed the report for the
  rejected 39-feature baseline; moot for the shipped single-feature scorer.

## 5. Fold choice (for the rejected GBM baseline)

`leave_one_group_out` defaults to **leave-one-EPITOPE-out** (`model.py:epitope_groups`, prefix of a
`target` column), per the report. The committed code uses antigen folds; switching is a one-arg change
(pass an antigen-keyed `groups`). The single-feature scorer's isotonic fit is a global monotone map,
so it does not depend on this choice.

## 6. Open items to confirm cluster-side (when the AF3 matrix exists)

1. **Build the AF3 feature matrix** on `$WS` (needs the design structures). `model.py` consumes it; it
   does not create it. Staged job still to be written.
2. **Label column** is `af3_n_clash_res` (stage05) / `af3_clashes` (shipped libraries);
   `resolve_label` accepts either.
3. **Epitope group key** for the baseline folds: confirm the matrix carries a `target` column, or pass
   an explicit epitope/patch id as `groups`.
4. **Stage column** (`stage`) if the matrix pools post-RFD3 and post-AF, for the equal-weight scheme.
5. **Re-check the single feature on AF3 DP3/DP4** (§2 above) — confirmatory, not gating.
