# AOMclash: the shipped scorer (decision)

**Settled 2026-11-13.** Code: `episcaf_analysis/aomclash/model.py`
(`fit_prior` / `score_prior`), tested by `tests/test_aomclash_model.py`.

## Decision: ship one interpretable feature

The AOMclash scorer is a **single feature, `s_p1_free`** — the antibody-occupancy map read at every
scaffold atom, summed over only the atoms outside native-antigen space (`free`, > 4 Å from any native
antigen atom). It is isotonic-calibrated onto the AF3 clashing-residue count, so the output reads on
the same scale as `af3_n_clash_res` (the real clash label on the antibody arm). For ranking designs
within a target — the deployment use — the calibration is cosmetic, so the raw sum ranks identically.

Why this and not the full feature bank:

- **It is the mechanism.** Occupancy at the free atoms is directly "how much of the design sticks into
  antibody territory beyond what mimicking the native antigen explains" — the thing a clash *is*.
- **Nothing to overfit.** One monotone feature; no learned multi-feature regressor, no per-target fit.
- **Reproducible from (map + native antigen) alone** — no dependence on the assembled AF3 matrix to
  *ship* (we still check it; see below).
- **The extra features don't pay.** John's AF2 ablation and the 8-Nov report put this single feature
  at ~0.87–0.89 and the full model only ~0.02 higher. ~38 highly-correlated features for +0.02, on a
  set this size, is overfitting risk, not signal.

## To confirm (does not gate shipping)

The margins above are from John's **AF2** campaign. When the AF3 feature matrix exists (cluster-side,
`$WS`), re-check that the single feature reproduces its result on our **AF3 DP3/DP4** labels — we are
reasonably certain it will, since the approach is mechanism-based and already validated, but the check
is cheap and worth recording. Shipping is **not** gated on the 39-feature model.

## What else we tried

The 39-feature gradient-boosted alternative, the full 111-column bank, the ablation numbers, the
leave-one-epitope-vs-antigen fold question, the branch reconciliation of John's code, and the cylinder
head-to-head are all recorded in **`results/aomclash/explorations/`**. That code is kept in
`model.py` (`FEATURES_39`, `gbm`, `fit_final`, `leave_one_group_out`) as the baseline; `feats=None`
restores the full bank. It is there for depth, not as the shipped path.
