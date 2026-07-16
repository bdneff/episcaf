# MIGRATION — what moved where (and what to double-check)

> **Historical snapshot** of the one-time code reorg (old → new file map). For the current pipeline
> see `docs/PIPELINE.md`; the `archive/` READMEs point at what replaced each superseded file.

This repo was assembled from the `episcaf_v2_bneff` workspace code (the `code_clean`
tarball: `episcaf_pipeline/`, `scripts/`, `legacy_*`, and top-level `.py`). Data was
**not** touched — it stays on `/tgen_labs`, referenced via `configs/paths.py`.

Guiding rule: **nothing was deleted.** Anything superseded, duplicated, or
exploratory went to `archive/` (tracked in git). Re-promoting something is one
`git mv` — so a wrong canonical/archive call is cheap to fix.

## Kept canonical

| from                                   | to                                    | why |
|----------------------------------------|---------------------------------------|-----|
| `episcaf_pipeline/` (whole package)    | `episcaf_pipeline/` (unchanged)       | already a clean package; has the stage05 coupling below |
| `legacy_steps/05_rmsd_vs_af3.py`       | `legacy_steps/05_rmsd_vs_af3.py`      | **live dependency** — `episcaf_pipeline/cli.py stage05` shells out to it at `<repo_root>/legacy_steps/05_rmsd_vs_af3.py` |
| `compute_metrics.py` (858 lines)       | `episcaf_analysis/compute_metrics.py` | README step-5 canonical (DP3 / no-MPNN 4-filter builder) |
| `scripts/native_cylinder_core.py`      | `episcaf_analysis/native_cylinder_core.py` | cylinder geometry core (self-tests pass) |
| `scripts/build_12mer_metrics.py`       | `episcaf_analysis/build_12mer_metrics.py` | 12-mer metrics builder; `import native_cylinder_core as C` (kept co-located) |
| `scripts/false_positive_check.py`      | `episcaf_analysis/false_positive_check.py` | recent validation tool |
| `scripts/composite_swap_validate.py`   | `episcaf_analysis/composite_swap_validate.py` | recent validation tool |
| `scripts/plot_*` (fp_reduction, composite_distribution, cylinder_before_after, filters_vs_composite) | `episcaf_analysis/viz/` | keeper plots |
| `scripts/01_*`,`02_submit_mpnn.py`,`03_emit_af3_jsons_mpnn.py`,`04_af3_array.sbatch` | `scripts/` | 12-mer/MPNN pipeline-branch steps (README pipeline) |
| `requirements.txt`                     | `requirements.txt` + `pyproject.toml` | deps; scipy/matplotlib added (used by core/viz) |
| `README_v2`                            | `docs/README_v2_original.md`          | original experiment writeup + 4-filter defs |

## New code (written for this repo)

- `episcaf_analysis/score.py` — the one config-driven composite scorer (gate →
  per-metric transform within scope → weighted sum → top-k per group).
- `episcaf_analysis/presets.py` — the scoring dials (`twelvemer`, `antibody`).
- `configs/paths.py` — absolute `/tgen_labs` data locations.
- `tests/test_scoring.py` — transform + selection unit tests (no data needed).
- `pyproject.toml`, `.gitignore`, `README.md`, `docs/REORG.md`.

These supersede the scoring/ranking sprawl: `rank_composite_12mer.py` (not present
in this snapshot), `apply_composite_filter.py`, `apply_filters.py`, `score_filters*.py`,
`filter_passes*.py`, and the `add_*`/`scan_*` cylinder one-offs. Their logic is
absorbed into the preset/transform model; the originals remain in `archive/`.

## Archived (preserved in `archive/`, not on the critical path)

- **dead twins:** `Xcompute_metrics.py`, `X02_submit_mpnn.py`, `X_score_af3_filters.py`.
- **backups:** every `*.bak` / `*.bak_YYYYMMDD_*`.
- **Lawson-replication variants:** `compute_rfd3_af3_metrics_lawson.py` (+baks),
  `compute_af3_clashes_like_lawson.py`, `lawson_rmsd.py`.
- **parameter sweeps / one-offs:** `scan_cylinder_params.py`, `scan_native_cylinder.py`,
  `scan_weighted_cylinder.py`, `add_clash_metrics*.py`, `add_cone_clash_metric.py`,
  `add_cylinder_metric.py`, `add_fab_probe_metric.py`, `add_mean_pae_af3.py`.
- **superseded filters/scorers:** `apply_composite_filter.py`, `apply_filters.py`,
  `filter_passes.py`, `filter_passes_dp2.py`, `score_filters.py`,
  `score_filters_with_epitope_rmsd.py`.
- **misc/legacy:** `legacy_steps/{02,03,04,old_03}`, `legacy_sbatch/`, `legacy_tools/`,
  `print_results.py`, `sbatch_metrics_my_run.sh`, `af3_cif_to_pdb.py`, `new_run.sh`,
  old `README.md`.

## ⚠ Verify against the data before trusting

These are **judgment calls** made without running against the real CSVs. Confirm,
then `git mv` out of `archive/` if any is actually canonical:

1. **`compute_metrics_mpnn.py` (22 K)** — the README mentions “`compute_metrics.py`
   *adapted for the MPNN pipeline*.” This may be the live DP3+MPNN metrics builder
   rather than archive. Check which one produced `runs/run_rfd3_mpnn/04_filter/*.csv`.
2. **`build_metrics_all.py`, `compute_decomposed_metrics.py`,
   `compute_rfd3_af3_metrics.py`, `compute_rfd3_af3_metrics_clean.py`,
   `compute_af3_clashes_gemmi.py`, `compute_cylinder_clash.py`** — overlapping
   metric builders; confirm none is the actual producer of a CSV you still use.
3. **`compute_sasa_retention.py`, `epitope_secondary_structure.py`** — standalone
   analyses you may still want; promote to `episcaf_analysis/` if so.
4. **`validate_dp2_rmsd_full.py` + `submit_validate_dp2.sh`** — dp2 ground-truth
   validation; useful but not core scoring.
5. **`presets.py` column names** — `cylinder_native_aware` (12mer) and
   `cylinder_ca_clashes` (antibody) must match the metrics CSVs; the scorer warns
   and renormalizes if a metric column is missing. Confirm `select.group="id"` is
   the right per-epitope key for the 12-mer set.
6. **MDAnalysis upgrade** — if your local `build_12mer_metrics.py` was the
   MDAnalysis-superposition version (validated identical to the hand-rolled Kabsch),
   make sure this snapshot's copy matches; re-run only for clean provenance.

## Optional phase 2 (do with Claude Code, against data, so it's tested)

- Unify under `src/episcaf/` with `pipeline/` as a subpackage — requires vendoring
  `legacy_steps/05_rmsd_vs_af3.py` into the package and updating `cli.py stage05`'s
  path logic (currently `<repo_root>/legacy_steps/05_rmsd_vs_af3.py`).
- Fold the validated metric builders into `episcaf_analysis/metrics/` as real modules
  (rather than placed scripts).
- Fit scorer weights via logistic regression on DP3 `is_pass` (the “backprop”).
