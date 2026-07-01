# Cylinder surrogate — geometry parameters (provenance and sweep)

The native-aware cylinder (the no-antibody accessibility surrogate, manuscript §`sec:noab`) has
four geometry dials. This file is the single record of what they are, where they came from, and
how they are chosen. It exists because they were previously set as bare constants with **no
recorded justification or sweep** — this fixes that.

## The parameters

| dial | symbol | current value | defined in |
|------|--------|---------------|------------|
| base offset | `OFFSET` | **−4.0 Å** (base = epitope centroid + OFFSET·normal) | `episcaf_analysis/native_cylinder_core.py:26`, `scripts/dp3_native_cylinder.py:37` |
| radius | `RADIUS` | **16.0 Å** | same |
| height | `HEIGHT` | **40.0 Å** | same |
| native carve distance | `exclude_dist` | **1.0 Å** (production) | passed per-run (`--exclude_dist`) |

Geometry definition (the in/out test) is in the manuscript §`sec:cyldef` and
`native_cylinder_core.py` (self-tested by `python episcaf_analysis/native_cylinder_core.py`).

## Provenance — honestly

These values were **inherited / hand-set and never swept against ground truth** on record. Two
concrete symptoms of the gap:

- **`exclude_dist` is inconsistent in the code.** `native_cylinder_core.py` defaults it to `4.0`
  (functions at lines ~56, ~102), while every actual DP3 / 12-mer run passes `1.0`
  (`build_12mer_metrics.py`, `dp3_native_cylinder.py`, the presets). The `1.0` runs are the ones
  the manuscript reports; the `4.0` default is stale. **Do not trust the core default — pass 1.0.**
- **`OFFSET = −4` puts the cylinder base 4 Å *below* the epitope**, which the 8pww probe
  (`scripts/cylinder_fp_probe.py`) showed scoops up near-epitope scaffold *below* the paratope,
  inflating the count with false positives (0/13 flagged atoms were within 4 Å of the real
  antibody for design DP2_0804). That is what prompted the sweep.

## The sweep (the reproducible basis going forward)

`scripts/cylinder_param_sweep.py` recomputes the cylinder count over a grid of
`(offset, radius, height, exclude_dist)` in one pass over the DP3 structures and reports, per
grid cell, the **AUC for predicting the real clash** (`af3_n_clash_res == 0`). The geometry is
then chosen on the whole DP3 ground truth, not inherited. Run on Gemini:

```bash
WS=/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff
python3 scripts/cylinder_param_sweep.py \
    --metrics_csv $WS/runs/run_rfd3_mpnn/04_filter/metrics_native_cyl_full.csv \
    --dp2_parquet $WS/datasets/dp2.parquet \
    --native_dir  /tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned \
    --offsets=-6,-4,-2,0,2,4 --radii 12,14,16,18,20 --heights 30,40,50 \
    --exclude_dists 1.0 --limit 20000 \
    --out_csv results/cylinder_param_sweep.csv
```

Output `results/cylinder_param_sweep.csv` has one row per grid cell (offset, radius, height,
variant plain/aware, AUC, median count for clash-free vs clashing). A second, independent check
is the **binding** correlation on the assayed 403 designs at the chosen geometry (within-antibody,
§`sec:whatpredicts`) — the geometry should hold up against real binding, not just the in-silico
clash.

## Chosen values — sweep result (2026-07-01)

Sweep: 20k random sample of the DP3 RFD3 set (`metrics_native_cyl_full.csv`), 180 grid cells
(offsets −6..4, radii 12..20, heights 30/40/50, carve 1.0), scored by AUC for clash-free
(`af3_n_clash_res==0`; 3.8% of the sample). `results/cylinder_param_sweep.csv`.

**Outcome: keep the inherited geometry `(OFFSET −4, RADIUS 16, HEIGHT 40, carve 1.0)`.** It is
near-optimal and no change is justified:

- **offset −4 is the best offset.** Higher offsets (0, +2, +4) do *not* appear in the top-12
  cells; −2 is slightly worse than −4. So despite the 8pww false positives (which motivated the
  sweep), raising the base *hurts* aggregate clash prediction — the 8pww over-count is a
  **local** artifact of that epitope's geometry, not a global mis-setting.
- **radius / height barely matter.** The whole top-12 spans AUC 0.9017–0.9004 (native-aware,
  carve 1.0). The current `(−4,16,40)` = 0.8987; the best cell `(−4,18,40–50)` = 0.9017, a
  **+0.003** gain that is within sample noise (754 positives) and not worth re-running the full
  DP3 + regenerating every figure. Height (30/40/50) is essentially flat.
- native-aware (carve 1.0) beats plain everywhere (`(−4,16,40)`: 0.899 vs 0.874), as expected.

(Absolute AUC here, ~0.90, is on the *ungated* full set; the manuscript's 0.935 is on the
well-predicted subset (epitope-chunk RMSD < 2.5) — different population, both valid. Only the
*relative* ranking across cells is used to choose the geometry.)

### Carve radius (exclude_dist) sweep — 2026-07-01

Same 20k sample, geometry fixed at (−4, 16, 40), sweeping only the native-carve exclusion radius
(`results/cylinder_carve_sweep.csv`):

| exclude_dist | 0.5 | 1.0 (current) | **1.5** | 2.0 | 2.5 | 3.0 | 4.0 (stale) |
|---|---|---|---|---|---|---|---|
| AUC vs clash | 0.880 | 0.899 | **0.908** | 0.903 | 0.890 | 0.879 | 0.866 |

Unlike the geometry (flat), this dial has a **real optimum at 1.5 Å** (AUC rises then falls; 0.042
range — the most consequential cylinder parameter, and it was never tuned). 1.5 beats the current
**1.0 by ~0.009** (clean monotone-then-peak, not noise); ~1.5 Å is a van-der-Waals contact with the
antigen, vs 1.0 carving only near-exact overlaps. The stale **4.0 is the worst of all** (0.866) —
confirms fixing it was right.

**Recommendation: adopt exclude_dist = 1.5 (pending confirmation).** Before changing 1.0→1.5:
(1) re-run the carve sweep on the full DP3 (drop `--limit`) and gated (epitope RMSD < 2.5); (2)
cross-check against binding on the assayed 403 — recompute the assayed cylinder at 1.5
(`assayed_native_cylinder.py --exclude_dist 1.5`) and re-run the within-antibody correlation
(§`sec:whatpredicts`). If both hold, update the default in `native_cylinder_core.py` (currently
1.0) + the preset, re-run the full DP3 native-aware at 1.5, and regenerate the manuscript's cylinder
numbers (the AUC 0.935 and the binding correlations are at 1.0).

**Geometry:** the sweep found offset/radius/height near-optimal at the inherited (−4, 16, 40), gain
noise-level, so **no geometry change**. The stale `exclude_dist=4.0` default in
`native_cylinder_core.py` is fixed to 1.0 (will move to 1.5 if the confirmations above hold).
