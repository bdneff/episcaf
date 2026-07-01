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
    --offsets -6,-4,-2,0,2,4 --radii 12,14,16,18,20 --heights 30,40,50 \
    --exclude_dists 1.0 --limit 20000 \
    --out_csv results/cylinder_param_sweep.csv
```

Output `results/cylinder_param_sweep.csv` has one row per grid cell (offset, radius, height,
variant plain/aware, AUC, median count for clash-free vs clashing). A second, independent check
is the **binding** correlation on the assayed 403 designs at the chosen geometry (within-antibody,
§`sec:whatpredicts`) — the geometry should hold up against real binding, not just the in-silico
clash.

## Chosen values

`[to record from the sweep]` — after the sweep, write the chosen `(OFFSET, RADIUS, HEIGHT,
exclude_dist)` here with the AUC that justifies each, update the constants in
`native_cylinder_core.py` and `dp3_native_cylinder.py`, fix the stale `exclude_dist=4.0`
defaults, re-run the full DP3 native-aware, and regenerate the cylinder figures/numbers in the
manuscript.
