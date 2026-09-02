# DP4 results analysis — the plan for when the assay data lands

**Status (2026-08-05):** the DP4 PepSeq library is ordered; synthesis + assay is a ~4-week
turnaround, so the binding data does not exist yet. This doc + `scripts/build_dp4_binding_join.py`
are staged now so that the day the data arrives we run one command and get the answer, instead of
writing the analysis from scratch under time pressure. It is deliberately narrow: **being ready to
read the DP4 results.** John's separate "broader epitope space" investigation (the 1,134-PDB island
landscape) is a different track — see the last section.

This mirrors what DP3 already did once: its assay data came back, got joined to the design metrics,
and produced the "what predicts binding" result. DP4 is the same loop, better positioned.

---

## TL;DR — the moment the data lands

1. Drop the assay CSV(s) into `data/dp4_binding/` (make the dir).
2. Confirm the **one unknown**: the assay column layout (see "The join", below).
3. Run the join:
   ```
   python scripts/build_dp4_binding_join.py --assay data/dp4_binding/<run>.csv
   ```
   → `results/dp4_binding_metrics.csv` (one row per assayed peptide, binding readout + every
   design metric + composite score attached).
4. Run the analysis — **staged and smoke-tested; no adaptation needed:**
   ```
   python scripts/analyze_dp4_binding.py
   ```
   → prints Q1–Q5 (per-antibody enrichment; pooled-vs-within-antibody correlation of binding with
   composite / rank / cylinder / each metric; 8VDL accessible-vs-blocked; scaffolded-vs-linear) and
   writes `manuscript/figures/dp4_binding_scatter.png` + `dp4_metric_binding.png`.

Two commands, not one — the honest version: the join gives you the table, the analysis gives you the
answers. Both are done; only the join's *input column mapping* waits for the data (below).

Validate the harness **today** with no data — this already passes:
```
python scripts/build_dp4_binding_join.py --self-check
```

**Local run note:** on the laptop use `/usr/bin/python3` (it has pandas/numpy/scipy/matplotlib);
the Homebrew `python3`/`python` on PATH does **not** have pandas. On the cluster use the
`rfd3_py312` conda env. The scripts are plain and need only those four libraries.

**Arm ↔ `category` string** (the plan's C-labels are these exact `category` values, which the code
keys on): C1 `scaffoldedAbEpitope` · C2 `scaffoldedSingleIsland` · C3 `scaffoldedPolyclonal` ·
C4 `tiled30mer` · C5 `metricSpaceTitration` · C6 `scaffoldedEpitopeControl` · 8VDL `scaffolded8VDL`
· minibinder. `composite`/metrics are populated on C1/C2/C5/8VDL and **blank on C3/C4/C6/minibinder**
(never folded / no antibody), so metric regressions run only where the axis is non-null.

---

## What DP3 found (the precedent, and why DP4 improves on it)

DP3 shipped binding intensities for 8 antibodies across two runs. Readout = **cognate
log-enrichment**, `log10(1+Ab) − log10(1+NoAb)`; a binder sits above the `y=x` line on a
`log10(1+NoAb)` vs `log10(1+Ab)` scatter. The join attached the AF3 metrics and asked which one
predicts binding. Result (manuscript `sec:whatpredicts`, `tab:metricbinding`, n=377 cognate
designs):

- The **native-aware cylinder** was the single best predictor and the **only** metric whose
  correlation survived deconfounding: pooled Pearson −0.24 → within-antibody **−0.22 (p=2e−5)**,
  same negative sign in every well-sampled antibody. Negative sign is mechanistically right —
  scaffold mass in the antibody's path means less binding.
- Every other metric collapsed within-antibody: epitope RMSD −0.41 → −0.14; overall RMSD −0.08;
  mean PAE +0.07 (wrong sign, noise).

**Three limits made DP3 a prior, not a fit:**
1. Everything assayed had already passed the 4-filter, so the filter metrics had **no variance**
   (overall_rmsd 0.35–1.95, mean_pae 3.0–5.0, clash all 0) — no failing designs to learn from.
2. The pooled signal was a between-antibody artifact; only within-antibody centering is honest.
3. Only two antibodies were assayed deeply (7ox3 n=193, 5fhx n=149).

**DP4 was built to break all three** (`sec:open`, `docs/DP4_LIBRARY.md`):
- **C5 `metricSpaceTitration` (2,715 designs)** deliberately samples the *full* metric space,
  **including designs the current filters reject.** That is the variance DP3 lacked — it turns the
  regression from a prior into an actual **fit of the composite dials** (`episcaf_analysis/presets.py`).
- The **8VDL arm (58 designs)** is a clean *designed* accessibility contrast: `8VDL_epitope`
  (full window, accessible, low clash) vs `8VDL_hotspots` (minimal graft, scaffold blocks the Fab,
  high clash). A prospective test of the accessibility hypothesis on one epitope.
- Many arms, one library — broader antibody coverage than DP3's 8.

---

## The join

DP4 is **much simpler than DP3**: DP3 round-tripped through `dp2.parquet` for metrics and ran a
separate cylinder job. In DP4, **every metric and the composite score are already columns in
`data/libraries/dp4_library.csv`** — so the join is a single merge on `library_member`. No parquet,
no cluster round-trip, no separate cylinder run.

**Verified now** (`--self-check`, against the shipped 36,000-row library):
- `library_member` is unique and contiguous `DP4_1..DP4_36000` — the master join key.
- `target` carries the cognate id for every antibody arm as `<pdb>_<N>P` (e.g. `2qqn_0P`),
  exactly like DP3's `Target`. `8VDL_epitope`/`8VDL_hotspots` and the EPCR minibinders
  (`fold_pfemp1_epcr_*`) sit on their own targets.
- All 15 metric/scoring columns the analysis needs are present: `epitope_rmsd, overall_rmsd,
  epitope_pae, scaffold_pae, mean_pae, ptm, af3_clashes, cylinder_clashes, composite,
  rank_in_group, is_global_pass, island_index, category, target, design_ID`.

**The ONE unknown until the data arrives** — the assay CSV's column layout:
- Does it keep a `library_member` column? If yes, the join is the trivial merge. If not, the script
  falls back to joining on the peptide sequence (`designedSequence` or the 103-mer `sequence`).
- The **NoAb baseline** column name, and the **per-antibody intensity** column names. DP3 used
  `NoAb…` for baseline and `X<pdb>_…` for each antibody; those regexes are the script defaults.
  Override with `--noab-re` / `--ab-re`, or supply an explicit `--ab-map` file (one `pdb,column`
  line per antibody) if DP4's naming differs.

That is the only thing to eyeball when the CSV shows up. Everything downstream is wired.

---

## The questions DP4 answers (the analysis, per arm)

Arm definitions are authoritative in `docs/DP4_LIBRARY.md`; counts below are the shipped library.

| Question | Arm(s) | What to compute |
|---|---|---|
| **Q1. Does our composite score select binders?** (the headline validation) | C1/C2 scaffolded antibody arms | enrichment vs `composite` and vs `rank_in_group`, **within-antibody centered**. Do high-composite designs bind more? |
| **Q2. Does the cylinder predict binding — prospectively?** | wherever there's cognate-antibody ground truth | re-run DP3's within-antibody regression of enrichment on `cylinder_clashes`. DP3 found it best in-sample; DP4 tests it out-of-sample and on more antibodies. |
| **Q3. Weight fit (not just a prior).** | **C5 `metricSpaceTitration` (2,715)** | regress enrichment on all metric axes (`epitope_rmsd, overall_rmsd, mean_pae, cylinder_clashes, af3_clashes`) with **real variance + failing designs**, within-antibody. This is what sets the dials in `presets.py`. |
| **Q4. Accessibility contrast.** | **8VDL (58: 29 epitope + 29 hotspots)** | do accessible `8VDL_epitope` designs bind and blocked `8VDL_hotspots` not? A prospective clash test. Small N — a contrast, not a fit. **Trap:** 8VDL's `target` is `8VDL_*` → cognate key `8vdl`, but its antibody is **C7** and the join's default regex is lowercase-only. If the assay names it `C7`, the join finds no 8VDL cognate — re-run the join with `--ab-map` mapping `8vdl,<C7 column>`. `analyze_dp4_binding.py` prints an explicit note when this happens. |
| **Q5. Does scaffolding beat a linear peptide?** | C1 vs **C4 `tiled30mer` (2,033)** | matched epitopes, scaffolded vs linear control. Compare enrichment for the same epitope presented both ways. (C4 rows carry **no** `composite`/metrics — this is an enrichment-only comparison, by design.) |
| **Q6. Per-epitope baseline.** | **C6 `scaffoldedEpitopeControl` (2,325)** | control set; use as the per-epitope reference when reading Q1/Q3. (Also metric-blank — enrichment only.) |
| **(separate readout) Minibinders.** | `minibinder` (21,759), `scaffoldedPolyclonal` (4,390) | these bind EPCR/PfEMP1 and the polyclonal tiled antigens (`1d2k`/`4wat`/`6m0j`), **not** the mAb panel — a different readout, keyed on their own targets, not the cognate-mAb enrichment above. Analyze separately. |

**The deconfounding rule (carry it over from DP3, don't skip it):** never trust a pooled
correlation. Subtract each antibody's mean from both the metric and the enrichment (fixed-effect
centering per `cognate_ab`) before correlating. Pooled numbers mix in each antibody's baseline
offset and manufacture signal — that is exactly how DP3's epitope-RMSD −0.41 evaporated to −0.14.

---

## Reproducibility

| Script | Role | Status |
|---|---|---|
| `scripts/build_dp4_binding_join.py` | assay CSV → `results/dp4_binding_metrics.csv` (join on `library_member`, attach metrics/composite, cognate log-enrichment) | **staged; `--self-check` passes** |
| `scripts/analyze_dp4_binding.py` | reads the join output → Q1–Q5 report + `dp4_binding_scatter.png` + `dp4_metric_binding.png` | **staged; smoke-tested on a synthetic join output (2026-08-05)** |

Both scripts are staged and verified. The analysis is **assay-layout-independent**: it reads only
the join's fixed-schema output (`cognate_ab`, `cognate_log_enrichment`, `composite`, the metric
columns), never the raw assay columns — so it needed no guessed column names and is panel-agnostic
(it groups by whatever antibodies appear). It was smoke-tested by feeding it a synthetic
`dp4_binding_metrics.csv`: it recovered the planted composite/cylinder correlations with the correct
signs, produced both figures, and correctly skipped the metric-blank arms. The DP3 originals
(`analyze_dp3_binding.py`, `plot_dp3_metric_binding.py`) remain as the reference the DP4 script
consolidates. The two figures do **not** exist yet — they are produced the first time
`analyze_dp4_binding.py` runs on real data.

**Known cosmetic data note (John, 2026-08-05):** a subset of minibinder rows carry a truncated
`target`/`lx_target` value `fold_pfemp1_epcr_mod` (should read `…_model_0`), e.g. `DP4_35177`. It
is a label-field truncation only — it does not touch the peptide or the oligo, and minibinders are
analyzed on their own EPCR readout, not the cognate-mAb join — so it does not affect the results
analysis. Flagged here so it doesn't cause confusion when reading the minibinder rows.

---

## Separate track — John's broader epitope-space investigation (NOT part of results-readiness)

John also asked (Slack, 2026-08-05) to use the wait to characterize the **broader epitope
landscape**: the 59 (≤2-island) subset was drawn from ~1,134 cleaned AbDb complexes, and he wants
to know how much our results generalize beyond ≤2 islands. This is a distinct piece of work from
reading the assay data; it is **not** built here, only pointed to so it isn't lost. Scoped questions
from the thread:

- Distribution of island count (`epitope_chunks`) across all 1,134 — and how it shifts if islands
  are "collapsed" by filling gaps of up to ~5 residues.
- Relationship between island count and number of contact residues (≤4 Å footprint).
- Filters worth applying: drop epitopes involving non-protein elements (glycosylation, modified
  amino acids, small molecules); drop epitopes too small to reliably scaffold/RMSD (e.g. 7a3t).
- Identify a **non-island-filtered subset (~50 PDBs)** as a second focus for a future design round.

**Where the data lives** (from the repo audit): island count is the pre-computed `epitope_chunks`
column in `dp2.parquet` (local copy at `known_antigen/analysis/full_run/dp2.parquet`), so the
island-count histogram and per-island spans run **fully locally** off that ledger
(`episcaf_analysis/dual_island_targets.py`, `episcaf_pipeline/build_dual_island_designs.py`). A
fresh contact footprint or a distance-based island definition across all 1,134 needs the raw
`ABDB_CLEANED_PDB_DIR` PDBs (cluster-only), using `dp4_8vdl/scripts/contact_epitope.py`'s 4 Å
heavy-atom + `contiguous_runs` logic as the template. **No non-protein/HETATM filter exists yet** —
that code would be new. This track should get its own plan doc when it's picked up.
