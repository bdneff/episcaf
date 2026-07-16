# results/

Small, **tracked** derived tables that figures and the manuscript depend on (shortlists, ledgers,
summary stats). Big raw data and the full rankings stay on `/tgen_labs` (or are gitignored locally)
and are regenerable; only the small cuts live here. Provenance index below — each file maps to the
script that produces it. Metric-CSV inputs (the `$D` sibling dirs, or the cluster) are documented in
`data/README.md` and `configs/paths.py`.

| file | produced by |
|------|-------------|
| `dp4_C1_whole_epitope_ranked.top20.csv` | `scripts/stage06_select.py --preset antibody --group id --topk 20` (C1-103 metrics) |
| `dp4_C2_single_island_ranked.top20.csv` | `scripts/stage06_select.py --preset antibody --group id,island_index --topk 20` (dual-island metrics, cluster) |
| `dp4_C3_12mer_ranked.top20.csv` | `scripts/stage06_select.py --preset twelvemer --group antigen,id --topk 20` (12-mer metrics) |
| `dp4_C{1,2,3,5}_scaffoldEPITOPE.csv` | `case_encode_selected.sbatch` / `case_encode_whole_epitope.py` / `case_encode_c2.py` / `case_encode_c3.py` (case-encode the selected designs) |
| `dp4_C5_titration.csv` | `scripts/stage06_sample_c5.py` (farthest-point titration over the scoring axes) |
| `dp4_C6_controls.csv` | `episcaf_pipeline/scaffolded_epitope_controls/build_c6_mutants.py` (island→Ala + scaffold disruption) |
| `dp4_8vdl_top10.csv` | `dp4_8vdl/scripts/07_consolidate.py` (PfEMP1 arm, top-10 per definition) |
| `dual_island_designs.csv` | `episcaf_pipeline/build_dual_island_designs.py` (C2 per-island design ledger) |
| `dual_island_targets.csv` | `episcaf_analysis/dual_island_targets.py` (92 island defs; `tab:dualisland`) |
| `dual_island_gate_summary.csv` | `scripts/stage05_summarize.py` (per-island four-filter breakdown; `tab:funnel`) |
| `whole_epitope_designs.csv` | `scripts/build_whole_epitope_designs.py` (C1 native-103 contig ledger, 2,206 contigs) |
| `dp3_binding_metrics.csv` | `scripts/build_dp3_binding_join.py` (DP3 assay binding joined to dp2 + metrics) |
| `assayed_native_cyl.csv`, `assayed_native_cyl_ed1.5.csv` | `scripts/assayed_native_cylinder.py` (cluster; native-aware cylinder on assayed designs, exclude_dist 1.0 and 1.5) |
| `assayed_cylinder_worklist.csv` | `scripts/build_assayed_cylinder_worklist.py` |
| `antigen_seq_vs_fasta.csv` | `scripts/verify_antigen_seq_vs_pdb_fasta.py` (antigen_seq vs PDB-FASTA check for 30-mer tiling) |
| `cylinder_fp/<id>/…` | `scripts/cylinder_fp_probe.py` (8pww false-positive visualization inputs) |

The gitignored full rankings (`results/dp4_*_ranked.csv`, ~200 MB) come from the same
`stage06_select.py` commands without the `--topk` cut. Full run order: `docs/PIPELINE.md`.
