# PIPELINE — the end-to-end run order

The single authoritative sequence for taking an epitope from a crystal structure to a synthesizable
DNA oligo. Each step names the canonical tool; where several tools exist for one step (different run
layouts), the alternatives are noted. Data locations live in `configs/paths.py`; the manuscript
counterpart is `sec:methods`.

**Paths.** Cluster steps assume `REPO=/scratch/bneff/episcaf` (the git checkout — disposable, it is a
clone) and `WS=/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff` (the durable workspace).
`/scratch` is ephemeral and swept on ~30 days; anything long-lived belongs under `$WS`. Full map:
memory `filesystem-map`.

**A note on the numbering.** Three numbering schemes coexist for historical reasons: the `stage0x`
script-name prefixes (`scripts/stage01…stage07`), the `episcaf_pipeline/stages/` package
(`stage01,02,04,05`), and the human 0–9 steps in `README.md`. They do **not** line up one-to-one
(e.g. the `stage03_mpnn_*` scripts are the MPNN sub-step of "sequences", while `stages/stage04` is the
AF3-direct variant). This table is the source of truth; the prefixes are just filenames.

| # | Step | Where | Canonical tool | Notes / alternatives |
|---|------|-------|----------------|----------------------|
| 0 | Init / snapshot dataset | local | `python -m episcaf_pipeline prep` (or `init`) | Snapshot `dp2.parquet`, set up the run dir. |
| 1 | Compile contigs | local | `episcaf_pipeline/stages/stage01_compile_contigs.py` | Turn each epitope's islands into RFD3 contig strings. |
| 2 | Emit RFD3 inputs | local | `episcaf_pipeline/stages/stage02_emit_rfd3_inputs.py` | Per-contig input files; needs the AbDb complex PDBs (`ABDB_CLEANED_PDB_DIR`). |
| 3 | RFD3 backbones | cluster GPU | `episcaf_pipeline/hpc/sbatch/rfd3_array.sbatch` | 8 backbones per contig. See `gemini-run-ops` for partitions. |
| 4 | Sequences | cluster GPU | **MPNN path** (canonical): `scripts/stage03_mpnn_fixed_pdbs.py` → `stage03_mpnn_submit.py` → `stage04_af3_emit_jsons.py`. **AF3-direct** alt: `stages/stage04_emit_af3_jsons.py`. | ProteinMPNN holds the epitope fixed, 8 seqs/backbone → the `*_fixed.pdb` naming contract. |
| 5 | AF3 predictions | cluster GPU | `episcaf_pipeline/hpc/sbatch/af3_array.sbatch` (or `scripts/stage04_af3_array.sbatch`) | Single-sequence, `--norun_data_pipeline`, seed 1. |
| 6 | Metrics | local/cluster | `scripts/stage05_extract_metrics.py` (dual-island); `episcaf_analysis/compute_metrics.py` (DP3/Lawson, `--validate`); `episcaf_pipeline/cli.py stage05` (→ `legacy_steps/05_rmsd_vs_af3.py`) | One builder per run layout — they differ by run-dir shape, not by result. RMSDs, PAE, clash, cylinder. |
| 7 | Score + select | local | `scripts/stage06_select.py --preset antibody_softgate\|twelvemer --group … --topk …` (wraps `episcaf_analysis.score`) | Rank, no hard gate. **`antibody_softgate` is the adopted scorer** (C1/C2 **and 8VDL** — all known-antibody, real Fab clash; epitope-PAE midpoint 2.5). `twelvemer` (percentile, cylinder accessibility) is **C3 only**. C5 is not composite-ranked at all — it's `stage06_sample_c5.py` farthest-point sampling over the C1 known-Ab pool. The bare `antibody` preset is the superseded percentile scorer. Dials in `presets.py`. |
| 8 | Assemble library | local | `scripts/stage06_assemble.py --depth 20` → `data/libraries/dp4_library.csv` (15,324 episcaf) | Concatenate the seven episcaf components; 56-mAb exclusion + depth cuts + global numbering. Ships the full 33-column schema. |
| 8b | Fold in minibinders | local | `dp4_8vdl/scripts/08_add_minibinders.py --lx <LX.csv>` → 37,083 rows | Append the 21,759 passing LX minibinders (with their `lx_*` columns) so the file is the whole DP4 library. Idempotent. |
| 9 | Export encoder input | local | `scripts/stage07_named_peptides.py` → `dp4_named_peptides.csv` | `name,seq`, no header, all 37,083 (add `--exclude-category minibinder` for episcaf-only). |
| 10 | Encode → order file | cluster | `oligo_encoding/encode_step1_generate.sbatch` → `encode_step2_select.sbatch`; then `scripts/stage07_order_file.py` | Peptide → DNA (LadnerLab, DP3 codon weights) → Twist-adapter order file, row-verified (349 nt, 20-mer adapters, every core translates back). **Pin `ADAPTER`** — the tool's default is the wrong 19-mer form (see `docs/DP4_LIBRARY.md` / memory `oligo-adapter-trap` for why 20-mers). |
| 11 | All-designs superset | cluster | `sbatch scripts/build_superset.sbatch` (one pass: builds C1/C2/C3 → `extend_superset.py` folds in 8VDL + minibinders + unions columns → gzip) | Optional/analysis, not on the deliverable path. Every candidate design across all arms (**357,789** rows × 36 cols), not just what shipped, → `$WS/dp4_superset.csv` and the committed `data/libraries/dp4_superset.csv.gz`. Reads AF3 chain A via `--af3-remap` at the durable `$WS` copies; needs the LX source on-cluster (fails loudly if missing). |

## Component-specific notes

- **C1 (whole epitope)** rebuilt natively at 103: `scripts/build_whole_epitope_designs.py` →
  `run_whole_epitope_rfd3.sh` → `run_whole_epitope_mpnn_af3.sh` (steps 1–5 in one wrapper each).
- **C2 (single island)**: the per-island run; metrics `metrics_dual_island.parquet` (cluster).
- **C3 (polyclonal 12-mer)**: no-antibody arm; accessibility via the cylinder surrogate (step 6/7).
- **C4 (linear 30-mer controls)**: no design — `episcaf_pipeline/build_dp4_tiled30mers_fasta.py` only.
- **C5 (titration)** / **C6 (mutant controls)**: derived from the C1 pool (steps 7/8), not fresh runs.
- **8VDL (PfEMP1 arm)**: self-contained under `dp4_8vdl/` (`01_contigs` → … → `07_consolidate.py`),
  reuses the shared MPNN/AF3 via the `*_fixed_dldesign_*` contract.

Reproduce commands for the DP4 deliverables specifically: `docs/DP4_LIBRARY.md`. Figure commands:
`manuscript/figures/FIGURES.md`. Result-file provenance: `results/README.md`.
