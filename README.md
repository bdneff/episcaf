# episcaf

Epitope scaffolding: design protein scaffolds that display a target epitope the way
it sits on the native antigen (RFdiffusion3 → ProteinMPNN → AlphaFold3), then score
designs on epitope fidelity and antibody accessibility.

Two settings:
- **antibody / DP3 set** — antigens with a known antibody; has ground-truth
  `af3_n_clash_res` and the 4-filter `is_pass`.
- **tiled 12-mer set** — no antibody; a native-antigen-aware *cylinder* surrogate
  stands in for the antibody-clash term.

DP3 has now been **assayed**: experimental binding for 8 mAbs (`data/dp3_binding/`)
shows the native-aware cylinder is the strongest predictor of binding among all our
metrics — the result that sets the scorer's priors (manuscript `sec:whatpredicts`).
The next library, **DP4**, is assembled from seven components (C1–C6 plus the 8VDL
PfEMP1 arm; manuscript `sec:dp4`) into `data/libraries/dp4_library.csv` (15,324 constructs)
and handed off to the LadnerLab oligo encoder (`episcaf_pipeline/oligo_encoding/`).

**DP4 status (2026-07-16): built, encoded, and synthesis-ready.** The library was selected under the
soft-gate scorer (`antibody_softgate`, manuscript `sec:composite`), encoded to DNA, and gated into
`data/libraries/dp4_order_file.csv` — 15,324 oligos, each 349 nt with the 20-mer adapters, every core
verified to translate back to exactly its own peptide. That file is what goes to Twist. An all-designs
superset (`$WS/dp4_superset.csv`, 334,750 rows) holds every candidate design, not just the ones that
shipped, for looking at the distributions the library was drawn from. Full reference:
`docs/DP4_LIBRARY.md`; step order: `docs/PIPELINE.md`.

> **Code lives in git; data and runs live on `/tgen_labs`.** Nothing under `runs/`,
> `run_12mer_scaffolding/`, or `datasets/` is committed (see `.gitignore`). The code
> finds the data through `configs/paths.py` — edit that one file per environment.

## Layout

```
episcaf/
├── episcaf_pipeline/     # generation package (RFD3→AF3): init, stage02–05, cli
│   └── oligo_encoding/   #   stage 07: peptide→DNA oligo encoding (LadnerLab tool + DP3 recipe)
├── legacy_steps/         # 05_rmsd_vs_af3.py — live dependency of pipeline stage05
├── episcaf_analysis/     # metrics, scoring, viz
│   ├── score.py          #   THE config-driven composite scorer
│   ├── presets.py        #   the scoring dials (one preset per dataset)
│   ├── native_cylinder_core.py   # cylinder geometry (native-aware clash surrogate)
│   ├── build_12mer_metrics.py    # 12-mer metrics builder (imports native_cylinder_core)
│   ├── compute_metrics.py        # DP3/no-MPNN 4-filter metrics builder
│   ├── false_positive_check.py / composite_swap_validate.py   # validation tools
│   └── viz/              #   distribution / fp-reduction / cylinder plots
├── scripts/              # 12-mer/MPNN pipeline-branch step scripts (01→04)
├── configs/paths.py      # absolute /tgen_labs data locations (the only path file)
├── data/                 # small tracked inputs: dp3_binding/ (assay CSVs), libraries/ (tiled sets)
├── results/              # small derived tables figures depend on (tracked)
├── manuscript/           # the living project log (LaTeX → main.pdf); build with `tectonic main.tex`
├── tests/                # scorer unit tests (no data needed)
├── archive/              # every superseded/duplicate/legacy script, preserved
└── docs/                 # REORG.md, MIGRATION.md, original experiment writeup
```

## Quickstart

```bash
# environment (cluster: use the rfd3_py312 env; it has MDAnalysis/gemmi/scipy)
pip install -e .                      # or: pip install -r requirements.txt

# point at the data
$EDITOR configs/paths.py              # set WORKSPACE / metrics CSV paths

# score (dials are in episcaf_analysis/presets.py)
python -m episcaf_analysis.score --preset twelvemer \
    --metrics-csv "$(python -c 'import configs.paths as p; print(p.METRICS_12MER)')" \
    --out /tgen_labs/.../run_12mer_scaffolding/06_score/composite_12mer_top5.csv

# tests
python tests/test_scoring.py
```

## Pipeline, end to end

One epitope → many scaffolded designs → metrics → the best few → a synthesis-ready
oligo. Stages live in a numbered **run directory** on `/tgen_labs`; code drives them via
the `episcaf_pipeline` CLI and a few cluster `sbatch` array jobs. `[cluster]` = runs on
Gemini; `[local]` = runs anywhere. Most commands take `--help`.

```bash
# 0–2. PREPARE  [local]  — init (snapshot dataset) + stage01 (contigs, seeds, reps)
#                          + stage02 (emit RFD3 JSON inputs). One call does all three:
python -m episcaf_pipeline prep \
    --dataset <dataset.parquet> --run_dir <run> --pdb_dir <cleaned_antigen_pdbs>
#   (or run `init` / `stage01` / `stage02` separately; see `... <cmd> --help`)

# 3. RFD3  [cluster GPU array]  — diffuse scaffold backbones around the fixed epitope
sbatch --array=1-N episcaf_pipeline/hpc/sbatch/rfd3_array.sbatch <run>   # N from 02_rfd3 manifest

# 4. SEQUENCES + AF3 INPUTS  [cluster]  — two branches:
#   (a) RFD3-direct (no MPNN):
python -m episcaf_pipeline stage04 --run_dir <run>           # emit AF3 JSONs straight from RFD3
#   (b) RFD3 → ProteinMPNN → AF3 (the canonical path); see each script's --help for args:
python scripts/stage03_mpnn_fixed_pdbs.py  --help            # epitope-fixed backbones for MPNN
python scripts/stage03_mpnn_submit.py      --help            # design sequences on the backbones
python scripts/stage04_af3_emit_jsons.py   --help            # emit AF3 JSONs from MPNN seqs

# 5. AF3  [cluster GPU array]  — predict each design's structure
sbatch --array=1-N episcaf_pipeline/hpc/sbatch/af3_array.sbatch <run>

# 6. METRICS  [cluster]  — epitope/overall RMSD, PAE, clash per design
python -m episcaf_pipeline stage05 --run_dir <run>            # (or scripts/stage05_metrics.sbatch)
#   no-antibody set also needs the accessibility surrogate:
python scripts/dp3_native_cylinder.py ...                    # cylinder_native_aware

# 7. SCORE + SELECT  [local]  — composite rank, top-k per group, NO hard gate
python -m episcaf_analysis.score --preset twelvemer|antibody \
    --metrics-csv <metrics.csv> --out <top5.csv>

# 8. ASSEMBLE  [local]  — concatenate the DP4 components into one ordered 8-column
#                         annotated peptide file (data/libraries/dp4_library.csv, 15,324)
python scripts/stage06_assemble.py --depth 20                # C1/C2 top-20, C3 top-10
python scripts/stage07_named_peptides.py \
    --library data/libraries/dp4_library.csv --out data/libraries/dp4_named_peptides.csv

# 9. ENCODE  [cluster]  — peptide → DNA oligo (LadnerLab tool, DP3 codon weights) → order file
#    Run both steps in the SAME working dir; step 1 is the long pole (hours). The adapters are
#    PINNED to the 20-mers inside encode_step2_select.sbatch: the tool's own --adapter default is
#    the 19-mer form, which silently yields 347-nt oligos. Never rely on that default.
sbatch episcaf_pipeline/oligo_encoding/encode_step1_generate.sbatch   # candidates
sbatch episcaf_pipeline/oligo_encoding/encode_step2_select.sbatch     # NN-pick the best
python scripts/stage07_order_file.py --best-encodings <rundir>/DP4_best_encodings \
    --peptides data/libraries/dp4_named_peptides.csv --out data/libraries/dp4_order_file.csv
#    -> 15,324 oligos, all verified (349 nt, 20-mer adapters, every core translates back). To Twist.

# 10. SUPERSET  [cluster, analysis]  — every candidate design (334,750), not just the 15,324 that
#     shipped, for looking at the distributions. Must run BEFORE the /scratch run dirs are deleted.
sbatch scripts/build_superset.sbatch                          # -> $WS/dp4_superset.csv
```

The full authoritative step order (with the `stage0x` script names) is in `docs/PIPELINE.md`.

Verify-as-you-go: `compute_metrics.py --validate` checks our metrics reproduce Lawson's
stored values before any selection is trusted; the manuscript records every figure's
exact command in `manuscript/figures/FIGURES.md`.

## Scoring model

The scorer is a single-layer perceptron with a per-metric activation and hand-set
weights ("feed-forward, no backprop"): each metric is transformed
(`percentile` / `minmax` / `zscore` / `sigmoid` / `identity`), oriented so higher =
better, then weighted-summed; the **top-k per group are kept by ranking — no hard gate**
(a threshold would only discard designs the ranking already buries). The accessibility
term is the real `af3_n_clash_res` when the antibody is known and the cylinder surrogate
when it is not. Weights are a **prior set from the DP3 binding data** (cylinder + epitope
RMSD strongest); the plan is to *fit* them against **experimental binding** once DP4 spans
the metric space (manuscript `sec:open`). Tune by editing `episcaf_analysis/presets.py`.

## DP4 library (what we are shipping)

Seven components, each answering a design question (manuscript `sec:dp4`): **C1** known-Ab
whole-epitope scaffolds; **C2** single-island scaffolds (87 contigs, the islands test);
**C3** polyclonal-tiling scaffolds; **C4** linear tiled-30mer controls; **C5** metric-space
sampling; **C6** scaffolded-epitope controls (island→Ala + scaffold-disruption); plus the **8VDL**
PfEMP1 conserved-epitope arm (`dp4_8vdl/`). All emit the constant 103-mer in the 8-column annotated
format, assemble into one ordered file (`data/libraries/dp4_library.csv`, 15,324 constructs), and encode
to DNA oligos (`episcaf_pipeline/oligo_encoding/`, stage 07). **Full reference — components, selection
math/weights, exclusions, the 104→103 trim, status: `docs/DP4_LIBRARY.md`.**

## Intended cleanups (deferred)

The code analog of the manuscript's open questions: structural improvements we've *chosen to defer*,
recorded so they aren't lost. None are bugs — they're friction noticed while extending the repo, and
we deliberately do not refactor while a run is in flight (the current C1-103 and 8VDL scaffolds).

- **Unify the contig generators.** `episcaf_pipeline/build_dual_island_designs.py`,
  `build_whole_epitope_designs.py`, and `dp4_8vdl/scripts/01_generate_contigs.py` each re-implement the
  same idea — distribute a scaffold budget into flanks/gaps around fixed islands and emit a contig
  string. Fold them into one parametric generator. Low priority / real regression risk; do it only if a
  *fourth* generator would otherwise be written.
- **Two input protocols, not one dialect.** Rather than collapse the two conventions, formalize them:
  a **standard protocol** for epitopes from IEDB / AbDb — already cleaned, antigen-on-chain-A +
  `epitope_resindices`, which is how the bulk of our epitopes arrive (we have many formatted this way) —
  plus a **custom-protocol** option for one-off structures that keep their own chain/numbering (like
  8VDL: chain-C + `scaffold_segs`, bridged today by `dp4_8vdl/scripts/04_make_fixed_pdbs.py`). The custom
  path also leaves 8VDL's run-dir layout split (`dp4_8vdl/02_rfd3/<run>/` vs `dp4_8vdl/runs/<run>/03_mpnn/`).
  Note only — not building now.
- **Group `results/` and `scripts/`** if they keep growing — both are flat and getting busy (many
  `dp4_*` CSVs; many `stage0*`/`build_*`/`case_encode_*`). Navigable by naming for now.
- **Manuscript reorganization** around the scientific questions (islands: one vs both; scaffold vs linear;
  conserved-epitope breadth) rather than pipeline chronology — deferred until the current runs land so the
  reorg is anchored to real results, not premature.

## Verify

```bash
python tests/test_scoring.py                      # scorer unit tests (no data)
python episcaf_analysis/native_cylinder_core.py   # cylinder geometry self-test
cd manuscript && tectonic main.tex                # build the living log
```

See `docs/MIGRATION.md` for exactly what moved where,
`docs/README_v2_original.md` for the original experimental background and the
4-filter definitions, and `docs/CYLINDER_PARAMS.md` for the cylinder geometry
parameters (provenance + the ground-truth sweep that chooses them).
