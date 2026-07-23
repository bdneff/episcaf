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
The next library, **DP4**, is assembled from seven episcaf components (C1–C6 plus the 8VDL
PfEMP1 arm; manuscript `sec:dp4`) — 14,241 episcaf constructs — and combined with a **21,759-row LX
PfEMP1/EPCR minibinder arm** into one file, `data/libraries/dp4_library.csv` (**36,000 rows**), handed
off whole to the LadnerLab oligo encoder (`episcaf_pipeline/oligo_encoding/`).

**DP4 status (2026-07-23): built, encoded, gated, and synthesis-ready.** The episcaf library was selected
under the soft-gate scorer (`antibody_softgate`, epitope-PAE midpoint 2.5, manuscript `sec:composite`); the
whole library (episcaf + minibinders, encoded together into one PepSeq assay) was encoded to DNA and gated
into `data/libraries/dp4_order_file.csv` — 36,000 oligos, each 349 nt with the 20-mer adapters, every core
verified to translate back to exactly its own peptide. (Assembled at 37,083, culled to 35,962 on
2026-07-21 — C6 controls to top-15, dedup of the picked-twice, drop of 60 no-accessibility designs — then
topped up to 36,000 on 2026-07-23 by deepening the 8VDL arm to top-29/run; see `docs/DP4_LIBRARY.md` →
*Cull to 35,962* and *Top-up to 36,000*.) That file is what goes to Twist, with
`data/libraries/dp4_quote_file.csv` as its 2-column vendor-quote view. An all-designs
superset (357,789 rows, every candidate arm: C1/C2/C3 + 8VDL + passing minibinders) holds every candidate
design, not just the ones that shipped. It's a true superset for the *candidate-pool* arms — every shipped
C1/C2/C3/8VDL/minibinder design is in it (the C4/C5/C6 controls aren't candidate designs, so they aren't).
The gzipped copy is committed at `data/libraries/dp4_superset.csv.gz` (~34 MB); the raw
`.csv` is regenerated on `$WS`. Full reference:
`docs/DP4_LIBRARY.md`; step order: `docs/PIPELINE.md`.

> **Code lives in git; data and runs live on `/tgen_labs`.** Nothing under `runs/`,
> `run_12mer_scaffolding/`, or `datasets/` is committed (see `.gitignore`). The code
> finds the data through `configs/paths.py` — edit that one file per environment.
>
> `$WS` throughout the docs is the durable cluster workspace,
> `/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff`; `$REPO` is the cluster checkout,
> `/scratch/bneff/episcaf`. **`/scratch` is ephemeral** (swept on ~30 days) — the checkout there is
> disposable, but anything long-lived must be copied to `$WS`.

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

# 8. ASSEMBLE  [local]  — concatenate the episcaf components into data/libraries/dp4_library.csv
#                         (33-col schema; 14,241 episcaf rows), then fold in the LX minibinders -> 36,000
python scripts/stage06_assemble.py --depth 20                # C1/C2 top-20, C3 top-10, 8VDL top-29 (14,241 episcaf)
python dp4_8vdl/scripts/08_add_minibinders.py --lx dp4_8vdl/data/LX_20260626.csv   # -> 36,000 rows
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
#    -> 36,000 oligos, all verified (349 nt, 20-mer adapters, every core translates back). To Twist.
#    Growing the library? Encode ONLY the new peptides: scripts/stage07_new_peptides.py, then merge.
#    The encoder input MUST be LF -- CRLF makes it silently encode nothing (see docs/DP4_LIBRARY.md).
python scripts/stage07_quote_file.py --order-file data/libraries/dp4_order_file.csv \
    --out data/libraries/dp4_quote_file.csv   # 2-col name,349mer for the vendor quote

# 10. SUPERSET  [cluster + local, analysis]  — every candidate design across all arms (357,789), not
#     just what shipped. build_superset.sbatch (C1/C2/C3) then extend_superset.py (+8VDL +minibinders).
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
PfEMP1 conserved-epitope arm (`dp4_8vdl/`). These are the **14,241 episcaf** constructs (after the 2026-07-21 cull and the 2026-07-23 top-up; see docs/DP4_LIBRARY.md). Folded in
alongside them is a **21,759-row LX PfEMP1/EPCR minibinder arm** (de-novo binders on the same antigen,
from a separate LatentX effort — not episcaf-scored, carried with their own `lx_*` metrics), so the
shipped `data/libraries/dp4_library.csv` is **36,000 rows** — two projects in one PepSeq assay. All rows
are the constant 103-mer; the file carries the full 33-column schema and encodes to DNA oligos
(`episcaf_pipeline/oligo_encoding/`, stage 07). **Full reference — components, selection math/weights,
exclusions, the minibinder arm, the column dictionary, status: `docs/DP4_LIBRARY.md`.**

## Where the run data lives (cluster)

The design runs behind this library are on the TGen cluster, not in git. The durable copies live under
the workspace `$WS = /tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff` (persistent); the working
checkout `$REPO = /scratch/bneff/episcaf` is disposable and `/scratch` is swept on ~30 days, so treat the
`$WS` paths below as the source of truth. Only three components (C1, C2, C3) plus the 8VDL arm have their
own RFD3→MPNN→AF3 output; C4/C5/C6 are derived from them (see `docs/DP4_LIBRARY.md`) and carry no
separate run.

| Component | Run directory (durable, under `$WS`) | Per-design metrics |
|---|---|---|
| C1 whole-epitope | `runs/whole_epitope_rfd3/` | `runs/whole_epitope_rfd3/05_analysis/metrics_whole_epitope.csv` |
| C2 single-island | `runs/dual_island_rfd3/` | `runs/dual_island_rfd3/05_analysis/metrics_dual_island.parquet` |
| C3 polyclonal 12-mer | `run_12mer_scaffolding/` | `run_12mer_scaffolding/06_score/metrics_12mer.csv` |
| 8VDL PfEMP1 arm | `dp4_8vdl/` | (top-10 per definition consolidated into `results/`) |
| Oligo encoding | `runs/dp4_encoding_full/` (`DP4_best_encodings`, `out_seqs`, `output_ratio`) | — |
| All-designs superset | raw `dp4_superset.csv` regenerated here (357,789 rows, all arms); **gzipped copy committed to git** at `data/libraries/dp4_superset.csv.gz` (~34 MB) | — |
| DP3 design table | `datasets/dp2.parquet` (Lawson's ledger; C5/C6 trace back to the C1 pool) | — |

C4/C5/C6 have no cluster run: **C4** is built from the antigen FASTA sequences, and **C5** (metric-space
titration) and **C6** (alanine-scan + scaffold-disruption controls) are both derived from the C1
`whole_epitope_rfd3` design pool. Their small selection files are tracked in `results/`. The shipped
deliverables — `dp4_library.csv`, `dp4_named_peptides.csv`, and the verified Twist order file
`dp4_order_file.csv` — live in git under `data/libraries/`, so they do not depend on the cluster at all.

## Intended cleanups (deferred)

The code analog of the manuscript's open questions: structural improvements we've *chosen to defer*,
recorded so they aren't lost. None are bugs — they're friction noticed while extending the repo, and
we deliberately did not refactor while the generation runs were in flight (the C1-103 and 8VDL scaffolds,
both now complete and shipped).

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
  conserved-epitope breadth) rather than pipeline chronology. The runs it was waiting on have landed, so
  this can now be anchored to real results whenever it's taken up.

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
