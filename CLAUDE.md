# CLAUDE.md — episcaf operating contract

Loaded at the start of every Claude Code session. This defines *how* to work in this
repo. For *what* the project is, see `README.md` and `docs/`.

## The project, in one paragraph
Epitope scaffolding: design protein scaffolds that present a target epitope the way it
sits on the native antigen (RFdiffusion3 → ProteinMPNN → AlphaFold3), then score designs
on epitope fidelity and antibody accessibility. Two settings: an **antibody / DP3 set**
(ground-truth `is_pass` + `af3_n_clash_res`) and a **tiled 12-mer set** (no antibody; a
native-aware cylinder surrogate stands in). PI: John Altin.

## How to work here (the contract)
- **One step at a time.** One change, verified, then the next. Don't batch unrelated
  changes into a single response or commit.
- **Reproducible by construction.** Every number and figure must be regenerable by a
  named script from written-down inputs. Make a figure → commit its script and record
  the exact command (in `manuscript/figures/FIGURES.md`).
- **Verify against data; don't assume.** Numbers come from running code on real files,
  never from memory or plausibility. Show the command and its output.
- **Defensible over fast.** Prefer the correct, checkable answer to the quick one. State
  assumptions explicitly; mark anything not yet verified as `[UNVERIFIED]`.
- **No silent magic.** Never invent column names, paths, or thresholds. If something is
  unknown, look it up in the data or ask — don't paper over it.
- **Plan before large or destructive actions.** Use plan mode for multi-file changes,
  refactors, or anything touching data. Never `rsync --delete` toward `/tgen_labs`;
  never `git init` inside a data directory.
- **Keep the manuscript current.** When a result changes, update `manuscript/` the same
  session: the section text, the figure, and the command that produced it. The
  manuscript is the living record, not an afterthought.
- **Check prior art before reconstructing.** Before reverse-engineering a mechanism from
  raw data, look for an existing script, a reference implementation, or the relevant
  `manuscript/sections/*.tex` — especially for anything carried over from DP3. Reach for
  the existing method before inventing one; a long data hunt is a signal the premise is wrong.

## Where things live
- Code: `episcaf_pipeline/` (generation, untouched), `episcaf_analysis/` (metrics,
  scoring, viz), `scripts/` (pipeline steps). Data paths live in `configs/paths.py` only.
- Data & big outputs: `/tgen_labs` (cluster, **persistent**), never in git (`.gitignore`).
  Durable cluster workspace: `$WS=/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff`
  (datasets, runs). `/scratch` is **ephemeral** — keep long-lived artifacts under `$WS`, not
  `/scratch`. Full map in memory `filesystem-map`.
- Small derived results that figures depend on: `results/` (tracked).
- Living manuscript: `manuscript/` (LaTeX → `manuscript/main.pdf`).
- Reference docs: `docs/REORG.md`, `docs/MIGRATION.md`, `docs/WORKFLOW.md`.

## Commands that verify
- `python tests/test_scoring.py`                      # scorer unit tests (no data)
- `python episcaf_analysis/native_cylinder_core.py`   # cylinder geometry self-test
- `python -m episcaf_analysis.score --preset twelvemer --metrics-csv <csv> --out <csv>`
- `cd manuscript && latexmk -pdf main.tex`            # build the manuscript

## Fixed facts
- Cluster env: `conda activate ~/rfd3/env/rfd3_py312` (has MDAnalysis/gemmi/scipy).
- DP3 4-filter pass: overall_rmsd<=2, epitope_chunk_rmsd<=1, mean_pae<5, af3_n_clash_res==0.
- Scorer weights/transforms are dials in `episcaf_analysis/presets.py` (not in code).
- Code lives in git; `/tgen_labs` holds data; they join only via `configs/paths.py`.
