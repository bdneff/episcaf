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
The next library, **DP4**, is assembled from six components (manuscript `sec:dp4`)
and handed off to the LadnerLab oligo encoder (`episcaf_pipeline/oligo_encoding/`).

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

Six components, each answering a design question (manuscript `sec:dp4`): **C1** known-Ab
whole-epitope scaffolds; **C2** single-island scaffolds (87 contigs, the islands test);
**C3** polyclonal-tiling scaffolds; **C4** linear tiled-30mer controls; **C5** metric-space
sampling; **C6** island alanine mutants. All emit the constant 103-mer in DP2 format,
assemble into one ordered file (`06_library`), and encode to DNA oligos
(`episcaf_pipeline/oligo_encoding/`, stage 07).

## Verify

```bash
python tests/test_scoring.py                      # scorer unit tests (no data)
python episcaf_analysis/native_cylinder_core.py   # cylinder geometry self-test
cd manuscript && tectonic main.tex                # build the living log
```

See `docs/MIGRATION.md` for exactly what moved where, and
`docs/README_v2_original.md` for the original experimental background and the
4-filter definitions.
