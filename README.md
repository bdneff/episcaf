# episcaf

Epitope scaffolding: design protein scaffolds that display a target epitope the way
it sits on the native antigen (RFdiffusion3 → ProteinMPNN → AlphaFold3), then score
designs on epitope fidelity and antibody accessibility.

Two settings:
- **antibody / DP3 set** — antigens with a known antibody; has ground-truth
  `af3_n_clash_res` and the 4-filter `is_pass`.
- **tiled 12-mer set** — no antibody; a native-antigen-aware *cylinder* surrogate
  stands in for the antibody-clash term.

> **Code lives in git; data and runs live on `/tgen_labs`.** Nothing under `runs/`,
> `run_12mer_scaffolding/`, or `datasets/` is committed (see `.gitignore`). The code
> finds the data through `configs/paths.py` — edit that one file per environment.

## Layout

```
episcaf/
├── episcaf_pipeline/     # generation package (RFD3→AF3): init, stage02–05, cli  [UNTOUCHED]
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
better, then weighted-summed; designs are gated first and the top-k per epitope are
kept. Tune it by editing `episcaf_analysis/presets.py`. To *learn* the weights
instead, fit a logistic regression of the transformed features against the DP3
`is_pass` label — see `docs/REORG.md`.

See `docs/MIGRATION.md` for exactly what moved where, and
`docs/README_v2_original.md` for the original experimental background and the
4-filter definitions.
