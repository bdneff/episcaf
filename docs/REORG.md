# REORG — design notes & deployment

> **Historical snapshot** (one-time reorg record). For the current layout see `README.md`, for the
> run order see `docs/PIPELINE.md`. Kept for the reasoning behind the structure; not a live guide.

The realized layout (see `README.md`) and the reasoning behind it, plus how to put
this repo on the cluster **without moving the 426 G of data**.

## Principles

1. **Code in git, data on `/tgen_labs`, joined only by `configs/paths.py`.** The repo
   contains no structures, AF3/MPNN outputs, or metrics CSVs.
2. **Nothing in the data dir moves.** Reorg is a code operation.
3. **One config-driven scorer** replaces the `rank_*`/`apply_*`/`add_*`/`scan_*` pile.
4. **Canonical vs archive is explicit and reversible** (everything preserved; re-promote with `git mv`).

## Why `episcaf_pipeline/` stayed at the root (not `src/episcaf/`)

The generation package is already clean, and `cli.py stage05` shells out to
`<repo_root>/legacy_steps/05_rmsd_vs_af3.py`. Moving the package under `src/` would
break that path and force edits to working generation code that can't be tested here.
So the package + its one legacy dependency were kept intact, and the messy part — the
metrics/scoring/viz scripts — got the clean home (`episcaf_analysis/`). Unifying under
`src/episcaf/` is a tested phase-2 task (see `docs/MIGRATION.md`).

## The scorer (the "feed-forward net, no backprop")

Dials live in `episcaf_analysis/presets.py`, one preset per dataset. Per metric:
`weight`, `better` (`low`/`high`), and a `transform`:

```
percentile : rank(pct=True)                              # assumption-light, population-relative
minmax     : (x - min) / (max - min)
zscore     : (x - mean) / std
sigmoid    : 1 / (1 + exp( s * k * (x - midpoint)))      # s=+1 if better="low" else -1
identity   : x
```

Each metric becomes one input neuron with its own activation; the weighted sum is the
output; `gate` drops failures first and `select` keeps top-k per epitope. `sigmoid`
gives an absolute, population-independent score with a tunable threshold (`midpoint`)
and steepness (`k`) — use it when you want scores that don't shift as the population
changes.

**v2 (the backprop):** the DP3 set has `is_pass` (4-filter ground truth). Fit a
logistic regression of the transformed features on `is_pass` to *learn* the weights;
carry them to the 12-mer preset where no antibody label exists. Leave the hand-set
weights until that's worth doing.

## Git + safe deploy (data never moves)

Local, after cleanup:
```bash
cd episcaf
git init && git add -A && git commit -m "reorg: clean structure + config-driven scorer"
git remote add origin git@github.com:<you>/episcaf.git   # durable off-cluster backup
git push -u origin main
```
If the cluster can't reach GitHub, push to a bare repo over SSH instead:
```bash
# once on cluster:  git init --bare ~/episcaf.git
git remote add cluster ssh://bneff@<cluster>/home/bneff/episcaf.git && git push cluster main
```

On the cluster — clone into a **new** directory, never onto the data dir:
```bash
# stage in home first to confirm, per your plan:
git clone <remote> ~/episcaf && cd ~/episcaf
$EDITOR configs/paths.py            # already points at episcaf_v2_bneff data
python tests/test_scoring.py
python -m episcaf_analysis.score --preset twelvemer \
    --metrics-csv "$(python -c 'import configs.paths as p; print(p.METRICS_12MER)')" \
    --out /tmp/composite_12mer_top5.csv     # compare to your last known top-5
```
Once verified, the repo can live wherever you like (home, or a sibling of
`episcaf_v2_bneff` on `/tgen_labs`). The old loose code in `episcaf_v2_bneff` can then
be deleted — it's tiny and now lives in git; **the data dirs stay**.

**Never:** `git init` inside `episcaf_v2_bneff`, or `rsync --delete` toward it.

## Validation checklist (on the cluster)

- [ ] `python -c "import configs.paths as p; print(p.METRICS_12MER, p.METRICS_ANTIBODY)"` resolves.
- [ ] `score.py --preset twelvemer` reproduces your current top-5 shortlist.
- [ ] `score.py --preset antibody` runs against the DP3 metrics.
- [ ] `python tests/test_scoring.py` passes; `python episcaf_analysis/native_cylinder_core.py` passes.
- [ ] Work through the “verify against the data” list in `docs/MIGRATION.md`.
- [ ] `git status` clean — no data/CSV/structure files staged.
