# WORKFLOW — cluster + local, one repo, living manuscript

Two machines, one git repo. Heavy compute on the cluster; figures and writing local;
everything syncs through git.

```
   cluster (gemini)                         local (laptop)
   ----------------                         --------------
   RFD3 / MPNN / AF3                         regenerate figures from results/
   metric builds over full data   <--git-->  update manuscript/ sections + figures
   big outputs -> /tgen_labs                 build manuscript/main.pdf
   small summaries -> results/               push
   push
```

## The loop
1. **Cluster** runs heavy jobs, produces metrics on `/tgen_labs`, and writes the *small*
   derived tables figures need into `results/` (e.g. a top-k shortlist, summary stats).
   Commit those and push.
2. **Local** pulls, regenerates every figure from `results/` via the viz scripts,
   updates the relevant `manuscript/` section, rebuilds `main.pdf`, commits, pushes.
3. Cluster pulls before the next run so both sides stay current.

## Data policy (keeps the repo small and reproducible)
- **Big data stays on `/tgen_labs`** and is never committed (`.gitignore`).
- **Only small derived artifacts** that a figure or table depends on go in `results/`
  (rule of thumb: < a few MB, itself reproducible from a cluster command).
- A figure is never hand-edited. It is produced by a script from `results/`, and that
  command is recorded in `manuscript/figures/FIGURES.md`. If you can't name the command
  that makes a figure, it doesn't go in the manuscript.

## Typical commands
```bash
# cluster: after a run, snapshot the small result the figure needs
python -m episcaf_analysis.score --preset twelvemer \
    --metrics-csv "$(python -c 'import configs.paths as p; print(p.METRICS_12MER)')" \
    --out results/composite_12mer_top5.csv
git add results/composite_12mer_top5.csv && git commit -m "results: refresh 12mer top5" && git push

# local: regenerate the figure + rebuild manuscript
git pull
python episcaf_analysis/viz/plot_fp_reduction.py --in results/... --out manuscript/figures/fp_reduction.png
cd manuscript && latexmk -pdf main.tex
git add manuscript && git commit -m "manuscript: refresh fp-reduction figure + text" && git push
```
