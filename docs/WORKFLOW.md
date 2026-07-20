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
cd manuscript && tectonic main.tex
git add manuscript && git commit -m "manuscript: refresh fp-reduction figure + text" && git push
```

## Cluster job sizing & walltime (size from throughput, with margin)

Always allocate generous SLURM walltime and size batch chunks to the *observed* per-item
throughput — never a guess. Prefer jobs that skip already-done outputs, so a timeout becomes
a cheap backfill rather than a silent hole.

**Worked example — the MPNN timeout (2026-06-21).** The first ProteinMPNN wave for the
dual-island run used `--time=2:00:00` with 500 backbones × 8 sequences per batch. MPNN ran
~215 backbones/h, so 500 needed ~2.3 h and every batch died at the wall ~430/500 — dropping
its tail, so only 97,786 of 111,360 designs survived (~12% lost). Fix: walltime 2 h → 8 h, batch
500 → 300, and a `--skip_done --tag redo` backfill mode in `scripts/stage03_mpnn_submit.py`
that re-runs only backbones missing their 8 outputs.

**Current per-stage walltimes (validate the first task before a full submit):**

| stage | sbatch | walltime | per task | watch |
|-------|--------|----------|----------|-------|
| RFD3 | `episcaf_pipeline/hpc/sbatch/rfd3_array.sbatch` | 10 min | 1 contig → 8 backbones | ok |
| MPNN | `scripts/stage03_mpnn_submit.py` | 8 h | 300 backbones × 8 seqs | fixed (was 2 h) |
| AF3  | `scripts/stage04_af3_array.sbatch` | 4 h | 100 JSONs | **validate** first task's timing on ~111k |

Rule of thumb: estimate runtime = (items × observed per-item time), then multiply by 2–3×.
