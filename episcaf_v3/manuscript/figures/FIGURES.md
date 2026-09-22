# episcaf_v3 figures — exact regeneration commands

Every figure in `manuscript/` is regenerable by the command below from committed inputs. Run from
`episcaf_v3/energetics/` with a python that has pandas / scipy / matplotlib.

## Section: Validating an energetic epitope definition

Inputs (all tracked): `md/3hfm/holo_ie_mean.csv` (per-residue interaction-energy means from
`holo_ie.py`), `md/3hfm/mmgbsa_perres.csv` (per-residue MM-GBSA ΔG from `mmgbsa_decomp_to_csv.py`),
`skempi_3hfm_ddg.csv` (SKEMPI antigen alanine-scan ground truth).

- `ie_vs_ddg_3hfm.png` — bare interaction energy (Coulomb + LJ) vs ΔΔG:
  ```
  python plot_ie_vs_ddg.py --ie md/3hfm/holo_ie_mean.csv --ddg skempi_3hfm_ddg.csv \
      --channel ab_total --out ../manuscript/figures/ie_vs_ddg_3hfm.png
  ```
- `ie_vs_ddg_3hfm_lj.png` — Lennard-Jones only (the LJ probe; referenced in prose/table, not shown):
  ```
  python plot_ie_vs_ddg.py --ie md/3hfm/holo_ie_mean.csv --ddg skempi_3hfm_ddg.csv \
      --channel ab_lj --out ../manuscript/figures/ie_vs_ddg_3hfm_lj.png
  ```
- `ie_vs_ddg_3hfm_mmgbsa.png` — desolvation-inclusive MM-GBSA ΔG vs ΔΔG:
  ```
  python plot_ie_vs_ddg.py --ie md/3hfm/mmgbsa_perres.csv --ddg skempi_3hfm_ddg.csv \
      --channel mmgbsa_dg --out ../manuscript/figures/ie_vs_ddg_3hfm_mmgbsa.png
  ```

The upstream `mmgbsa_perres.csv` is itself regenerable — see `energetics/mmgbsa/README.md` and
`docs/DECISIONS.md` (D2): run `mmgbsa/run_mmgbsa.sbatch` on the 3HFM trajectory, `grep '^TDC'` the
complex/ligand sander mdouts into `md/3hfm/{complex,ligand}_tdc.txt`, then `mmgbsa_decomp_to_csv.py`.


## Tamarind pilot and comparison (2026-09-08)

Run from repository root:

```bash
python3 episcaf_v3/energetics/tamarind/pilot.py verify
python3 episcaf_v3/energetics/tamarind/report_pilot.py
MPLCONFIGDIR=/tmp/episcaf-mpl /usr/bin/python3 episcaf_v3/energetics/tamarind/plot_comparison.py
```

`report_pilot.py` reads `energetics/tamarind/3hfm_pilot/results.zip`, frozen request,
manifest and completed status, plus `energetics/skempi_3hfm_ddg.csv`. It generates
`tamarind_pilot_table.tex`, comparison.csv and summary.json. `plot_comparison.py` reads
that comparison and the existing holo_ie_mean.csv/mmgbsa_perres.csv plus the experimental
CSV to generate `tamarind_vs_simulation_3hfm.png` and `.pdf`. Exact plotted values and
correlations are in 3hfm_pilot/plot_points.csv and plot_statistics.json. ML n=2; simulations
n=13. No classification threshold or two-point ML correlation is computed.
Remote preparation/validation/submission and retrieval are recorded in
`episcaf_v3/energetics/tamarind/README.md`; rerunning local analysis does not submit jobs.


Full 13-mutation update (same figure filenames; pilot raw results retained):

```bash
python3 episcaf_v3/energetics/tamarind/report_pilot.py --run-dir episcaf_v3/energetics/tamarind/3hfm_all13 --table-name tamarind_all13_table.tex
MPLCONFIGDIR=/tmp/episcaf-mpl /usr/bin/python3 episcaf_v3/energetics/tamarind/plot_comparison.py --run-dir episcaf_v3/energetics/tamarind/3hfm_all13
```

Reads the separate all13 frozen run and writes its comparison/summary/plot data into
`3hfm_all13/`. Do not pool duplicate predictions from the initial pilot.
