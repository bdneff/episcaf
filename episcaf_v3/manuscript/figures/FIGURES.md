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
