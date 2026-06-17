# episcaf_pipeline

A reproducible, manifest-driven pipeline for epitope-scaffold generation and evaluation:

**Design ledger (Parquet)** → **expanded contigs** → **RFDiffusion3 inputs/outputs** → **AlphaFold3 inputs/outputs** → **RMSD + QC metrics**

This repo is organized to be:
- **reproducible** (runs are snapshots + deterministic artifacts),
- **easy to use** (single CLI entrypoint),
- **cluster-friendly** (manifests + SLURM array templates),
- **collaborator-friendly** (clear directory layout + schema contract).

---

## Repository layout

```
datasets/                     # canonical source-of-truth parquets (not run-specific)
runs/                         # generated run artifacts (snapshot + outputs)
episcaf_pipeline/             # pipeline code (stages + CLI)
legacy_steps/                 # preserved, proven scripts (reference + compatibility)
```

Each run uses a fixed layout:

```
runs/<run_name>/
  00_input/      designs.parquet              # snapshot copy of dataset
  01_design/     contigs.parquet              # expanded contigs (seed/rep + design_id)
  02_rfd3/       inputs/ + inputs_manifest.csv + outputs/ + logs/
  03_af3/        inputs/ + inputs_manifest.csv + outputs/ + logs/
  04_analysis/   rmsd_vs_af3_all.csv + rmsd_vs_af3_best_per_run.csv
```

---

## Install / environment

Minimum for stages 02–04:
- Python 3.10+
- pandas + pyarrow (parquet IO)

For stage 05 (RMSD):
- `gemmi`
- `MDAnalysis`
- `numpy`

Example:

```bash
pip install -r requirements.txt
```

---

## Quickstart

### 1) Create a run snapshot (copies the dataset parquet into the run)

```bash
python -m episcaf_pipeline init \
  --dataset datasets/dp2.parquet \
  --run_dir runs/run_20260220_120000
```

### 2) Compile contigs (seed/rep expansion + stable design_id)

```bash
python -m episcaf_pipeline stage02 \
  --run_dir runs/run_20260220_120000 \
  --seeds 0,1,2,3 \
  --reps 1
```

### 3) Emit RFD3 JSON inputs + manifest

```bash
python -m episcaf_pipeline stage03 \
  --run_dir runs/run_20260220_120000 \
  --input_pdb /abs/path/to/antigen_clean.pdb
```

### 4) Run RFD3 on SLURM (array)

```bash
# N = number of rows in 02_rfd3/inputs_manifest.csv (minus header)
sbatch --array=1-N episcaf_pipeline/hpc/sbatch/rfd3_array.sbatch runs/run_20260220_120000
```

### 5) Emit AF3 JSON inputs from the RFD3 outputs

```bash
python -m episcaf_pipeline stage04 --run_dir runs/run_20260220_120000
```

### 6) Run AF3 on SLURM (array)

```bash
sbatch --array=1-N episcaf_pipeline/hpc/sbatch/af3_array.sbatch runs/run_20260220_120000
```

### 7) RMSD analysis (RFD3 vs AF3)

```bash
python -m episcaf_pipeline stage05 --run_dir runs/run_20260220_120000
```

---

## Schema contract

The pipeline enforces a strict separation:
- **SAFE columns** (non-null across the design ledger): may be used by stages 02–04
- **RESULT columns** (nullable): only stage 05+ may consume/produce these

See: `episcaf_pipeline/schema.py`

---

## Notes for collaborators

- Runs are intended to be **reproducible snapshots**; never hardcode absolute paths inside Python stages.
- SLURM scripts read their work items from manifest CSVs so array indices remain stable.
- If you add columns to the dataset parquet, update `schema.py`.

---

## License / internal use

(Insert your preferred license / internal notice.)
