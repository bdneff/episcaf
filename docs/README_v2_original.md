# RFD3 + MPNN Scaffold Design Experiment

## Overview

This repository contains the pipeline and analysis for a controlled experiment comparing two epitope scaffolding strategies:

1. **RFD3-only** (`run_test_rfd3_nompmn`): RFD3 all-atom designs submitted directly to AlphaFold3 for validation
2. **RFD3+MPNN** (`run_rfd3_mpnn`): The same RFD3 backbones re-sequenced with ProteinMPNN before AF3 validation

The central scientific question is whether the 14× lower pass rate observed in the RFD3-only pipeline relative to Lawson's RFD1+MPNN pipeline is due to **backbone geometry** (RFD3 designs bad scaffolds) or **sequence distribution** (RFD3 sequences are out-of-distribution for AF3). If the RFD3+MPNN pass rate is comparable to Lawson's, the sequence distribution hypothesis is confirmed.

---

## Background

### Lawson's Pipeline (baseline)
```
RFD1 backbone → ProteinMPNN (8 seqs/backbone) → AF3 → 4-filter validation
```
- 2,360 contigs × 8 RFD1 designs × 8 MPNN sequences = ~151,000 designs
- Pass rate: ~2.7% (403/~15,000 passing the 4 filters)

### RFD3-only Pipeline
```
RFD3 all-atom design (8 designs/contig) → AF3 → 4-filter validation
```
- 2,360 contigs × 8 RFD3 designs = 18,880 designs
- Pass rate: ~0.19% (36/18,880) — approximately 14× lower than Lawson

### This Experiment (RFD3+MPNN)
```
RFD3 backbone (stripped to backbone) → ProteinMPNN (8 seqs/backbone) → AF3 → 4-filter validation
```
- 2,360 contigs × 8 RFD3 backbones × 8 MPNN sequences = ~143,000 designs
- Expected: if sequence distribution is the bottleneck, pass rate should approach Lawson's

---

## Theoretical Motivation

The empirical pass rate is:

$$\hat{P}_{\text{pass}} = \frac{1}{N} \sum_{i=1}^{N} \mathbf{1}[f(s_i, b) = 1]$$

where $f(s, b) = 1$ indicates sequence $s$ passes all four filters given backbone $b$.

ProteinMPNN was trained explicitly on the **inverse folding task**: given a fixed backbone, sample sequences likely to fold into it. This shapes $p_{\text{MPNN}}(s \mid b)$ to concentrate mass on physically plausible sequences for a given geometry.

RFD3 learns the **joint distribution** $p_{\text{RFD3}}(s, b)$, and the conditional $p_{\text{RFD3}}(s \mid b)$ is a byproduct of that training rather than its explicit objective. Crucially, both MPNN and AF3 are trained predominantly on natural PDB proteins, meaning they share a common training distribution. RFD3 sequences, being more generative/de novo, are more out-of-distribution relative to AF3's expectations.

Analysis of passing designs in the Lawson dataset supports this interpretation. For `7ox3_0P` (highest pass rate), passing designs were distributed roughly uniformly across all 8 RFD1 backbones (per-backbone pass rates 13–25%), suggesting scaffold geometry is not the primary determinant of success. MPNN's 424 sequences per backbone drives the productive search.

---

## 4-Filter Validation Criteria

All four filters must pass for a design to be considered successful:

| Filter | Metric | Threshold | Description |
|--------|--------|-----------|-------------|
| 1 | `overall_rmsd` | ≤ 2 Å | Backbone RMSD of full scaffold: designed structure vs AF3 prediction |
| 2 | `epitope_chunk_rmsd` | ≤ 1 Å | Backbone RMSD of epitope chunk only: designed vs AF3 |
| 3 | `mean_pae` | < 5 | Mean predicted aligned error from AF3 confidence JSON |
| 4 | `af3_n_clash_res` | == 0 | No antibody residues within 4 Å of non-epitope scaffold residues |

Clash detection: Kabsch-align AF3 epitope CA onto true epitope CA, then check non-epitope scaffold heavy atoms vs antibody atoms.

---

## Repository Structure

```
repo_refactored/
├── datasets/
│   └── dp2.parquet                    # Lawson's ground truth (151,232 rows)
├── runs/
│   ├── run_test_rfd3_nompmn/          # RFD3-only pipeline run
│   │   ├── 02_rfd3/outputs/           # RFD3 all-atom CIF.gz outputs
│   │   └── 04_filter/
│   │       └── metrics_full.csv       # Computed metrics for all 18,880 designs
│   └── run_rfd3_mpnn/                 # RFD3+MPNN pipeline run (this experiment)
│       ├── 01_fixed_pdbs/             # RFD3 CIFs converted to backbone PDB w/ FIXED remarks
│       ├── 02_mpnn_pdbs/              # ProteinMPNN all-atom PDB outputs
│       │   └── batch_XXXX/            # 38 batches of ~500 backbones each
│       └── 03_af3/
│           ├── inputs/                # AF3 JSON inputs (~143,000 files)
│           └── outputs/               # AF3 CIF predictions (pending)
├── scripts/
│   ├── 01_rfd3_cif_to_fixed_pdb.py   # Step 1: Convert RFD3 CIFs → backbone PDBs with FIXED remarks
│   ├── 01_rfd3_cif_to_fixed_pdb.sh   # SLURM wrapper for step 1
│   ├── 02_submit_mpnn.py             # Step 2: Batch and submit ProteinMPNN jobs
│   ├── 03_emit_af3_jsons_mpnn.py     # Step 3: Generate AF3 JSON inputs from MPNN PDBs
│   ├── 04_af3_array.sbatch           # Step 4: SLURM array job for AF3 inference
│   └── compute_metrics.py            # Step 5: Compute 4-filter metrics on AF3 outputs
├── logs/
│   └── rfd3_mpnn/
│       ├── 01_cif_to_pdb_*.log
│       ├── 02_mpnn/
│       └── 03_af3/
└── README.md                          # This file
```

---

## Pipeline Steps

### Step 1: Convert RFD3 CIFs to Fixed PDBs

Converts all 18,880 RFD3 all-atom CIF.gz outputs to backbone-only PDB files with `REMARK PDBinfo-LABEL` FIXED annotations for the epitope residues. Epitope indices are looked up per-token from dp2.

```bash
sbatch scripts/01_rfd3_cif_to_fixed_pdb.sh
```

**Output:** `runs/run_rfd3_mpnn/01_fixed_pdbs/{token}_pred{N}_fixed.pdb`
**Result:** 18,880 PDBs, 0 failures, completed in ~30 seconds

---

### Step 2: Run ProteinMPNN

Submits 38 SLURM GPU jobs, each processing 500 fixed PDBs with ProteinMPNN (`dl_interface_design.py`). Uses Lawson's exact parameters: `temperature=0.1`, `relax_cycles=0`, `seqs_per_struct=8`.

```bash
python scripts/02_submit_mpnn.py \
    --fixed_pdb_dir runs/run_rfd3_mpnn/01_fixed_pdbs \
    --outdir        runs/run_rfd3_mpnn/02_mpnn_pdbs \
    --batch_size    500
```

Uses Lawson's exact ProteinMPNN installation:
- Script: `/tgen_labs/altin/alphafold3/workspace/dl_binder_design/mpnn_fr/dl_interface_design.py`
- Env: `/tgen_labs/altin/alphafold3/miniconda3/envs/proteinmpnn_binder_design`
- Weights: `ProteinMPNN/vanilla_model_weights/v_48_020.pt`

**Output:** `runs/run_rfd3_mpnn/02_mpnn_pdbs/batch_XXXX/{token}_pred{N}_fixed_dldesign_{M}.pdb`
**Result:** 143,052 all-atom PDBs (18,880 backbones × 8 sequences, minus ~0.01% PyRosetta failures)

---

### Step 3: Generate AF3 JSON Inputs

Extracts sequences from MPNN PDB files (chain A) and writes AF3 JSON input files. Uses seed=1 matching Lawson's pipeline.

```bash
python scripts/03_emit_af3_jsons_mpnn.py \
    --mpnn_pdb_dir runs/run_rfd3_mpnn/02_mpnn_pdbs \
    --out_dir      runs/run_rfd3_mpnn/03_af3/inputs \
    --seed         1
```

**Output:** `runs/run_rfd3_mpnn/03_af3/inputs/{token}_pred{N}_fixed_dldesign_{M}.json`
**Result:** ~143,000 JSON files

---

### Step 4: Run AlphaFold3

Submits a SLURM array job. Each task processes 100 JSONs sequentially on one A100 GPU, amortizing singularity startup overhead.

```bash
# Calculate N first
python3 -c "
import math, pathlib
n = len(list(pathlib.Path('runs/run_rfd3_mpnn/03_af3/inputs').glob('*.json')))
print(f'Total JSONs: {n}')
print(f'Array size:  {math.ceil(n/100)}')
"

sbatch --array=1-N scripts/04_af3_array.sbatch runs/run_rfd3_mpnn
```

Uses Lawson's exact AF3 setup:
- Container: `/tgen_labs/altin/alphafold3/containers/alphafold_3.0.1.sif`
- Model dir: `/ref_genomes/alphafold/alphafold3/models`
- DB dir: `/ref_genomes/alphafold/alphafold3/`
- `--norun_data_pipeline` (stub MSA, single-sequence mode)
- `--num_diffusion_samples 1`

Jobs are resumable — a `_DONE` flag is written on success, allowing safe resubmission of the same array command to recover from failures.

**Output:** `runs/run_rfd3_mpnn/03_af3/outputs/{pred_id}/` containing CIF + confidences JSON

---

### Step 5: Compute Metrics and Compare (TODO)

Run `compute_metrics.py` adapted for the MPNN pipeline to compute the 4 filters on AF3 outputs, then compare pass rates:

| Pipeline | Designs | Expected Pass Rate |
|----------|---------|-------------------|
| Lawson RFD1+MPNN | ~151,000 | ~2.7% |
| RFD3-only | 18,880 | ~0.19% |
| RFD3+MPNN (this run) | ~143,000 | TBD |

If RFD3+MPNN ≈ Lawson, the sequence distribution hypothesis is confirmed.

---

## Key File Paths (Cluster)

| Resource | Path |
|----------|------|
| Lawson dp2 ground truth | `datasets/dp2.parquet` |
| Lawson MPNN PDBs | `/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/sourced_antibody_v1/no_antibody/proteinmpnn/` |
| Lawson AF3 predictions | `/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/sourced_antibody_v1/no_antibody/af3_predictions/` |
| True complex PDBs | `/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned/` |
| dl_binder_design | `/tgen_labs/altin/alphafold3/workspace/dl_binder_design/` |
| ProteinMPNN env | `/tgen_labs/altin/alphafold3/miniconda3/envs/proteinmpnn_binder_design` |
| AF3 container | `/tgen_labs/altin/alphafold3/containers/alphafold_3.0.1.sif` |
| AF3 models | `/ref_genomes/alphafold/alphafold3/models` |

---

## Notes

- **Epitope indices**: `scaffolded_epitope_chunk_resindices` in dp2 vary per contig (different scaffold lengths place the epitope at different positions). Always look up per-token, never hardcode.
- **MPNN FIXED remarks**: 1-based residue numbering for PyRosetta compatibility.
- **PAE scale**: AF3 confidences JSON values are in raw units (multiply × 0.01 to match Lawson's h5-stored values — confirmed empirically).
- **Clash alignment**: We align the AF3 epitope onto the true epitope (not MPNN epitope as Lawson did). This is actually more consistent for the RFD3 pipeline since no MPNN PDB is available. Validated on 50 designs passing the other 3 filters: 50/50 agreement with Lawson's stored values.
- **RFD3 seed**: All RFD3 designs used `seed0, rep0`. The token (MD5 hash) encodes the contig identity.
