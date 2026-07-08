#!/usr/bin/env bash
# run_whole_epitope_rfd3.sh -- stage the C1 (whole-epitope) RFD3 run at 103, then print the sbatch.
#
# This is the C1 redo at 103 (see docs/DP4_LIBRARY.md, manuscript sec:production): our first C1
# pool reproduced Lawson's DP3 run and inherited his 104-residue contigs, but the assay length is
# 103 and DP3 saw weaker signal on 104->103-truncated peptides, so we regenerate C1 natively at
# 103 rather than trim after the fact. The ledger (results/whole_epitope_designs.csv) is Lawson's
# exact whole-epitope contigs with one scaffold residue removed (build_whole_epitope_designs.py);
# C5 and C6 rebuild off this new pool once it is scored.
#
# The RFD3 step only generates backbones; MPNN -> AF3 -> scoring come later
# (scripts/run_whole_epitope_mpnn_af3.sh). This does the cheap login-node prep and hands you the
# array command; it does NOT submit, so you stay in control of the cluster.
#
# Usage (on gemini, after `git pull` and `conda activate ~/rfd3/env/rfd3_py312`):
#   bash scripts/run_whole_epitope_rfd3.sh [RUN_DIR]
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root

# --- knobs ---------------------------------------------------------------------------------
RUN_DIR="${1:-runs/whole_epitope_rfd3}"
# RFD3 emits 8 backbones per task (n_batches=1) -- the established 8-designs/contig depth, so ONE
# task per contig, not 8 seeds. Array size = #contigs x (#seeds x reps); the 56-mAb 103 ledger is
# 2,206 contigs x 1 x 1 = 2,206 tasks (each making 8 backbones -> 141,184 designs after 8x MPNN).
SEEDS="${SEEDS:-0}"
REPS="${REPS:-1}"
# Per-epitope antigen PDBs (named <id>.pdb). configs.paths ABDB_CLEANED_PDB_DIR.
PDB_DIR="${PDB_DIR:-/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned}"
LEDGER_CSV="results/whole_epitope_designs.csv"          # tracked, shipped with the repo
LEDGER_PARQUET="results/whole_epitope_designs.parquet"  # materialized here (gitignored)
# -------------------------------------------------------------------------------------------

echo ">> materializing ledger parquet from $LEDGER_CSV (no dp2 needed)"
python3 - "$LEDGER_CSV" "$LEDGER_PARQUET" <<'PY'
import sys, ast, pandas as pd
src, dst = sys.argv[1], sys.argv[2]
df = pd.read_csv(src)
df["epitope_resindices"] = df["epitope_resindices"].apply(ast.literal_eval)
df.to_parquet(dst, index=False)
print(f"   {len(df)} contigs -> {dst}")
PY

# init (00_input) + stage01 (01_design contigs) + stage02 (02_rfd3 inputs). Run explicitly rather
# than via `prep` so each stage's output is visible.
echo ">> init $RUN_DIR"
python3 -m episcaf_pipeline init --dataset "$LEDGER_PARQUET" --run_dir "$RUN_DIR" --force
echo ">> stage01 (compile contigs, expand seeds x reps)  [seeds=$SEEDS reps=$REPS]"
python3 -m episcaf_pipeline stage01 --run_dir "$RUN_DIR" --seeds "$SEEDS" --reps "$REPS"
echo ">> stage02 (emit RFD3 inputs)  [pdb_dir=$PDB_DIR]"
python3 -m episcaf_pipeline stage02 --run_dir "$RUN_DIR" --pdb_dir "$PDB_DIR"

MANIFEST="$RUN_DIR/02_rfd3/inputs_manifest.csv"
N=$(($(wc -l < "$MANIFEST") - 1))   # minus header
CHUNK="${CHUNK:-1000}"              # SLURM MaxArraySize is often ~1001; keep chunks under it
THROTTLE="${THROTTLE:-200}"         # cap tasks running at once (the %N suffix)
echo
echo ">> staged $N RFD3 input JSONs in $RUN_DIR/02_rfd3/inputs"
if [ "$N" -le "$CHUNK" ]; then
  echo ">> launch the RFD3 array with:"
  echo
  echo "   sbatch --array=1-$N%$THROTTLE episcaf_pipeline/hpc/sbatch/rfd3_array.sbatch $RUN_DIR"
else
  echo ">> $N tasks exceeds one array (SLURM MaxArraySize ~$CHUNK). Submit in chunks:"
  echo
  echo "   for s in \$(seq 1 $CHUNK $N); do e=\$((s+$CHUNK-1)); [ \$e -gt $N ] && e=$N; \\"
  echo "     sbatch --array=\$s-\$e%$THROTTLE episcaf_pipeline/hpc/sbatch/rfd3_array.sbatch $RUN_DIR; done"
  echo
  echo "   (check your limit first: scontrol show config | grep MaxArraySize ; set CHUNK below it)"
fi
echo
echo ">> after RFD3 finishes:  bash scripts/run_whole_epitope_mpnn_af3.sh $RUN_DIR"
echo