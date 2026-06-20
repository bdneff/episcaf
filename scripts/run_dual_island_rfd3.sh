#!/usr/bin/env bash
# run_dual_island_rfd3.sh -- stage the dual-island per-island RFD3 run, then print the sbatch.
#
# The RFD3 step is the FIRST step of the RFD3+MPNN pipeline; it only generates backbones.
# MPNN -> AF3 -> scoring come later (see docs / the per-island design plan, manuscript
# sec:perisland). This script does the cheap login-node prep and hands you the array command;
# it does NOT submit, so you stay in control of the cluster.
#
# Usage (on gemini, after `git pull` and `conda activate ~/rfd3/env/rfd3_py312`):
#   bash scripts/run_dual_island_rfd3.sh [RUN_DIR]
# Override the knobs inline if you want, e.g.:
#   SEEDS=0,1,2,3 REPS=2 bash scripts/run_dual_island_rfd3.sh runs/dual_island_$(date +%Y%m%d)
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root

# --- knobs ---------------------------------------------------------------------------------
RUN_DIR="${1:-runs/dual_island_rfd3}"
# 8 RFD3 designs per island-contig = the established protocol (8 designs/contig). Array size
# is 87 islands x (#seeds x reps). Default: 8 seeds x 1 rep = 696 tasks.
SEEDS="${SEEDS:-0,1,2,3,4,5,6,7}"
REPS="${REPS:-1}"
# Per-epitope antigen PDBs (named <id>.pdb). configs.paths ABDB_CLEANED_PDB_DIR.
PDB_DIR="${PDB_DIR:-/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned}"
LEDGER_CSV="results/dual_island_designs.csv"          # tracked, shipped with the repo
LEDGER_PARQUET="results/dual_island_designs.parquet"  # materialized here (gitignored)
# -------------------------------------------------------------------------------------------

echo ">> materializing ledger parquet from $LEDGER_CSV (no dp2 needed)"
python3 - "$LEDGER_CSV" "$LEDGER_PARQUET" <<'PY'
import sys, ast, pandas as pd
src, dst = sys.argv[1], sys.argv[2]
df = pd.read_csv(src)
df["epitope_resindices"] = df["epitope_resindices"].apply(ast.literal_eval)
df.to_parquet(dst, index=False)
print(f"   {len(df)} island rows -> {dst}")
PY

# NB: the all-in-one `prep` subcommand is currently broken (cmd_prep references an
# undefined args.in_parquet), so we run init + stage02 + stage03 separately -- same result.
echo ">> init $RUN_DIR"
python3 -m episcaf_pipeline init --dataset "$LEDGER_PARQUET" --run_dir "$RUN_DIR" --force
echo ">> stage02 (expand seeds x reps)  [seeds=$SEEDS reps=$REPS]"
python3 -m episcaf_pipeline stage02 --run_dir "$RUN_DIR" --seeds "$SEEDS" --reps "$REPS"
echo ">> stage03 (emit RFD3 inputs)  [pdb_dir=$PDB_DIR]"
python3 -m episcaf_pipeline stage03 --run_dir "$RUN_DIR" --pdb_dir "$PDB_DIR"

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
