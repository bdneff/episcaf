#!/usr/bin/env bash
# run_dual_island_mpnn_af3.sh -- stage 03_mpnn + 04_af3 for the dual-island run.
#
# Run this AFTER the 02_rfd3 array has finished (all backbones in 02_rfd3/outputs). It does the
# cheap login-node fixed-PDB step, then prints the two GPU sbatch waves in order. It does NOT
# submit -- you run the printed sbatch lines, and MPNN (wave 1) must finish before AF3 (wave 2),
# since AF3 reads the MPNN sequences.
#
# Pipeline:  02_rfd3/outputs
#   -> [login] stage03_mpnn_fixed_pdbs  : CIF -> backbone PDB + FIXED epitope  (03_mpnn/fixed_pdbs)
#   -> [gpu  ] stage03_mpnn_submit      : ProteinMPNN, 8 seqs/backbone         (03_mpnn/mpnn_pdbs)
#   -> [login] stage04_af3_emit_jsons   : MPNN PDB -> AF3 single-seq JSON       (04_af3/inputs)
#   -> [gpu  ] stage04_af3_array        : AlphaFold3                            (04_af3/outputs)
#
# NOTE: the ProteinMPNN tool path/env (stage03_mpnn_submit.py) and the AF3 container
# (stage04_af3_array.sbatch) are cluster-specific and have NOT been verified off-cluster.
set -euo pipefail
cd "$(dirname "$0")/.."   # repo root

RUN_DIR="${1:-runs/dual_island_rfd3}"
LEDGER_CSV="${LEDGER_CSV:-results/dual_island_designs.csv}"
N_WORKERS="${N_WORKERS:-8}"
BATCH_SIZE="${BATCH_SIZE:-500}"
THROTTLE="${THROTTLE:-200}"

RFD3_OUT="$RUN_DIR/02_rfd3/outputs"
FIXED_DIR="$RUN_DIR/03_mpnn/fixed_pdbs"
MPNN_DIR="$RUN_DIR/03_mpnn/mpnn_pdbs"
AF3_IN="$RUN_DIR/04_af3/inputs"

n_cif=$(find "$RFD3_OUT" -name '*.cif.gz' 2>/dev/null | wc -l | tr -d ' ')
echo ">> 02_rfd3/outputs has $n_cif model CIFs"
[ "$n_cif" -gt 0 ] || { echo "   no CIFs yet -- has the RFD3 array finished?"; exit 1; }

echo ">> stage03_mpnn_fixed_pdbs: CIF -> FIXED backbone PDB (login node)"
python3 scripts/stage03_mpnn_fixed_pdbs.py \
    --rfd3_outputs_dir "$RFD3_OUT" \
    --ledger           "$LEDGER_CSV" \
    --outdir           "$FIXED_DIR" \
    --n_workers        "$N_WORKERS"
n_fixed=$(find "$FIXED_DIR" -name '*_fixed.pdb' | wc -l | tr -d ' ')
echo "   wrote $n_fixed fixed PDBs -> $FIXED_DIR"

cat <<EOF

============================================================================
Next, the two GPU waves (run on gemini; MPNN must finish before AF3):

  # WAVE 1 -- ProteinMPNN (submits a per-batch sbatch array, 8 seqs/backbone):
  python3 scripts/stage03_mpnn_submit.py \\
      --fixed_pdb_dir $FIXED_DIR \\
      --outdir        $MPNN_DIR \\
      --batch_size    $BATCH_SIZE
  #   (add --dry_run first to stage batches without submitting)

  # WAVE 2 -- after WAVE 1 completes, emit AF3 JSONs then submit AF3:
  python3 scripts/stage04_af3_emit_jsons.py --mpnn_pdb_dir $MPNN_DIR --out_dir $AF3_IN --seed 1
  N=\$(find $AF3_IN -maxdepth 1 -name '*.json' | wc -l)
  TASKS=\$(( (N + 99) / 100 ))   # stage04 packs 100 JSONs per array task
  #   chunk under MaxArraySize if TASKS is large (here it can be ~1100):
  for s in \$(seq 1 1000 \$TASKS); do e=\$((s+999)); [ \$e -gt \$TASKS ] && e=\$TASKS; \\
    sbatch --array=\$s-\$e%$THROTTLE scripts/stage04_af3_array.sbatch $RUN_DIR; done
============================================================================
EOF
