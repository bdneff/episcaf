# Paired computational alanine scanning (staged 2026-09-08)

Purpose: estimate mutation-induced binding changes, rather than using WT residue
decomposition as a surrogate. For each matched WT snapshot, gmx_MMPBSA constructs an
alanine mutant and evaluates both endpoint binding estimates:

`ddG(frame) = (Gcomplex-Gantibody-Gantigen)mut - (Gcomplex-Gantibody-Gantigen)WT`.

Positive weakens binding. This is a fixed-geometry, single-trajectory endpoint
approximation. It does not sample the mutant or independent unbound ensemble, perform
alchemical integration, or include a configurational entropy estimate. Added/rebuilt
alanine atoms follow the tool's mutation construction; retained protein geometry is
inherited from WT. Lower variance from paired snapshots cannot remove model bias.

## Frozen scope

All 13 antigen alanine mutations from `../skempi_3hfm_ddg.csv`; identity checked against
the committed 3HFM PDB. `stage_scan.py` generates one config per mutant and manifest.json.
Use frames 1,21,...,1981 (100 frames), matching the staged denser GB resample; actual
frame count must be checked after execution. GB igb=5, salt=0.150 M, temperature=300 K,
mbondi2 radii (PBRadii=3), intdiel=1, extdiel=78.5, surften=0.0072, surfoff=0.
These explicit settings follow the current GB setup and documented defaults, not a
new fitted model. `mutant_only=0` calculates both WT and mutant for each task;
`cas_intdiel=0` prevents mutation-specific dielectric reassignment. No decomposition.
WT recomputation costs extra but keeps each initial paired calculation self-contained.

## Execution boundary and current verification

**Staged, not executed.** Brandon runs Gemini; no cluster has been accessed. Locally:

```bash
python3 episcaf_v3/energetics/alanine_scan/stage_scan.py
bash -n episcaf_v3/energetics/alanine_scan/run_scan.sbatch
```

All 13 reference selections were checked. The runner/CSV parser still need verification
against the installed gmx_MMPBSA and actual output. The existing topology/protonation
state is reused; preserve generated WT/mutant Amber topologies and tool logs. Check that
only the intended antigen residue is mutated, that chain Y numbering matches the
reference, and that the number/order of residues matches the MD topology. The supplied
crystal reference establishes chain IDs and residue numbering, not trajectory coordinates.

## Commands for Brandon (Gemini)

From the repository root, after pulling this commit, set the actual local run paths:

```bash
export ALASCAN_CODE="$PWD/episcaf_v3/energetics/alanine_scan"
export ALASCAN_SOURCE=/scratch/bneff/episcaf_run/episcaf_v3/energetics/md/3hfm/out
export ALASCAN_OUT=/scratch/bneff/3hfm_alanine_scan_20260908
sbatch "$ALASCAN_CODE/run_scan.sbatch"
```

The source path is the previously documented example; point it at the real existing run.
Required: md.tpr, md_whole.xtc, topol.top + all topology includes, gbsa_index.ndx with
unique groups named `antibody` and `antigen`. The runner resolves group numbers from names.
If md_whole.xtc is absent, generate it in the source directory using the established
`gmx trjconv -s md.tpr -f md.xtc -o md_whole.xtc -pbc whole` with group 0, then verify the
complex is whole. Never treat a partially written trajectory as complete.

The array processes tasks 0–12, one at a time, serial within each task. Partition and GPU
allocation follow the working GB recipe (GPU idles); 12 h is a per-task cap, not a runtime
prediction. New directories are required and existing results are never deleted. Relative
topology files are copied; installed force-field includes remain environment dependencies.
The runner records hashes of topology inputs, trajectory and TPR, copies the config and
reference, and captures GROMACS/gmx_MMPBSA versions plus explicit conda packages. Retain
all generated Amber topologies and MD inputs under persistent $WS before scratch cleanup.

## Analysis and return files

For each completed mutation, e.g. K96A:

```bash
python "$ALASCAN_CODE/analyze_scan.py" --run "$ALASCAN_OUT/Y96_ALA"
```

The parser reads WT and mutant **Delta Energy Terms** from FRAME_ENERGIES.csv, pairs by
frame ID, requires all 100 planned frames, and reports mutant-minus-WT values and per-term
differences. It rejects missing/mismatched frames. Block means and SEM are reported for
5,10,20,25 consecutive sampled frames per block; this is sensitivity analysis, not proof
of independence or an uncertainty interval covering model error. Upstream CSV energies
are rounded to 0.01 kcal/mol; preserve raw mdouts for higher precision if needed.

Bring back each mutation's FRAME_ENERGIES.csv, FINAL_RESULTS.dat, run.log,
run_manifest.json, versions.txt, paired_ddg.csv, paired_summary.json, fixed complex PDB,
and WT/mutant Amber topologies. Compare all 13 to experiment and StaB-ddG without tuning
on the intended held-out set. PB/alternative dielectric/relaxation are later decisions.

## Method sources (accessed 2026-09-08)

- https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/examples/Alanine_scanning/
- https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/input_file/
- CSV implementation: GMXMMPBSA/output_file.py and amber_outputs.py in
  https://github.com/Valdes-Tresanco-MS/gmx_MMPBSA (inspected current master).
