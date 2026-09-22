# Tamarind binding-ddG pilot

One StaB-ddG job, authorized and submitted 2026-09-08:
`episcaf-v3-3hfm-stabddg-K96A-R73A-20260908`.
Completed: 2026-09-08 20:06:59 (API timestamp). Execution 30 s after 183 s in queue;
0.01 weighted hours (not a dollar cost). K96A = +1.7630341 and R73A = +0.48058623
kcal/mol. See `3hfm_pilot/comparison.csv`; this selected pair is not a validation set.

Uses the existing raw 3HFM PDB without local coordinate changes, chains H/L versus Y.
KY96A and RY73A are separate single mutants, not a double mutant. Their experimental
values in `../skempi_3hfm_ddg.csv` are +6.49 and -0.33 kcal/mol respectively.
MC samples = 20 (service default); seed = 42 (explicit reproducibility setting).
Tamarind automatically renumbers chains and reports the remapped mutation strings.
Check those strings and the output sign/units before interpretation. Checkpoint and
SKEMPI training overlap remain unverified; this is a diagnostic, not validation.

From repository root, commands actually run:

```bash
python3 episcaf_v3/energetics/tamarind/pilot.py prepare
python3 episcaf_v3/energetics/tamarind/pilot.py validate
python3 episcaf_v3/energetics/tamarind/pilot.py submit
python3 episcaf_v3/energetics/tamarind/pilot.py status
```

Retrieve later with `pilot.py status`, then `pilot.py download` when Complete.
The script refuses to repeat preparation or submission for this frozen run. An
interrupted submission must be resolved by checking status, not by deleting the guard.
`3hfm_pilot/` saves the input checksum, inline PDB request, authenticated schema,
validation response, normalized submitted request, submission receipt and status.
Credentials are read from `TAMARIND_API_KEY` or `~/.config/tamarind/api_key` and never
written into these artifacts. Download URLs are not persisted or printed.

Sources: https://app.tamarind.bio/api/tools/stabddg/schema (authenticated snapshot),
https://app.tamarind.bio/llms-full.txt and https://docs.tamarind.bio/tamarind/api.


## Reproduce the analysis locally (no API calls)

From the repository root:

```bash
python3 episcaf_v3/energetics/tamarind/pilot.py verify
python3 episcaf_v3/energetics/tamarind/report_pilot.py
MPLCONFIGDIR=/tmp/episcaf-mpl /usr/bin/python3 episcaf_v3/energetics/tamarind/plot_comparison.py
```

The plotting command uses this Mac's Python with pandas/scipy/matplotlib; another Python
with those dependencies works too. The report uses only the Python standard library.
The downloaded archive remains intact; the report checks input checksum, mutation mapping,
and returned settings, saves original small CSV/log files, writes comparison.csv and
summary.json, and regenerates the manuscript table. Plotting saves PNG/PDF plus exact
plot_points.csv and plot_statistics.json. No ML correlation is calculated for two points.
Sign/units source: https://github.com/LDeng0205/StaB-ddG (README, accessed 2026-09-08).
The log names model_ckpts/stabddg.pt, without a checkpoint hash. Upstream describes SKEMPI
fine-tuning; overlap with this deployed model is unresolved.

## Repeat the remote calculation without overwriting this run

This incurs a new Tamarind job. Choose a unique name and use the SAME new directory in
all four commands. Example from repository root:

```bash
python3 episcaf_v3/energetics/tamarind/pilot.py prepare --run-dir /tmp/3hfm-repeat --job-name episcaf-v3-3hfm-stabddg-repeat-UNIQUE
python3 episcaf_v3/energetics/tamarind/pilot.py validate --run-dir /tmp/3hfm-repeat
python3 episcaf_v3/energetics/tamarind/pilot.py submit --run-dir /tmp/3hfm-repeat
python3 episcaf_v3/energetics/tamarind/pilot.py status --run-dir /tmp/3hfm-repeat
```

Hosted model/software changes can alter predictions even with the same seed. The saved
schema/request/outputs document this run; they do not freeze Tamarind's execution image.


## Full measured antigen scan (2026-09-08)

After the two-point API check, the user requested all measured antigen mutations.
One new job contains all 13 as separate single mutants, using the same PDB, chain
partners, mcSamples=20 and seed=42. Native identities were checked for every mutation.
The initial pilot is preserved. Commands from repository root:

```bash
python3 episcaf_v3/energetics/tamarind/pilot.py prepare --all-measured --run-dir episcaf_v3/energetics/tamarind/3hfm_all13 --job-name episcaf-v3-3hfm-stabddg-all13-20260908
python3 episcaf_v3/energetics/tamarind/pilot.py validate --run-dir episcaf_v3/energetics/tamarind/3hfm_all13
python3 episcaf_v3/energetics/tamarind/pilot.py submit --run-dir episcaf_v3/energetics/tamarind/3hfm_all13
python3 episcaf_v3/energetics/tamarind/pilot.py status --run-dir episcaf_v3/energetics/tamarind/3hfm_all13
```

Retrieval and local analysis:

```bash
python3 episcaf_v3/energetics/tamarind/pilot.py download --run-dir episcaf_v3/energetics/tamarind/3hfm_all13
python3 episcaf_v3/energetics/tamarind/report_pilot.py --run-dir episcaf_v3/energetics/tamarind/3hfm_all13 --table-name tamarind_all13_table.tex
MPLCONFIGDIR=/tmp/episcaf-mpl /usr/bin/python3 episcaf_v3/energetics/tamarind/plot_comparison.py --run-dir episcaf_v3/energetics/tamarind/3hfm_all13
```

The last command updates the comparison PNG/PDF with all 13 ML predictions. The pilot
comparison.csv and result archive remain unchanged. Repeated mutations may change across
batch compositions even with the same seed; do not silently combine the two runs.


Full-run outcome: Complete (created 20:16:21, started 20:18:57, completed 20:19:41,
2026-09-08 API timestamps). Execution 44 s; queue 156 s; 0.01 weighted hours.
StaB-ddG Spearman 0.74553, Pearson 0.84542, MAE 1.17412, RMSE 1.87421 kcal/mol.
See comparison.csv for unrounded predictions and plot_statistics.json for metrics.
These are same-system descriptive results; SKEMPI training overlap is unresolved.
