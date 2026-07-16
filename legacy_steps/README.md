# legacy_steps/ — one live dependency, kept here for that reason

This directory holds a single script, `05_rmsd_vs_af3.py`, carried over from the pre-reorg pipeline.
It is **not dead code**: `episcaf_pipeline/cli.py` (the `stage05` metrics step) shells out to it at
`cli.py:121`. Leave it here until that call is ported into `episcaf_pipeline/` proper.

If you move or rename it, update `episcaf_pipeline/cli.py` in the same commit. See `docs/PIPELINE.md`
for where stage05 sits in the run order and `docs/MIGRATION.md` for the reorg history.
