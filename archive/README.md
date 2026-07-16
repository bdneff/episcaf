# archive/ — superseded code, kept on purpose

This holds earlier versions of scripts and pipeline steps that have been replaced. Nothing here is
on a live code path; it is kept so the history of how the pipeline was built stays legible (and so a
result produced by an old script can still be traced). The reorganization that created it is recorded
in `docs/MIGRATION.md`.

Do not import from here. For the current pipeline see `docs/PIPELINE.md`; for what replaced a given
file, `docs/MIGRATION.md` maps old to new. Everything here is also in git history, so it can be pruned
later without losing anything.

Subdirs: `legacy_sbatch/`, `legacy_steps/`, `legacy_tools/`, `scripts/`, `top_level/`. The `.bak`
and dated `.bak_YYYYMMDD_*` files are point-in-time backups from the migration.
