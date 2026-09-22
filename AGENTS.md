# AGENTS.md

Operating contract for any coding agent in this repo (Codex and others). This project is
model-agnostic: the same rules and the same living record apply whoever is driving.

**The contract lives in `CLAUDE.md` (repo root). Read it; it governs how to work here.** When
working under `episcaf_v3/`, also read `episcaf_v3/CLAUDE.md` (the v3-specific contract). This
file does not restate them; it points to them so they don't drift.

**To get oriented, read `docs/HANDOFF.md`.** It maps what this project is, what "v3" is and
why, where everything lives, and the current open work. Start there, then follow its pointers.

The three rules that catch people out:
1. **One step at a time, logged.** In v3, nothing is a standard until it's a dated entry in
   `episcaf_v3/docs/DECISIONS.md` with a rationale and a reproducible check.
2. **Reproducible by construction.** Every number/figure regenerable by a named committed
   script; verify against real files, never assert from memory.
3. **Do not drive the cluster.** The user runs all Gemini/SLURM jobs (GROMACS, RFdiffusion,
   AlphaFold3, MM-GBSA). You stage inputs + configs, commit, push, and hand off exact commands.
