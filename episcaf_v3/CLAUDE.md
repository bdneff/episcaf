# CLAUDE.md — episcaf_v3 operating contract (the energetic-epitope campaign)

Loaded when working under `episcaf_v3/`. This scopes *how* to work on the v3 campaign. The shipped
**v2** work lives at the repository root (its own `CLAUDE.md`, the DP4 library, the production
pipeline); v2 and v3 may run as **two separate campaigns for two separate agents**, so this file is
meant to stand on its own for the v3 one. For *what* v3 is, read `README.md`, the manuscript
(`manuscript/main.pdf`), and `docs/DECISIONS.md` — in that order — before acting.

## The project, in one paragraph
Energetic epitope design. v2 defines the epitope by 3D **geometry** (every antibody contact) and
fixes all of it in the RFdiffusion contig. v3 defines it by **energy** — which residues carry the
antigen–antibody binding free energy — and lets the design constraints follow: **Category 1**
(energetic / hot-spot) residues keep identity *and* shape; **Category 2** (structural / occluding)
residues keep shape but not identity; **Category 3** (the rest of the scaffold) is free, subject to
not clashing with or competing for the binding site. A cheap RFdiffusion **backbone pre-filter**
(epitope RMSD to native + antibody clash) drops doomed designs before ProteinMPNN/AlphaFold3. The
whole thing rests on getting the per-residue energy right, including solvent — so it is validated
against experimental ΔΔG before it drives any design. PI: John Altin.

## How to work here (the contract)
- **One decision at a time, logged.** Every methodological choice — the energy method, the MD
  protocol, the classification thresholds, the contig strategy, the filters — is a dated entry in
  `docs/DECISIONS.md`: options considered, the choice, why, and the reproducible check behind it.
  **Nothing is a "standard" in v3 until it is a logged entry with provenance.** This is the core
  discipline of the campaign.
- **Reproducible, informed, step by step.** Every number and figure regenerable by a named script
  from written-down inputs. Verify against data; never assert from memory or plausibility. Mark
  anything unverified `[UNVERIFIED]`. Prefer the correct, checkable answer to the quick one.
- **Reference, don't inherit.** The sibling `bcell_epitope` project is a **template** for the MD +
  per-residue interaction-energy machinery (its `campaign2_holo_epitopeness` computes almost exactly
  the interface decomposition we need), *not* a protocol to copy. Decide v3's energy definition,
  filters, and design strategy deliberately and record *why* they match or diverge from it.
- **Freeze inputs per run (MD).** Each simulation owns its frozen `.mdp`/config/scripts; never run
  against a moving shared file (adopted from `bcell_epitope`'s hard rule).
- **The manuscript is the living record.** When a decision or result changes, update `manuscript/`
  the same session (`tectonic main.tex`, not latexmk). It reads top to bottom.
- **Save scripts, don't paste them.** Inspection/analysis code that informs a decision becomes a
  committed named script (e.g. `energetics/skempi_overlap.py`), not a throwaway.

## Boundaries (carried from the repo; do not cross)
- **Do not ssh into or drive the cluster.** The user runs the GROMACS / RFdiffusion / AlphaFold3
  jobs. Stage inputs + frozen configs here, commit, and hand off; never submit or interface with
  SLURM yourself.
- Never `rsync --delete` toward `/tgen_labs`; never `git init` inside a data directory.
- Manuscript voice: flow, not staccato; **no "John said/wanted"** in the PDF — project voice.
- Tuned constants need recorded provenance **and** a reproducible check (the v2 cylinder-geometry
  lesson) — this is what `docs/DECISIONS.md` enforces.

## Where things live
- `manuscript/` — the living record (`main.tex` → `main.pdf`).
- `docs/DECISIONS.md` — the decision log (spine D1–D6).
- `energetics/` — MD + per-residue energy decomposition + Cat 1/2/3 labels **(to build)**;
  `skempi_overlap.py` is the D1 pilot-selection evidence.
- `contigs/` — the three-tier contig builder **(to build)**.
- `filters/` — the RFdiffusion backbone pre-filter **(to build)**.
- Reuse the **v2 method at the repository root** (`..`) by relative path — structures, scoring,
  AlphaFold3 stages — rather than forking it. Data and big outputs live on the cluster
  (`$WS` / `/tgen_labs`), never in git; `/scratch` trajectories are ephemeral.

## Current state (update as it moves)
- **v2 campaign:** DP4 library ordered; waiting on the assay binding data. The plan for reading it
  is at the repo root: `docs/DP4_RESULTS_ANALYSIS.md` + `scripts/build_dp4_binding_join.py` +
  `scripts/analyze_dp4_binding.py` (turnkey, self-checked).
- **v3 campaign:** D1 (pilot = **3HFM**) and D3 (MD protocol, from bcell) decided; the 3HFM holo MD
  is staged at `energetics/md/3hfm/` for the user to run on Gemini. Next: D2 (the per-residue energy
  quantity, once there's a trajectory). D6/D7 (the post-RFD3 filter + escalating-scale generation)
  have empirical backing from John's throughput pilots — see `docs/DECISIONS.md`. Mostly still to build.

## Commands that verify
- `cd manuscript && tectonic main.tex`                       # build the living record
- `python energetics/skempi_overlap.py --skempi <csv>`       # reproduce the D1 evidence

## Fixed facts
- The three categories: **Cat 1** energetic (fix identity + shape) · **Cat 2** structural/O-ring
  (fix shape, free identity) · **Cat 3** scaffold (free, must not clash with or compete for the site).
- bcell MD stack (**reference, not adopted**): AMBER99SB-ILDN + TIP3P, cubic box 1.2 nm, 0.15 M NaCl,
  PME, V-rescale + Parrinello-Rahman, 2 fs / LINCS h-bonds, ~20 ns holo production, GROMACS on GPU.
- Validate before design: the per-residue energy must recover known hot spots (experimental ΔΔG /
  SKEMPI) before the definition is allowed to drive any design.
