# Handoff — picking up episcaf (any model)

This repo is heavily documented on purpose: the living record *is* the handoff, so a new
agent (or a new model) can pick up cold. This file is the **map** — it tells you where to
look, not everything there is to know. Read the pointers below in order; they are the source
of truth, and this file defers to them wherever they disagree.

## Read these first, in order
1. **`CLAUDE.md`** (repo root) — the operating contract: *how* to work here (one step at a
   time, reproducible by construction, verify against data, don't drive the cluster). Applies
   to every session regardless of model. `AGENTS.md` points here too.
2. **`README.md`** (root) + **`docs/PIPELINE.md`** — what the shipped **v2** method is and the
   authoritative end-to-end run order (RFdiffusion3 → ProteinMPNN → AlphaFold3 → scoring →
   PepSeq oligos).
3. **`docs/DP4_LIBRARY.md`** and **`docs/DP4_RESULTS_ANALYSIS.md`** — the current v2 deliverable
   (the DP4 PepSeq library, ordered) and the staged plan for reading the assay data when it
   returns. This is v2's open loop (see Current state).
4. **`episcaf_v3/README.md`** — what the **v3** avenue is and why (summary below).
5. **`episcaf_v3/docs/DECISIONS.md`** — v3's decision spine (D1–D7). Every methodological
   choice is a dated entry with its rationale and a reproducible check. **Start here to see
   exactly where v3 stands** — it is the fastest way to get current.
6. **`episcaf_v3/manuscript/main.pdf`** (source in `manuscript/sections/*.tex`) — v3's living
   manuscript. Read top to bottom to get up to speed on the science.
7. **`docs/WORKFLOW.md`** — the two-machine git/compute loop (below).

Deeper background notes (read-only, optional; written by the prior Claude sessions) live at
`~/.claude/projects/-Users-bneff-Desktop-projects-episcaf-episcaf-v2/memory/` — start at
`MEMORY.md` (the index). The in-repo docs above are the authority; the memory is color.

## What v3 is, and why
v2 defines an epitope by **geometry**: every residue the antibody contacts (a heavy-atom
distance cutoff), and it hands RFdiffusion a contig that fixes the identity and backbone of
all of them. But those contacts don't contribute equally — a few carry most of the binding
free energy and the rest mostly fill space. **v3 defines the epitope by energy instead**:
which residues actually carry the binding, with the design constraints following from that —
**Category 1** (energetic / hot-spot) residues keep identity *and* shape; **Category 2**
(structural / occluding) residues keep shape but not identity; **Category 3** (the rest) is
free. Fixing fewer residues should give RFdiffusion more room to fold a good scaffold while
still holding the residues that do the binding. v3 does **not** replace v2 yet — v2 is the
source of truth and the shipped method; v3 is an exploratory subproject that takes over only
if it proves out, and writes down why if it doesn't.

The whole thing rests on getting the per-residue energy right (including solvent), so it is
**validated against experimental ΔΔG before it drives any design**. That validation is the
current v3 work (D2).

## Current state (2026-09-08)
- **v2 / DP4:** library is ordered; waiting ~4 weeks (from mid-Aug order) on the assay binding
  data — expected roughly late Sept / early Oct 2026. The turnkey analysis harness is staged and
  self-checked (`docs/DP4_RESULTS_ANALYSIS.md`, `scripts/build_dp4_binding_join.py`,
  `scripts/analyze_dp4_binding.py`). Nothing to do until the data lands.
- **v3 / D2 (active):** validating the per-residue energy on the 3HFM pilot (HyHEL-10 Fab /
  hen egg lysozyme) against SKEMPI alanine-scan ΔΔG. Criterion is **precision, not recall**
  (Jacob's): no residue the method flags as high-energy may have a low measured ΔΔG; false
  negatives are fine (they're Category-2 structural residues). Three methods tried so far —
  bare interaction energy (over-calls charged residues, e.g. Arg73), LJ-only (best rank
  correlation but most false positives), and per-residue MM-GBSA (desolvation demotes Arg73
  correctly — mechanism confirmed — but 20-frame single-replica GB is too noisy: it
  over-desolvates the biggest hot spot Lys96 and ρ collapses to 0.05). Full write-up in
  `episcaf_v3/manuscript/sections/energetic_validation.tex` and `DECISIONS.md` D2.
  - **Open loop right now:** a **denser MM-GBSA resample** (interval 100→20, i.e. 20→100 frames
    of the same trajectory) is staged and being submitted on Gemini by Brandon. When it
    finishes he brings back `complex_tdc.txt` and `ligand_tdc.txt`; the next steps are then
    `energetics/mmgbsa_decomp_to_csv.py` → `energetics/plot_ie_vs_ddg.py --channel mmgbsa_dg`
    → update the manuscript table/figure and `DECISIONS.md` D2 → commit. The exact recipe is in
    `episcaf_v3/energetics/mmgbsa/README.md` §4 and `manuscript/figures/FIGURES.md`.
- **v3 / Tamarind pilot (2026-09-08):** completed StaB-ddG on separate K96A and R73A
  mutations: +1.763 and +0.481 kcal/mol versus experimental +6.49 and -0.33. Frozen
  requests/results and reproducible comparison are in `episcaf_v3/energetics/tamarind/`.
  See its README and manuscript learned-ddG subsection. Two-point diagnostic only;
  checkpoint hash and training overlap unresolved. No estimator or threshold adopted.
- **v3 still to build:** D4 (Cat 1/2/3 classification thresholds), D5 (the three-tier contig
  strategy — decided fresh for v3, *not* inherited), D6/D7 (RFD3 backbone pre-filter +
  escalating-scale generation; have empirical backing from John's throughput pilots). See the
  decision spine in `DECISIONS.md`.

## How to work here (the three hard rules; full contract in CLAUDE.md)
- **One step at a time, logged.** In v3 especially: nothing is a "standard" until it's a dated
  entry in `docs/DECISIONS.md` with a rationale and a reproducible check. One change, verified,
  then the next — don't batch unrelated work.
- **Reproducible by construction.** Every number/figure regenerable by a named, committed script
  from written-down inputs; record the exact command (`manuscript/figures/FIGURES.md`). Verify
  against real files; never assert from memory. Mark unverified things `[UNVERIFIED]`.
- **Do not drive the cluster.** Brandon runs all GROMACS / RFdiffusion / AlphaFold3 / MM-GBSA
  jobs on Gemini. You **stage** inputs + frozen configs, commit, push, and hand off the exact
  commands; you never ssh in or submit SLURM yourself.

## Collaborating with Brandon (clarified 2026-09-08)

Brandon Neff has a PhD in chemistry specializing in molecular dynamics; immunology is the newer
domain for him. Treat MD, statistical mechanics, sampling, and energetics as his technical
foundation, and explain immunology-specific assumptions and assay interpretation explicitly.
He sees energetic approaches as a contribution he brings to the TGen project. Relevant work:
[Fast Sampling of Protein Conformational Dynamics](https://arxiv.org/abs/2411.08154) and
[Protein-Water Energy Transfer via Anharmonic Low-Frequency Vibrations](https://arxiv.org/abs/2601.02699).
The v3 static/rigid starting approximation is deliberate scope control, not an assumption that
proteins actually lack dynamics. See the 2026-09-08 entry in `episcaf_v3/docs/DECISIONS.md`.

## The git/compute loop (see docs/WORKFLOW.md)
Local (here): edit, stage configs/scripts, commit, `git push`. Cluster (Gemini, Brandon runs):
`git pull`, run the job, bring results back. Data and large outputs never go in git — they live
on `/tgen_labs` (persistent) and `/scratch` (ephemeral); `$WS =
/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff` is the durable workspace. The repo and the
cluster join only through paths, not by moving data into git.

## Build / verify commands
- `cd episcaf_v3/manuscript && tectonic main.tex` — build v3's living record (not latexmk here).
- `cd manuscript && tectonic main.tex` — build the v2 manuscript.
- `python tests/test_scoring.py` — v2 scorer unit tests (no data needed).
- `python episcaf_v3/energetics/plot_ie_vs_ddg.py --ie <csv> --ddg episcaf_v3/energetics/skempi_3hfm_ddg.csv --channel <ab_total|ab_lj|mmgbsa_dg> --out <png>` — the D2 validation plot.
