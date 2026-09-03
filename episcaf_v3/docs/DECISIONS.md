# episcaf_v3 decision log

The record of every methodological choice in v3. Each is an entry below: the options
considered, what we chose, why, and the check that backs it. Nothing is a "standard" in v3 until it
appears here with a rationale and a check. This mirrors the lesson from v2 that a tuned constant
needs recorded provenance and a reproducible test behind it, applied from day one.

Format for each decision:

> **Dn — <title>** (YYYY-MM-DD)
> *Question:* …
> *Options:* …
> *Decision:* …
> *Why:* …
> *Check / provenance:* the command, dataset, or reference that backs it (and how to reproduce).
> *Status:* open | decided | revisit-if …

---

## The decision spine (order we intend to work through)

These are **open** until they have an entry. They are listed so the path is visible; listing is not
deciding.

- **D1 — Pilot complex + validation anchor.** *Status: **DECIDED** 2026-09-02 (see Decisions below).*
  **3HFM** (HyHEL-10 / hen egg lysozyme) is the pilot, **1VFB** (D1.3 / lysozyme) a second check. Its
  holo MD is staged at `energetics/md/3hfm/`.
- **D2 — What "energetic" means.** *Status: **IN PROGRESS** — method being fit/validated on 3HFM.*
  The per-residue quantity that ranks residues, and how solvent enters it. **Approach (empirical
  fit):** rather than argue the method a priori, fit it to experimental ΔΔG. Start with the
  cheap direct interaction-energy decomposition (`energetics/holo_ie.py`, adapted from
  `bcell_epitope`: reaction-field Coulomb + Lorentz–Berthelot LJ, per antigen residue vs. antibody,
  averaged over the trajectory), compare to SKEMPI's 3HFM alanine ΔΔG (`skempi_3hfm_ddg.csv`,
  `plot_ie_vs_ddg.py`), then test generalization on a few more complexes.
  **Success criterion (Jacob) — precision, not recall:** the method must have **no false positives**
  — no residue it flags as high interaction energy may have a *low measured* ΔΔG (residues with no
  measured value are exempt). **False negatives are acceptable**: a real hot spot the per-residue
  energy misses is expected — it is likely critical for structure/orientation/entropy, i.e. a
  **Category-2** residue, handled by shape-preservation, not identity-locking. So we tune the energy
  method (cutoff, terms, whether the solvent channel enters) to empty the false-positive zone.
  If direct interaction energy alone leaves false positives, that is the evidence we need the fuller
  desolvation/water-mediated treatment (MM-GBSA / computational alanine scan).
  **First result (2026-09-03, 20 ns single replica):** raw interaction energy (`ab_total`) FAILS the
  precision criterion. It recovers the top hot spots (K97 fav_ie 467, K96 126, both real; Y20, D101
  positive), but produces a clear false positive at **ARG73 — fav_ie 186 kJ/mol (2nd-strongest of all
  129 residues) yet experimental ΔΔG −0.33** — plus weaker charged false positives (R21). This is the
  predicted failure: raw Coulomb over-calls charged residues because it ignores the desolvation they
  pay. Spearman ρ=+0.46. **LJ-only probe:** dropping Coulomb (`--channel ab_lj`) demotes Arg73 and
  lifts ρ to 0.73 (top-3 by LJ all real hot spots), but then over-calls well-packed non-critical
  residues (Leu75) — which brackets the problem: electrostatics over-calls charged residues, packing
  over-calls buried ones, so the desolvation-aware *net* ΔG is the fix. **MM-GBSA setup staged**
  (`energetics/mmgbsa/`: config + recipe) — single-trajectory per-residue decomposition on the
  existing 3HFM trajectory; the test is whether the GB desolvation term empties the false-positive
  zone (demote Arg73 + Leu75/Arg21, keep K96/K97/Y20/D101).
- **D3 — MD protocol.** *Status: **DECIDED** 2026-09-02 (starting template; see Decisions below).*
  Adopt `bcell_epitope`'s frozen GROMACS stack as the starting protocol (AMBER99SB-ILDN + TIP3P, cubic
  box 1.2 nm, 0.15 M NaCl, PME, V-rescale + C-rescale(equil)/Parrinello-Rahman(prod), 2 fs / LINCS
  h-bonds), 20 ns holo production, frozen per run under `energetics/md/<pdb>/`; revisit if the
  energetics call for it. Staged for 3HFM.
- **D4 — Classification rule.** The threshold(s) that split Cat 1 (energetic) / Cat 2 (structural) /
  Cat 3 (rest). *Reference (not adopted):* `bcell_epitope` calls a residue "enthalpic" when its
  antibody interaction energy is ≤ −3σ of the repulsive tail — one concrete rule to weigh.
- **D5 — Three-tier contig strategy.** How each category maps to RFdiffusion (fix backbone / diffuse)
  and ProteinMPNN (lock identity to native / free), plus the randomization ranges (scaffold length,
  epitope register in the construct, native-spacing gaps for the contact shell). **This is the part
  we explicitly decide fresh for v3, not inherit.**
- **D6 — RFdiffusion backbone pre-filter.** The epitope-RMSD-to-native cutoff and the antibody-clash
  test applied at the backbone stage, before ProteinMPNN. Candidate to reuse v2's cylinder/clash
  geometry rather than define new. *Empirical backing (John, 2026-09; see Reference material):* the
  clash count on the RFD3 backbone (which also comes with a draft sequence) is highly predictive of
  final post-AlphaFold3 success, so filtering there skips AF3 — the dominant cost — on doomed backbones.
- **D7 — Escalating-scale generation with a stopping criterion (John's second idea).** Run RFdiffusion
  at escalating scales per target (e.g. 100 → 1,000 → 10,000) with a composite-score stopping rule, so
  compute concentrates on the hard cases and stops early on the easy ones. Open.
- **Cross-cutting — platform-agnostic pipeline.** The workflow should run on Gemini *or* Tamarind
  (John is piloting on Tamarind, agnostic between them), and possibly directly on AWS, chosen by
  throughput/economics. Keep the stages portable across schedulers, not wired to one. Open.

---

## Reference material on file (inputs to the decisions above, not decisions themselves)

- **`bcell_epitope` MD + energy stack** — the sibling project's `campaign2_holo_epitopeness` computes
  a per-residue interface interaction-energy decomposition (reaction-field Coulomb + Lorentz–Berthelot
  LJ, vs. antibody and vs. solvent channels) plus a single-water-bridge channel, from ~20 ns GROMACS
  trajectories; reusable scripts include `analysis/scripts/holo_ie.py`, `holo_ie_bridge.py`,
  `holo_mindist.py`, and the frozen mdp stack under `campaign1_apo_features/md/1AKI/apo/configs/`.
  Its "epitopeness" is a persistence-weighted interaction energy, explicitly **not** ΔΔG_bind — a key
  thing D2 must decide whether to keep or extend. Treated as template and reference throughout.

---

- **SKEMPI 2.0 overlap with our structures (D1 evidence, 2026-09-02).** Ran
  `energetics/skempi_overlap.py` against SKEMPI 2.0 (7,085 mutation records). **None of our 59
  scaffolded structures appear in SKEMPI** — our set is almost entirely recent (2020s) PDBs and
  SKEMPI 2.0 is a 2019 database, so we cannot validate on one of our own complexes. The fallback is
  the standard antibody:antigen (AB/AG) alanine-scan complexes in SKEMPI. Top candidates by number of
  alanine measurements and hot-spot content:
  - **3HFM** — HyHEL-10 / hen egg lysozyme: 18 Ala, ΔΔG up to +6.5 kcal/mol, 7 hot spots. Textbook
    system, **and lysozyme is the system `bcell_epitope` already runs MD on** → reuses validated
    machinery. *Recommended primary pilot.*
  - **1VFB** — D1.3 / lysozyme: 16 Ala; another classic lysozyme benchmark. *Recommended second check.*
  - Others: 1DVF (D1.3/E5.2), 1N8Z (Herceptin/HER2), 3NGB (VRC01/gp120, HIV), 1JRH (32 Ala, largest).
  Reproduce: `curl -sSL -o skempi_v2.csv https://life.bsc.es/pid/skempi2/database/download/skempi_v2.csv`
  then `python energetics/skempi_overlap.py --skempi skempi_v2.csv`.

- **John's throughput pilots (2026-09; PDF pending his GitHub access).** Relayed via Slack, summarized
  from a PDF not yet in the repo. Two findings: (1) the RFD3 backbone clash count is highly predictive
  of post-AlphaFold3 success, so a post-RFD3 / pre-MPNN clash filter avoids running AF3 (the dominant
  cost) on doomed backbones; (2) escalating-scale RFD3 (100 / 1,000 / 10,000) with a composite-score
  stopping rule concentrates compute on hard cases. Plus: stay platform-agnostic (Gemini / Tamarind /
  AWS). Empirical basis for D6 and D7 — fold in the exact numbers/plots once the PDF is in the repo.

---

## Decisions (most recent first)

### D1 — Pilot complex + validation anchor (2026-09-02)
- *Question:* which complex do we validate the per-residue energy method on before it drives design?
- *Options:* one of our own 59 structures (none are in SKEMPI); a standard SKEMPI antibody:antigen
  alanine-scan complex.
- *Decision:* **3HFM** (HyHEL-10 / hen egg lysozyme) as the primary pilot; **1VFB** (D1.3 / lysozyme)
  as a second check.
- *Why:* 3HFM has 18 experimental alanine ΔΔG in SKEMPI (hot spots up to +6.5 kcal/mol), it is a
  textbook hot-spot system, and lysozyme is already the system `bcell_epitope` runs MD on — so we
  reuse validated machinery. None of our own structures overlap SKEMPI (2020s PDBs vs. a 2019 database).
- *Check / provenance:* `energetics/skempi_overlap.py` reproduces the overlap and the candidate list.
- *Status:* decided.

### D3 — MD protocol (starting template) (2026-09-02)
- *Question:* what GROMACS protocol produces the trajectory the energy decomposition runs on?
- *Options:* reconstruct fresh; adopt `bcell_epitope`'s proven stack.
- *Decision:* adopt bcell's frozen stack as the *starting* protocol — AMBER99SB-ILDN + TIP3P, cubic
  box 1.2 nm, 0.15 M NaCl, PME, V-rescale + C-rescale (equil) / Parrinello-Rahman (prod), 2 fs / LINCS
  h-bonds — with 20 ns holo production, frozen per run under `energetics/md/<pdb>/`.
- *Why:* a validated, documented protocol from the same lab, on lysozyme systems; no reason to
  reinvent it for the pilot. Marked a starting protocol, to revisit if the energetics require it
  (e.g. longer runs for the interaction-energy averages to converge).
- *Check / provenance:* `energetics/md/3hfm/` — `configs/*.mdp` copied verbatim from bcell
  (`campaign1_apo_features/md/1AKI/apo/configs/`); the run README documents the one holo departure
  (Fab H–L interchain disulfide / `pdb2gmx -merge`).
- *Status:* decided (starting template); revisit-if convergence or the energetics require.
