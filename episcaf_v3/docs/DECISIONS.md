# episcaf_v3 decision log

The reproducibility backbone of v3. Every methodological choice is an entry below: the options
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

- **D1 — Pilot complex + validation anchor.** Pick 1–2 antibody:antigen complexes with experimental
  hot-spot / ΔΔG data (SKEMPI 2.0) so every later choice is testable against known numbers before it
  drives any design. *Status: evidence gathered (see below), a pilot proposed, awaiting confirmation.*
  *Proposed:* **3HFM** (HyHEL-10 / hen egg lysozyme) as the primary pilot, **1VFB** (D1.3 / lysozyme)
  as a second lysozyme check — see the SKEMPI finding under Reference material.
- **D2 — What "energetic" means.** The per-residue quantity that ranks residues, and how solvent
  enters it (a direct interaction-energy decomposition à la `bcell_epitope`, vs. a fuller
  ΔΔG-of-binding that includes desolvation/entropy; implicit solvent vs. an explicit water-mediated
  channel). This is the definitional core of v3.
- **D3 — MD protocol.** Force field, water model, box, ions, equilibration, sampling and run length.
  *Reference (not adopted):* `bcell_epitope` runs AMBER99SB-ILDN + TIP3P, cubic box 1.2 nm clearance,
  0.15 M NaCl, PME, V-rescale + Parrinello-Rahman, 2 fs / LINCS h-bonds, 20 ns holo production on
  GROMACS (a100/v100), frozen mdp per run. We decide v3's protocol against this template, with
  reasons for any deviation.
- **D4 — Classification rule.** The threshold(s) that split Cat 1 (energetic) / Cat 2 (structural) /
  Cat 3 (rest). *Reference (not adopted):* `bcell_epitope` calls a residue "enthalpic" when its
  antibody interaction energy is ≤ −3σ of the repulsive tail — one concrete rule to weigh.
- **D5 — Three-tier contig strategy.** How each category maps to RFdiffusion (fix backbone / diffuse)
  and ProteinMPNN (lock identity to native / free), plus the randomization ranges (scaffold length,
  epitope register in the construct, native-spacing gaps for the contact shell). **This is the part
  we explicitly decide fresh for v3, not inherit.**
- **D6 — RFdiffusion backbone pre-filter.** The epitope-RMSD-to-native cutoff and the antibody-clash
  test applied at the backbone stage, before ProteinMPNN. Candidate to reuse v2's cylinder/clash
  geometry rather than define new.

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

---

## Decisions (most recent first)

*(none yet — D1 is proposed above and awaiting confirmation)*
