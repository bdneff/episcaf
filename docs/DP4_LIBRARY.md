# DP4 library — components, selection, and status

Living reference for the DP4 PepSeq library: what each component is, how designs are selected, and
where each deliverable lives. Manuscript counterpart: `manuscript/sections/dp4_library.tex`
(`sec:dp4`). Related: `docs/CASE_ENCODING.md`, `docs/CYLINDER_PARAMS.md`.

Status: all six components selected/built. **Only the ranked components (C1/C2/C3, marked †) have a
selectable depth** — their counts are shown at **top-20 per group** and scale with the chosen shipping
depth (top-*n*, set at assembly from the budget). C4/C5 are fixed-size; C6 (‡) is derived from C1's
top-20 base, so it scales with C1's depth.

## Components

| Comp | Component | What it is | Selection | Constructs | Deliverable | Status |
|---|---|---|---|---|---|---|
| **C1** | known-Ab, whole epitope | best-*n* scaffolds per mAb (comparator) | ranked, top-*n* per mAb | 1,180 †§ | `results/dp4_C1_whole_epitope_ranked.top20.csv` | ranked |
| **C2** | known-Ab, single island | best-*n* per (mAb, island); 87 island contigs | ranked, top-*n* per island | 1,740 † | `results/dp4_C2_single_island_ranked.top20.csv` | ranked |
| **C3** | polyclonal 12-mer tiles | best-*n* per window, no antibody | ranked, top-*n* per window | 8,780 † | `results/dp4_C3_12mer_ranked.top20.csv` | ranked |
| **C4** | linear 30-mer controls | bare tiled peptides (no scaffold) | exhaustive tiling (fixed) | 2,174 | `data/libraries/dp4_tiled30mers_fasta.csv` | built |
| **C5** | metric-space titration | designs spread across metrics (calibration) | farthest-point sample (fixed) | 3,000 §  | `results/dp4_C5_titration.csv` | sampled |
| **C6** | scaffolded-epitope controls | island→Ala + scaffold-disruption | all C1 top-*n* base × flavors | 3,071 ‡§ | `results/dp4_C6_controls.csv` | built |

† **Ranked selection** — count shown at **top-20 per group**; scales with the shipped depth (top-*n*,
elastic — see *Budget & depth*). At top-5 these shrink ~4×.
‡ **Derived** from C1's top-20 base, so scales with C1's depth (island1→A + island2→A dual-only + scaffoldMutX4).
§ **104-mer designs** — C1 reused Lawson's 104-residue contigs, and C5/C6 derive from the C1 pool, so
these three are 104-mers **trimmed to the 103-mer assay ceiling at assembly** (epitope-preserving; see
*104→103 truncation*). C2/C4 are natively 103; C3's length is not yet confirmed.

Full composite rankings (`results/dp4_*_ranked.csv`) are regenerable and gitignored; only the top-*n*
cuts + case-encoded sequences are tracked.

## How "best 20" is defined (C1 / C2 / C3)

We **rank, we don't gate** — no hard pass/fail cutoff. Every design gets one composite score; we rank
within each group and take the top *n* per group. Tooling: `scripts/stage06_select.py`; dials in
`episcaf_analysis/presets.py`; scorer `episcaf_analysis/score.py`.

**The math.** Each of four metrics is converted to a **percentile** within its population and oriented
so higher = better (all four are lower-is-better), then weighted-summed (weights sum to 1):

```
composite = 0.35 · accessibility
          + 0.35 · epitope_RMSD
          + 0.15 · overall_RMSD
          + 0.15 · epitope_PAE
```

- **Accessibility term** depends on what's known: **real AF3 clash** (`af3_n_clash_res`) for C1/C2
  (antibody known); **native-aware cylinder** surrogate for C3/C5 (no antibody).
- **Percentile scope:** pooled across the mAb set (C1/C2); per-antigen for the 12-mer set (C3).
- **Grouping for "top *n*":** per mAb / `id` (C1); per `(id, island_index)` (C2); per `(antigen, id)` (C3).

**Weight provenance.** Set from the DP3 binding data: accessibility and epitope-RMSD were the
strongest within-antibody predictors of experimental enrichment (~0.35 each); overall-RMSD and PAE
carried little binding signal (0.15 each). This is a hand-set prior from an all-passing set — **C5 is
designed to span the metric space so these weights can be re-fit on real DP4 binding** (manuscript Q2).

## Exclusions — the canonical 56-mAb set (John, DP4)

Apply **one** exclusion set across the known-Ab components so counts are consistent. From the 59 DP3
known-Ab epitopes, drop three → **56**:
- `4xwo_5P` — low assay yield
- `7a3t_0P` — 4-residue epitope (smallest in DP3, too small to present)
- `2h32_0P` — "not a standard antibody case" (John) — *confirm the exact reason and record it here*

Applies to the **known-Ab components: C1, C2, C5, C6**. **C3** (polyclonal 1D2K/6M0J/4WAT) is
unaffected (different antigens). **C4** (linear tiled controls) is an **open question**: drop these
three antigens too, or keep them as linear controls? — decide with John.

Tooling supports it: `stage06_select.py --drop-ids 2h32,4xwo,7a3t`, `stage06_sample_c5.py` `DROP_IDS`
(updated to include 2h32), `build_c6_mutants.py --drop-targets 2h32,4xwo,7a3t`. The currently-committed
deliverables predate the `2h32` drop (C1/C2/C3 are 59-set rankings; C5/C6 dropped only `4xwo`/`7a3t` = 57),
so **the canonical 56-set is applied uniformly at the final assembly cut** (together with the chosen
depth) — no premature re-runs.

## Case-encoded `designedSequence` for visualization (John, DP4)

John asked that outputs carry the design sequence with **epitope UPPERCASE, scaffold lowercase** — which
is exactly the `scaffoldEPITOPE` column we already produce (`docs/CASE_ENCODING.md`). Status: **C1 and C5
done** (`results/dp4_C{1,5}_scaffoldEPITOPE.csv`); **C2 and C3 need their case-encode run** (same
mechanism; C3 is a different run). In the assembled DP2 file this casing is carried as the
`designedSequence` column, so John's visualization ask is satisfied by the assembled output.

## Component notes

- **C5 — metric-space titration** (`scripts/stage06_sample_c5.py`). Farthest-point (max–min) spread
  over the four standardized scoring axes (cylinder as accessibility), ~53 per mAb over the 57 mAbs.
  3,000 designs spanning 89–97% of each axis's full range — deliberately including the low-quality tail
  the filters reject, so binding read off the spread calibrates the scorer.
- **C6 — scaffolded-epitope controls** (`episcaf_pipeline/scaffolded_epitope_controls/`). Base = C1
  top-20 over **57 mAbs** (dropped `4xwo` low-yield, `7a3t` 4-residue epitope). Not new scaffolding —
  string substitution on the case-encoded sequence (port of John's DP3 R code). Flavors: every-residue
  island1→Ala, island2→Ala (dual-island only), and `scaffoldMutX4` (a `PPDDGG` hexamer in 4 scaffold
  windows, each ≥4 residues from the epitope, seeded). **89/1,140 (7.8%)** designs can't fit 4 hexamers
  → X4 arm covers 1,051/1,140; alanine arms cover all.
- **C4 — linear tiled-30mer controls** (`data/libraries/dp4_tiled30mers_fasta.csv`, John-approved).
  The **full antigen sequences** of the ~60 proteins (57 mAb targets + 3 tiled antigens 1D2K/6M0J/4WAT),
  taken from the **FASTA files (no unresolved-gap holes, unlike the PDBs)**, chopped into overlapping
  **30-mers, step = 6**. No RFD/MPNN/AF3. Each 30-mer goes at the END of a constant construct:
  `GSGAGSGA…GSGA` filler + `ENLYFQGA` (TEV protease site) + `[30-mer]` → a constant **103-mer** that is
  cleaved to the 30-mer in the final step (matches the linear-epitope assays). Use the `_fasta` files,
  not the earlier PDB-derived ones.
- **Case-encoding** (`scripts/case_encode_selected.py`, `docs/CASE_ENCODING.md`). C1/C5 designs' epitope
  positions were recovered (token → `dp2.assay_scaffolded_epitope_id` → `scaffolded_epitope_chunk_resindices`)
  and written as case-encoded sequences (`results/dp4_C{1,5}_scaffoldEPITOPE.csv`), feeding C6 + assembly.

## Assembly format (DP2 annotated, John-approved)

The final synthesis file is the **8-column DP2 annotated format** — already settled and approved on C4
(`dp4_tiled30mers_fasta.csv` is the reference instance):

| column | meaning |
|---|---|
| `library_member` | global id, `DP4_<N>` |
| `sequence` | the constant **103-mer** actually synthesized |
| `category` | component type (`tiled30mer`, `scaffoldedAbEpitope`, `metricSpaceTitration`, …) |
| `model` | `RFD` for designs, `(none)` for linear controls |
| `designedSequence` | the payload (the 30-mer for C4; the scaffold for C1/C2/C3/C5/C6) |
| `designedSequenceLength` | len(`designedSequence`) |
| `design_ID` | within-category id |
| `target` | antigen / mAb id |

Per-category `sequence` construction: **linear controls (C4)** = filler + `ENLYFQGA` + 30-mer;
**scaffolded designs (C1/C2/C3/C5/C6)** = the design's own 103-mer directly. Assembly concatenates all
components in this schema with global `library_member` numbering.

**104→103 truncation — the whole-epitope (C1) family only.** Not all components are 104-mers. **C1
*reproduced* Lawson's whole-epitope run, reusing his contigs** (`contig_length "104-104"`), so **C1 —
and C5 (sampled from C1's pool) and C6 (built from C1) — are 104-residue** proteins (`af3_window_end=104`)
that must be trimmed to the 103-mer assay ceiling. **C2 is natively 103** (we *generated* new contigs
at 103, `build_dual_island_designs.py`, correcting Lawson's 104→103) — no trim. **C3 (12-mer): length
not yet confirmed** (manuscript approximates ~104; check the 12-mer run's contig length before assuming).
Trim rule for the 104 family: drop one residue from whichever terminus is **scaffold** (default
C-terminal; N-terminal when the C-terminus is an epitope residue), so no epitope residue is lost. In C1,
2 designs (`6qb6_0P`) have an epitope C-terminus and are handled by the N-terminal trim; C5 had none.
The `3ux9_1P` epitope-at-both-termini case is dropped for its **rank-21** replacement (kept at 20/epitope).

*Why the trim is safe (structural rationale, recorded so we remember it).* A single **terminal** residue
is the least structurally committed position in a fold — termini fray, make few tertiary contacts, and
don't anchor the hydrophobic core — so its contribution to stability is small, especially for hyperstable
designed miniproteins. We cut the **scaffold** terminus, distal from the **preserved** epitope, so the
epitope's environment/geometry (set by the intact core) is untouched. Residual risk is limited to the rare
design whose trimmed residue caps a helix or forms a specific contact — a small *local* destabilization,
not fold-breaking. NOTE the DP3 binding data was measured on the **104**mers, so the trim is reasoned from
folding principles, **not re-validated**; a spot re-fold of a 103mer sample (confirm epitope RMSD unchanged)
would close that gap cheaply. Full detail in manuscript `sec:production`.
**One exception that cannot be trimmed either way:** `3ux9_1P` **rank 9** (token
`0ab98e18e5c6a6de6dc3f9a25881ee10`, mpnn_id 2) — its 24-residue epitope reaches **both** termini, so any
trim to 103 clips one epitope residue. Decision: **drop it and ship the rank-21 design for `3ux9_1P`**
(keeps 20/epitope; the full ranking is regenerable, so rank-21 is well-defined). Apply this at assembly.

## Reproduce (exact commands)

Every deliverable is regenerable by a named script; the numbers reported (counts, coverage, the 89 X4
skips) are printed by these scripts, not hand-computed. Metric CSVs live in the local sibling data dirs
(`$D = /Users/bneff/Desktop/projects/episcaf`, see `filesystem-map`); C2 + case-encoding run on Gemini.

```bash
# C1 (local)
python scripts/stage06_select.py --preset antibody \
  --metrics-csv $D/known_antigen/analysis/data/metrics_native_cyl_full.csv \
  --group id --topk 20 --out results/dp4_C1_whole_epitope_ranked.csv

# C2 (Gemini, at the dual-island run's metrics)
python scripts/stage06_select.py --preset antibody \
  --metrics-csv runs/dual_island_rfd3/05_analysis/metrics_dual_island.parquet \
  --group id,island_index --topk 20 --out results/dp4_C2_single_island_ranked.csv

# C3 (local)
python scripts/stage06_select.py --preset twelvemer \
  --metrics-csv $D/12mer_tiling/analysis/data/metrics_12mer.csv \
  --group antigen,id --topk 20 --out results/dp4_C3_12mer_ranked.csv

# C4 (local; defaults -> data/libraries/dp4_tiled30mers_fasta.csv)
python -m episcaf_pipeline.build_dp4_tiled30mers_fasta

# C5 (local; deterministic FPS)
python scripts/stage06_sample_c5.py \
  --metrics-csv $D/known_antigen/analysis/data/metrics_native_cyl_full.csv \
  --total 3000 --out results/dp4_C5_titration.csv

# Case-encode C1 + C5 (Gemini SLURM) -> results/dp4_C{1,5}_scaffoldEPITOPE.csv
sbatch scripts/case_encode_selected.sbatch

# C6 (local; seeded, reproducible)
python episcaf_pipeline/scaffolded_epitope_controls/build_c6_mutants.py \
  --input results/dp4_C1_scaffoldEPITOPE.csv \
  --id-col token --target-col target --seq-col scaffoldEPITOPE \
  --drop-targets 4xwo,7a3t --out results/dp4_C6_controls.csv
```

Scorer weights/transforms are config, not magic numbers: `episcaf_analysis/presets.py` (provenance
above). C5 and C6 are deterministic (FPS is seed-free deterministic; C6 seeds its RNG).

## Budget & depth

DP4 = a 36k library that includes all minibinders → **~10–15k slots for Episcaf designs**
(`memory: dp4-budget`). Fixed core (C1–C4, C6) is built; **C5 + the selection depth (top-*n*) are the
elastic buffer**, sized last once the final minibinder count is known. At top-5 the ranked counts
shrink ~4× from the top-20 figures above.

## Pending — assembly & encoding

Format is settled (above), so assembly is **not** blocked on a template anymore — only on John's okay
on the components and the chosen depth.

1. **Assembly (`06_library`)** — build each component's rows in the 8-column DP2 schema (C4 already is),
   concatenate, assign global `library_member` numbering. Needs: John's **okay on components/counts**
   and the **shipping depth** (top-*n*, sized from the budget). Then it's a mechanical build.
2. **Oligo encoding** — LadnerLab `oligo_encoding` + DP3 codon weights
   (`episcaf_pipeline/oligo_encoding/`), then the order-file step (confirm with Erin).
