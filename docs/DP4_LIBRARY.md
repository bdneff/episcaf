# DP4 library — components, selection, and status

Living reference for the DP4 PepSeq library: what each component is, how designs are selected, and
where each deliverable lives. Manuscript counterpart: `manuscript/sections/dp4_library.tex`
(`sec:dp4`). Related: `docs/CASE_ENCODING.md`, `docs/CYLINDER_PARAMS.md`.

**Status: ASSEMBLED AND SHIPPED (2026-07-13).** All seven components (C1–C6 + the 8VDL arm) are
selected/built, case-encoded, and concatenated into the final synthesis file
**`data/libraries/dp4_library.csv` — 12,251 constructs**, every one a 103-mer with a unique
`library_member` and `design_ID`. The shipping depth is **top-20 per group for C1/C2** and **top-3 per
window for C3**; C4/C5/C6 and 8VDL are fixed-size. Only **oligo encoding** remains (see *Pending*).

## Components (as shipped)

| Comp | Component | What it is | Selection | Constructs | Deliverable |
|---|---|---|---|---|---|
| **C1** | known-Ab, whole epitope | best scaffolds per mAb (comparator) | ranked, **top-20** per mAb | **1,120** § | `results/dp4_C1_whole_epitope_ranked.top20.csv` |
| **C2** | known-Ab, single island | best per (mAb, island); 87 island contigs | ranked, **top-20** per island | **1,660** | `results/dp4_C2_single_island_ranked.top20.csv` |
| **C3** | polyclonal 12-mer tiles | best per window, no antibody | ranked, **top-3** per window † | **1,317** | `results/dp4_C3_12mer_ranked.top20.csv` |
| **C4** | linear 30-mer controls | bare tiled peptides (no scaffold) | exhaustive tiling (fixed) | **2,034** | `data/libraries/dp4_tiled30mers_fasta.csv` |
| **C5** | metric-space titration | designs spread across metrics (calibration) | farthest-point sample (fixed) | **3,000** § | `results/dp4_C5_titration.csv` |
| **C6** | scaffolded-epitope controls | island→Ala + scaffold-disruption | all C1 top-20 base × flavors | **3,100** ‡§ | `results/dp4_C6_controls.csv` |
| **8VDL** | PfEMP1 conserved epitope | two epitope definitions, scaffolded | ranked, **top-10** each | **20** | `results/dp4_8vdl_top10.csv` |
| | | | | **12,251 total** | `data/libraries/dp4_library.csv` |

† **C3 is top-3, not top-20** — the tiling steps by 6 residues so neighbouring windows overlap heavily
(adjacent tiles cover nearly the same epitope space), so best-3 per window suffices; deeper would just
inflate the count with near-duplicates. C1/C2 get the full top-20 the ranked files hold.
‡ **Derived** from C1's top-20 base (island1→A + island2→A dual-only + scaffold disruption), so C6
tracks C1's depth; built at top-20, so depth-20 needs no C6 rebuild.
§ **Native 103 (redo landed 2026-07-11).** C1 originally reused Lawson's 104-residue contigs; we
regenerated it natively at 103 (see *C1 redo at 103*) — RFD3→MPNN→AF3 done, metrics extracted (140,716
designs, all `status==ok`), C1 re-selected, C5/C6 rebuilt off the new pool. **All components are native
103-mers, so the 104→103 assembler trim is a no-op** (kept only as a guard).

Full composite rankings (`results/dp4_*_ranked.csv`) are regenerable and gitignored; only the top-*n*
cuts + case-encoded sequences + the final library are tracked.

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

## Exclusions — the canonical 56-mAb set

Apply **one** exclusion set across the known-Ab components so counts are consistent. From the 59 DP3
known-Ab epitopes, drop three → **56**:
- `4xwo_5P` — low assay yield
- `7a3t_0P` — 4-residue epitope (smallest in DP3, too small to present)
- `2h32_0P` — **not a typical antibody:antigen structure — it is the pre-B cell receptor (pre-BCR)**,
  which binds quite differently; not a valid mAb test case (recorded 2026-02-22). It passed the inclusion
  criteria by accident. All other reviewed DP3 epitopes are well-behaved.

Applies to the **known-Ab components: C1, C2, C5, C6**. **C3** (polyclonal 1D2K/6M0J/4WAT) is
unaffected (different antigens). **C4** (linear tiled controls) **also drops the three** (decided 2026-07-07) so the linear controls stay
consistent with the scaffold arms → 56 mAb + 3 polyclonal = 59 antigens, **2,034** tiles.

Tooling supports it: `stage06_select.py --drop-ids 2h32,4xwo,7a3t`, `stage06_sample_c5.py` `DROP_IDS`
(updated to include 2h32), `build_c6_mutants.py --drop-targets 2h32,4xwo,7a3t`. C4, C5, and C6 are rebuilt on the 56-set; the C1/C2 **rankings** still include all epitopes (a ranking, not a
cut) and C1's `scaffoldEPITOPE` still includes all epitopes (C5 is now on the 56-set), so the 56-set is applied to those at
the final assembly cut (together with the chosen depth).

## Case-encoded `designedSequence` for visualization

Outputs carry the design sequence with **epitope UPPERCASE, scaffold lowercase** (for visualization) — which
is exactly the `scaffoldEPITOPE` column we already produce. Status by component:
- **C1** done (`dp4_C1_scaffoldEPITOPE.csv`, dp2-token route, 104-mer → trim at assembly).
- **C3** done (`dp4_C3_scaffoldEPITOPE.csv`, LOCAL — `scripts/case_encode_c3.py`, 12-mer in `design_seq`, 103-mer).
- **C2** done (`dp4_C2_scaffoldEPITOPE.csv`, `scripts/case_encode_c2.py`, af3 seq + local `n_flank`/`island_size`, 103-mer).
- **C5** done (`dp4_C5_scaffoldEPITOPE.csv`, dp2-token route, 56-set, 104-mer → trim).
- **C4** carries its 30-mer as `designedSequence` already.

In the assembled library file this casing is carried as the
`designedSequence` column, so the visualization casing is carried through to the assembled output.

## Component notes

- **C5 — metric-space titration** (`scripts/stage06_sample_c5.py`). Farthest-point (max–min) spread
  over the standardized scoring axes, ~54 per mAb over the 56 mAbs, 3,000 designs total — deliberately
  including the low-quality tail the filters reject, so binding read off the spread calibrates the scorer.
  **Accessibility is split evenly:** half the sample spans the **real AF3 clash** (`af3_n_clash_res`,
  available because the antibody is known), half the **native-aware cylinder** surrogate, so the titration
  calibrates both the ground-truth accessibility we have here *and* the surrogate we must rely on where no
  antibody is known (the two halves are disjoint and sum to 3,000). Axis coverage (sample range / pool
  range): epitope RMSD 65% (its pool max is a few unfolded outliers the sample doesn't chase), epitope PAE
  100%, overall RMSD 98%, AF3 clash 95%, cylinder 93%. Coverage plot:
  `manuscript/figures/c5_titration_coverage.png` (`episcaf_analysis/viz/plot_c5_titration.py`).
- **C6 — scaffolded-epitope controls** (`episcaf_pipeline/scaffolded_epitope_controls/`). Base = C1
  top-20 over **56 mAbs** (dropped `2h32` pre-BCR, `4xwo` low-yield, `7a3t` 4-residue epitope). Not new
  scaffolding — string substitution on the case-encoded sequence (port of the DP3 mutation-control R code).
  Flavors: every-residue island1→Ala, island2→Ala (dual-island only), and scaffold disruption (`PPDDGG`
  hexamers in scaffold windows, each ≥4 residues from the epitope, seeded). Scaffold disruption is **X4
  with a graceful fallback 4→3→2→1** (`--scaffold-min 1`) rather than dropping the control when 4 don't
  fit. **Shipped build (native-103 C1 pool):** 1,980 island-alanine mutants + 1,120 scaffold-disruption
  controls = **3,100**. Scaffold-disruption fallback distribution: X4 1,034 + X3 60 + X2 25 + X1 1 =
  **1,120/1,120 bases covered** (86 fell back below 4 windows, none dropped). Alanine arms cover all bases.
- **C4 — linear tiled-30mer controls** (`data/libraries/dp4_tiled30mers_fasta.csv`).
  The **full antigen sequences** of 59 antigens (56 mAb targets + 3 tiled antigens 1D2K/6M0J/4WAT; the three excluded mAbs 2h32/4xwo/7a3t are dropped here too for consistency with the 56-mAb set),
  taken from the **FASTA files (no unresolved-gap holes, unlike the PDBs)**, chopped into overlapping
  **30-mers, step = 6**. No RFD/MPNN/AF3. Each 30-mer goes at the END of a constant construct:
  `GSGAGSGA…GSGA` filler + `ENLYFQGA` (TEV protease site) + `[30-mer]` → a constant **103-mer** that is
  cleaved to the 30-mer in the final step (matches the linear-epitope assays). Use the `_fasta` files,
  not the earlier PDB-derived ones.
- **Case-encoding** (`scripts/case_encode_selected.py`, `docs/CASE_ENCODING.md`). C1/C5 designs' epitope
  positions were recovered (token → `dp2.assay_scaffolded_epitope_id` → `scaffolded_epitope_chunk_resindices`)
  and written as case-encoded sequences (`results/dp4_C{1,5}_scaffoldEPITOPE.csv`), feeding C6 + assembly.

## 8VDL — PfEMP1 conserved-epitope arm (`dp4_8vdl/`)

A self-contained seventh arm, John's request, scaffolding a **conserved** epitope from a different
antigen: the EPCR-binding surface of the *Plasmodium falciparum* PfEMP1 CIDRα1.4 domain (crystal
`8VDL`; Reyes et al., *Nature* 636:182–189, 2024). PfEMP1 escapes immunity by antigenic variation, but
the residues its CIDR domain uses to grip the host receptor EPCR cannot vary freely — so presenting that
conserved contact surface is a candidate for eliciting **broadly** reactive, variant-transcending
antibodies. The crystal contains the cognate C7 antibody (chains H/L), so this is a known-antibody target
and the **real clash** term applies (as for C1/C2).

**Two epitope definitions, top-10 each (20 total):**
- **`epitope`** — the whole contiguous contact window **C652–C673** (22 residues; spans all 13 residues
  with a heavy atom within 4 Å of the Fab). The strong constraint, the direct analog of C1.
- **`hotspots`** — only the three functional residues **F655 / F656 / E666** fixed at their native crystal
  coordinates, design builds everything else around them (a minimal "hotspot graft").

**Result — a clean, testable contrast** (from `07_consolidate.py`, which aligns each predicted epitope
onto the native chain-C frame and scores the real H/L clash there): the whole-epitope designs are
antibody-**accessible** (top-10 epitope RMSD 0.97–1.39 Å, only 1–4 clashing residues), while the minimal
hotspot grafts recover the three residues almost exactly (epitope RMSD 0.01–0.22 Å) but **bury** them
under scaffold that would block the antibody (22–52 clashing residues). Comparing the two in the assay
asks whether the minimal hotspot cluster suffices or the whole epitope is needed.

Pipeline: `dp4_8vdl/scripts/` (`01_generate_contigs` → `02_emit_rfd3_inputs` → `03_rfd3_array.sbatch` →
`04_make_fixed_pdbs` → reuses episcaf MPNN/AF3 via the `*_fixed_dldesign_*` naming contract →
`07_consolidate.py`). Its 20 designs are merged into the library at assembly (`category=scaffolded8VDL`).

## Assembly format (the 8-column annotated format)

The final synthesis file is the **8-column annotated format** (the column schema of the earlier PepSeq
library annotation; reference file `episcaf_pipeline/scaffolded_epitope_controls/reference_dp3/DP3_annot.csv`)
— already settled and validated on C4
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

**104→103 truncation — the whole-epitope (C1) family only.** *(INTERIM — superseded once the C1-at-103
rerun lands; see* ***C1 redo at 103*** *below. Kept as the record of what the current 104 pool needs and
why native 103 is cleaner. The assembler still applies this until the new pool is scored.)*
Not all components are 104-mers. **C1
*reproduced* Lawson's whole-epitope run, reusing his contigs** (`contig_length "104-104"`), so **C1 —
and C5 (sampled from C1's pool) and C6 (built from C1) — are 104-residue** proteins (`af3_window_end=104`)
that must be trimmed to the 103-mer assay ceiling. **C2 is natively 103** (we *generated* new contigs
at 103, `build_dual_island_designs.py`, correcting Lawson's 104→103) — no trim. **C3 (12-mer): natively
103** (confirmed — `design_seq` is 103 for all 8,780) — no trim.
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
`0ab98e18e5c6a6de6dc3f9a25881ee10`, mpnn_id 2) — a **two-island** epitope (24 res = a 15-res island at
positions 0–14 and a 9-res island at 95–103) whose islands sit **flush against both termini** (n_flank =
c_flank = 0), so there's no scaffold residue to trim at either end and any trim to 103 clips one epitope
residue. Decision: **drop it and ship the rank-21 design for `3ux9_1P`**
(keeps 20/epitope; the full ranking is regenerable, so rank-21 is well-defined). Apply this at assembly.

## C1 redo at 103 (native 103, supersedes the trim) — DONE 2026-07-11

**Why.** John observed that in DP3, DP2 104-mers truncated to 103 gave generally weaker binding
signal — likely assay run-to-run variation, but enough to make a single-residue truncation a
liability we'd rather not carry. Since C1 is 104 and **C5/C6 derive from C1's pool, 3 of the 4
antibody components ride on it**. So instead of trimming post-hoc we regenerate C1 natively at 103.

**How (`episcaf_pipeline/build_whole_epitope_designs.py`).** Take Lawson's *exact* whole-epitope
contigs from `dp2` (each `N-N/A…/spacer/A…/C-C`, always summing to 104) and drop **one scaffold
residue** — the larger terminal flank, or the largest interior spacer when the islands are flush at
both termini. Every epitope residue and Lawson's inter-island-spacing sweep are preserved; only the
length changes, so C1 stays a faithful DP3 comparator (this is why we *edit* his contigs rather than
resample fresh ones like C2). **Native 103 dissolves the `3ux9_1P` "can't-trim" case** — its spacer
is shortened, both islands kept — so there is **no drop and no rank-21 substitution** anymore.

**Ledger:** `results/whole_epitope_designs.csv` — **2,206 contigs**, 56 mAbs (drops 2h32/4xwo/7a3t),
**all 103**, → 2,206 × 8 RFD3 × 8 MPNN = **141,184 designs = 141,184 AF3 structures** (~1.27× the
single-island run). Verified: 0 island edits, exactly one scaffold residue dropped per contig; `init`
+ `stage01` compile cleanly to `contig_rfd3` at `103-103`.

**Run (Gemini):** `bash scripts/run_whole_epitope_rfd3.sh` (init→stage01→stage02, prints the chunked
RFD3 `sbatch`), then after RFD3 finishes `bash scripts/run_whole_epitope_mpnn_af3.sh runs/whole_epitope_rfd3`
(MPNN wave → AF3 wave). Then re-run stage05 metrics + `stage06_select` for the new C1, re-case-encode
(token→dp2), and rebuild C5/C6 off the new pool. (**8VDL** is run separately as its own arm, `dp4_8vdl/`
— see the *8VDL* section above.)

## Reproduce (exact commands)

Every deliverable is regenerable by a named script; the numbers reported (counts, coverage, the 93 X4
skips) are printed by these scripts, not hand-computed. Metric CSVs live in the local sibling data dirs
(`$D = /Users/bneff/Desktop/projects/episcaf`, see `filesystem-map`); C2 + case-encoding run on Gemini.

```bash
# C1 redo at 103 (local: build the ledger; Gemini: run RFD3->MPNN->AF3)
python scripts/build_whole_epitope_designs.py --drop-targets 2h32,4xwo,7a3t \
  --out results/whole_epitope_designs.parquet      # -> 2,206 contigs, all 103
#   (module form: python -m episcaf_pipeline.build_whole_epitope_designs ...)
bash scripts/run_whole_epitope_rfd3.sh             # Gemini: init->stage01->stage02 + RFD3 sbatch
bash scripts/run_whole_epitope_mpnn_af3.sh runs/whole_epitope_rfd3   # Gemini: after RFD3 done

# C1 selection (local; on the NEW 103 metrics once scored)
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
  --drop-targets 2h32,4xwo,7a3t --out results/dp4_C6_controls.csv

# 8VDL arm (Gemini: RFD3->MPNN->AF3; then consolidate top-10 per definition)
python dp4_8vdl/scripts/07_consolidate.py --out results/dp4_8vdl_top10.csv

# Assemble the final library (local) -> data/libraries/dp4_library.csv, 12,251 constructs
python scripts/stage06_assemble.py --depth 20   # C1/C2 top-20; C3 top-3 (--c3-depth default)
```

Scorer weights/transforms are config, not magic numbers: `episcaf_analysis/presets.py` (provenance
above). C5 and C6 are deterministic (FPS is seed-free deterministic; C6 seeds its RNG).

## Budget & depth (decided)

DP4 = a 36k library that includes all minibinders → **~10–15k slots for Episcaf designs**
(`memory: dp4-budget`). **Depth is settled at top-20 for C1/C2** — the maximum the ranked files hold
(they are top-20 cuts), and the depth C6 was built at, so no C6 rebuild is needed. With **C3 held at
top-3** (neighbouring-tile overlap, see above), the whole library lands at **12,251 constructs**, right
in the 10–15k budget. Per-component at the shipped depth: C1 1,120 · C2 1,660 · C3 1,317 · C4 2,034 ·
C5 3,000 · C6 3,100 · 8VDL 20.

## Pending

0. **C1 redo at 103 — DONE (2026-07-11).** RFD3→MPNN→AF3 on the 2,206-contig 103 ledger → metrics
   (140,716 designs, all `status==ok`) → re-selected C1 top-20 (1,120), re-case-encoded
   (`case_encode_whole_epitope.py`), and rebuilt C5 (3,000) + C6 (3,100) off the new pool. **C1/C5/C6 are
   native 103; the 104→103 trim is now a no-op.**
1. **Assembly (`06_library`) — DONE (2026-07-13).** `scripts/stage06_assemble.py --depth 20` concatenated
   all seven components into `data/libraries/dp4_library.csv` (**12,251** constructs), applying the
   56-exclusion (C1/C2), the top-20/top-3 depth cuts, and global `library_member` numbering. Verified:
   `library_member` and `design_ID` both unique, all sequences 103-mers. (C4's `design_ID`, a per-antigen
   tile index, is namespaced `C4_<target>_t<pos>` so every `design_ID` is globally unique/traceable.)
2. **Oligo encoding — NEXT.** LadnerLab `oligo_encoding` + DP3 codon weights
   (`episcaf_pipeline/oligo_encoding/`), then the order-file step (confirm with Erin).
