# DP4 library — components, selection, and status

Reference for the DP4 PepSeq library: what each component is, how designs were selected, and where
each file lives. Manuscript counterpart: `manuscript/sections/dp4_library.tex` (`sec:dp4`). Related:
`docs/CASE_ENCODING.md`, `docs/CYLINDER_PARAMS.md`.

Status (2026-07-14): assembled. The seven components (C1–C6 plus the 8VDL arm) are selected, built,
case-encoded, and concatenated into `data/libraries/dp4_library.csv` — 15,324 constructs, each a
103-mer with a unique `library_member` and `design_ID`. The shipping depth is top-20 per group for
C1/C2 and top-10 per window for C3; C4/C5/C6 and 8VDL are fixed size. Oligo encoding is the only
remaining step (see Pending).

## Output files

The deliverable and the input to the next step:

- **`data/libraries/dp4_library.csv`** — the library. All seven components merged into one file, 15,324
  rows, in the 8-column PepSeq annotated format **plus 5 scoring columns** (schema below). This is the
  file to hand off. `designedSequence` is the full 103-mer in EPITOPEscaffold casing (epitope uppercase /
  scaffold lowercase) for every row; `sequence` is the plain uppercase 103-mer that gets synthesized.
- **`data/libraries/dp4_named_peptides.csv`** — the oligo-encoder input. A two-column `name,seq` slice of
  the library (no header), for the DNA encoding step. Not a separate result, just a reformat.
  **Regenerate after any library change** (`scripts/stage07_named_peptides.py`).

Behind the library, each scaffolded component has two intermediate files in `results/` (the provenance
the library is assembled from; the depth cut and the 56-mAb exclusion are applied at assembly, so these
hold the full top-20 pools and are larger than the shipped counts):

- **`*_ranked.top20.csv`** — the selection output: which designs were chosen and their scores (plain
  sequence). For C5 this is `dp4_C5_titration.csv`, for C6 `dp4_C6_controls.csv`, for 8VDL `dp4_8vdl_top10.csv`.
- **`*_scaffoldEPITOPE.csv`** — the case-encoded sequences for that component (from the `case_encode_*`
  scripts). C6 encodes by substitution on C1's. C4 and 8VDL have no such file — C4 is case-encoded at
  assembly (tile upper / filler lower) and 8VDL by `07_consolidate` from its contig positions.

## Components

| Comp | Component | What it is | Selection | Constructs | File |
|---|---|---|---|---|---|
| C1 | known-Ab, whole epitope | best scaffolds per mAb (comparator) | ranked, top-20 per mAb | 1,120 § | `results/dp4_C1_whole_epitope_ranked.top20.csv` |
| C2 | known-Ab, single island | best per (mAb, island); 87 island contigs | ranked, top-20 per island | 1,660 | `results/dp4_C2_single_island_ranked.top20.csv` |
| C3 | polyclonal 12-mer tiles | best per window, no antibody | ranked, top-10 per window † | 4,390 | `results/dp4_C3_12mer_ranked.top20.csv` |
| C4 | linear 30-mer controls | bare tiled peptides, no scaffold | exhaustive tiling | 2,034 | `data/libraries/dp4_tiled30mers_fasta.csv` |
| C5 | metric-space titration | designs spread across metrics | farthest-point sample | 3,000 § | `results/dp4_C5_titration.csv` |
| C6 | scaffolded-epitope controls | island→Ala + scaffold disruption | C1 top-20 base × flavors | 3,100 ‡§ | `results/dp4_C6_controls.csv` |
| 8VDL | PfEMP1 conserved epitope | two epitope definitions | ranked, top-10 each | 20 | `results/dp4_8vdl_top10.csv` |
| | | | | 15,324 total | `data/libraries/dp4_library.csv` |

† C3 is shipped deep (top-10). Its windows are **12-mers stepping by 2 residues**, so neighbouring tiles
are highly redundant (adjacent windows share 10 of 12 residues) — which argues for a shallow cut. We
ship top-10 anyway because that redundancy does not protect against the failure mode here: C3 has the
weakest clash distribution of any component, so the per-design success probability is low and
overlapping tiles can fail together. It is also the arm with the most to gain (a few hits in the
no-antibody setting would be the first evidence the approach works without a known antibody), so the
spare budget capacity is spent here (John, 2026-07-14). *(Not to be confused with C4, which tiles
**30-mers at step 6**.)*

‡ C6 is derived from the C1 top-20 base (island1→Ala + island2→Ala for dual-island epitopes + scaffold
disruption), so it tracks C1's depth. It was built at top-20, so depth-20 needs no C6 rebuild.

§ Native 103 (redo landed 2026-07-11). C1 originally reused Lawson's 104-residue contigs; I regenerated
it natively at 103 (see "C1 redo at 103"). RFD3→MPNN→AF3 is done, metrics extracted (140,716 designs,
all `status==ok`), C1 re-selected, and C5/C6 rebuilt on the new pool. Every component is now native 103,
so the 104→103 assembler trim is a no-op (kept only as a guard).

The full composite rankings (`results/dp4_*_ranked.csv`) are regenerable and gitignored; only the
top-*n* cuts, the case-encoded sequences, and the final library are tracked.

## How "best 20" is defined (C1 / C2 / C3)

Designs are ranked, not gated: there is no hard pass/fail cutoff. Each design gets one composite score,
and I take the top *n* per group. Tooling: `scripts/stage06_select.py`; weights and transforms in
`episcaf_analysis/presets.py`; scorer in `episcaf_analysis/score.py`.

Each of four metrics is converted to a percentile within its population and oriented so higher is better
(all four are lower-is-better), then weighted-summed (weights sum to 1):

```
composite = 0.35 · accessibility
          + 0.35 · epitope_RMSD
          + 0.15 · overall_RMSD
          + 0.15 · epitope_PAE
```

- Accessibility is the real AF3 clash (`af3_n_clash_res`) for C1/C2, where the antibody is known, and
  the native-aware cylinder surrogate for C3/C5, where it is not.
- Percentiles are pooled across the mAb set for C1/C2 and taken per-antigen for the 12-mer set (C3).
- Groups for "top *n*": per mAb / `id` (C1); per `(id, island_index)` (C2); per `(antigen, id)` (C3).

The weights come from the DP3 binding data: accessibility and epitope RMSD were the strongest
within-antibody predictors of enrichment (~0.35 each), while overall RMSD and PAE carried little signal
(0.15 each). This is a hand-set prior from a set where every design already passed the filters. C5 is
built to span the metric space so these weights can be re-fit on the real DP4 binding data (manuscript Q2).

## Exclusions — the 56-mAb set

One exclusion set is applied across the known-Ab components so the counts stay consistent. From the 59
DP3 known-Ab epitopes I drop three, leaving 56:
- `4xwo_5P` — low assay yield.
- `7a3t_0P` — 4-residue epitope, the smallest in DP3 and too small to present.
- `2h32_0P` — not an antibody:antigen structure. It is the pre-B cell receptor (pre-BCR), which binds
  differently, and it passed the DP3 inclusion criteria by accident (recorded 2026-02-22).

This applies to C1, C2, C5, and C6. C3 (polyclonal 1D2K/6M0J/4WAT) is unaffected, since those are
different antigens. C4 also drops the three (decided 2026-07-07) so the linear controls match the
scaffold arms: 56 mAb + 3 polyclonal = 59 antigens, 2,034 tiles.

The tooling carries it: `stage06_select.py --drop-ids 2h32,4xwo,7a3t`, `stage06_sample_c5.py` `DROP_IDS`,
`build_c6_mutants.py --drop-targets 2h32,4xwo,7a3t`. C4, C5, and C6 are built on the 56-set. The C1/C2
rankings still include all epitopes (they are rankings, not cuts), so the 56-set is applied to those at
the final assembly cut, along with the depth.

## Case-encoded `designedSequence`

Each design carries its sequence with the epitope in uppercase and the scaffold in lowercase, for
visualization. This is the `scaffoldEPITOPE` column the pipeline already produces, and it is carried
into the assembled library as the `designedSequence` column. By component:
- C1 — `dp4_C1_scaffoldEPITOPE.csv` (native-103 route, `case_encode_whole_epitope.py`, from contig positions).
- C2 — `dp4_C2_scaffoldEPITOPE.csv` (`case_encode_c2.py`, AF3 sequence + local island positions, 103-mer).
- C3 — `dp4_C3_scaffoldEPITOPE.csv` (`case_encode_c3.py`, 12-mer in `design_seq`, 103-mer).
- C5 — `dp4_C5_scaffoldEPITOPE.csv` (native-103 pool, 56-set).
- C4 carries its 30-mer as `designedSequence` directly.

## Component notes

- C5 — metric-space titration (`scripts/stage06_sample_c5.py`). A farthest-point (max–min) spread over
  the standardized scoring axes, about 54 designs per mAb across the 56 mAbs, 3,000 total. It includes
  the low-quality tail the filters reject, so binding read off the spread can calibrate the scorer.
  Accessibility is split evenly: half the sample spans the real AF3 clash (`af3_n_clash_res`, available
  because the antibody is known) and half the native-aware cylinder surrogate, so the titration calibrates
  both the ground-truth accessibility we have here and the surrogate we rely on where no antibody is
  known. The two halves are disjoint and sum to 3,000. Axis coverage (sample range / pool range):
  epitope RMSD 65% (its pool maximum is a few unfolded outliers the sample does not chase), epitope PAE
  100%, overall RMSD 98%, AF3 clash 95%, cylinder 93%. Coverage plot:
  `manuscript/figures/c5_titration_coverage.png` (`episcaf_analysis/viz/plot_c5_titration.py`).
- C6 — scaffolded-epitope controls (`episcaf_pipeline/scaffolded_epitope_controls/`). Base is the C1
  top-20 over the 56 mAbs. These are not new scaffolds; they are string substitutions on the case-encoded
  sequence (a port of the DP3 mutation-control R code). Flavors: every-residue island1→Ala, island2→Ala
  (dual-island only), and scaffold disruption — `PPDDGG` hexamers inserted into scaffold windows, each at
  least 4 residues from the epitope, seeded for reproducibility. Scaffold disruption is X4 with a fallback
  4→3→2→1 (`--scaffold-min 1`) rather than dropping the control when four windows do not fit. Shipped
  build (native-103 C1 pool): 1,980 island-alanine mutants + 1,120 scaffold-disruption controls = 3,100.
  Fallback distribution: X4 1,034 + X3 60 + X2 25 + X1 1 = 1,120 of 1,120 bases covered (86 fell back
  below four windows, none dropped). The alanine arms cover all bases.
- C4 — linear tiled-30mer controls (`data/libraries/dp4_tiled30mers_fasta.csv`). The full antigen
  sequences of 59 antigens (56 mAb targets + 3 tiled antigens 1D2K/6M0J/4WAT; the three excluded mAbs are
  dropped here too), taken from the FASTA files rather than the PDBs to avoid unresolved-gap holes, and
  chopped into overlapping 30-mers at step 6. No RFD/MPNN/AF3. Each 30-mer goes at the end of a constant
  construct: `GSGAGSGA…GSGA` filler + `ENLYFQGA` (TEV site) + `[30-mer]`, giving a 103-mer that is cleaved
  to the 30-mer in the assay. Use the `_fasta` files, not the earlier PDB-derived ones.

## 8VDL — PfEMP1 conserved-epitope arm (`dp4_8vdl/`)

A separate arm (John's request) that scaffolds a conserved epitope from a different antigen: the
EPCR-binding surface of the *Plasmodium falciparum* PfEMP1 CIDRα1.4 domain (crystal `8VDL`; Reyes et
al., *Nature* 636:182–189, 2024). PfEMP1 escapes immunity by antigenic variation, but the residues its
CIDR domain uses to bind the host receptor EPCR are constrained, so a scaffold that presents that
conserved surface is a route to broadly reactive antibodies. The crystal contains the cognate C7
antibody (chains H/L), so this is a known-antibody target and the real clash term applies, as for C1/C2.

Two epitope definitions, top-10 each (20 total):
- `epitope` — the contiguous contact window C652–C673 (22 residues, covering all 13 residues with a
  heavy atom within 4 Å of the Fab). This is the strong constraint, the analog of C1.
- `hotspots` — only F655, F656, and E666, fixed at their native crystal coordinates, with the design
  building around them (a minimal hotspot graft).

Scoring is by `07_consolidate.py`, which aligns each predicted epitope onto the native chain-C frame and
counts the real H/L clash there. The two definitions separate cleanly: the whole-epitope designs stay
accessible (top-10 epitope RMSD 0.97–1.39 Å, 1–4 clashing residues), while the hotspot grafts recover
the three residues almost exactly (epitope RMSD 0.01–0.22 Å) but bury them under scaffold that would
block the antibody (22–52 clashing residues). The assay can then test whether the minimal hotspot cluster
is enough or the whole epitope is needed.

Pipeline: `dp4_8vdl/scripts/` (`01_generate_contigs` → `02_emit_rfd3_inputs` → `03_rfd3_array.sbatch` →
`04_make_fixed_pdbs`, then the shared episcaf MPNN/AF3 via the `*_fixed_dldesign_*` naming contract, then
`07_consolidate.py`). Its 20 designs are merged into the library at assembly (`category=scaffolded8VDL`).

## Assembly format (8-column annotated format)

The synthesis file uses the 8-column annotated format (the column schema of the earlier PepSeq library
annotation; reference `episcaf_pipeline/scaffolded_epitope_controls/reference_dp3/DP3_annot.csv`),
validated on C4 (`dp4_tiled30mers_fasta.csv` is the reference instance):

| column | meaning |
|---|---|
| `library_member` | global id, `DP4_<N>` |
| `sequence` | the 103-mer synthesized (plain uppercase) |
| `category` | component type (`tiled30mer`, `scaffoldedAbEpitope`, `metricSpaceTitration`, …) |
| `model` | `RFD` for designs, `(none)` for linear controls |
| `designedSequence` | the full 103-mer in **EPITOPEscaffold** casing — epitope UPPER, scaffold lower |
| `designedSequenceLength` | len(`designedSequence`) — 103 for every row |
| `design_ID` | per-design id (globally unique) |
| `target` | antigen / mAb id |
| `epitope_rmsd`, `overall_rmsd`, `epitope_pae`, `af3_clashes`, `cylinder_clashes` | the **5 scoring columns** — the metrics each design was selected on; left **blank** where a design has no such value, never imputed |

`sequence` construction: for the linear controls (C4) it is filler + `ENLYFQGA` + 30-mer; for the
scaffolded designs it is the design's own 103-mer. Assembly concatenates all components with global
`library_member` numbering (`scripts/stage06_assemble.py`).

**EPITOPEscaffold casing applies to every row**, including the linear controls: C4's `designedSequence`
is the full 103-mer with the 30-mer tile uppercase (it is the epitope) and the `GSGA…` filler +
`ENLYFQGA` TEV site lowercase (the scaffold). 8VDL is case-encoded from its fixed contig positions.

**Which scoring columns are populated:**

| component | epitope_rmsd | overall_rmsd | epitope_pae | af3_clashes | cylinder_clashes |
|---|---|---|---|---|---|
| C1, C2, C5 | ✓ | ✓ | ✓ | ✓ | ✓ |
| C3 | ✓ | ✓ | ✓ | — (no antibody) | ✓ |
| 8VDL | ✓ | ✓ | ✓ | ✓ | — (cylinder not computed) |
| C4, C6 | — | — | — | — | — (never folded) |

C4 (linear tiles) and C6 (mutants) never went through AF3, so they have none of the five; C3 has no
antibody, so its real clash is undefined; 8VDL is a known-antibody target scored on the real clash and
never needed the cylinder surrogate. Blank cells are honest gaps, not missing work.

Two format notes carried through assembly: (1) the 104→103 trim is now a no-op because every component is
native 103 (background below); (2) C4's `design_ID` was a per-antigen tile-start index that repeated
across antigens, so it is namespaced `C4_<target>_t<pos>` to make every `design_ID` in the library
unique and traceable.

### Background: the 104→103 trim (now a no-op)

This is kept as a record of what the old 104-mer C1 pool needed and why native 103 is cleaner; the
assembler still carries the trim as a guard. C1 originally reproduced Lawson's whole-epitope run and
reused his contigs (`contig_length "104-104"`), so C1 — and C5 (sampled from C1) and C6 (built from C1)
— were 104-residue proteins that had to be trimmed to the 103-mer assay ceiling. C2 and C3 were already
native 103. The trim dropped one residue from whichever terminus is scaffold (C-terminal by default,
N-terminal when the C-terminus is epitope), so no epitope residue was lost. One case could not be trimmed
either way: `3ux9_1P` rank 9, a two-island epitope whose islands sit flush against both termini, where
any trim to 103 clips an epitope residue; the plan was to drop it and ship the rank-21 design. Native 103
(below) removes this case entirely.

The trim is reasoned from folding principles, not re-validated: a single terminal residue is the least
structurally committed position in a fold, and we cut the scaffold terminus, away from the epitope. The
DP3 binding data was measured on the 104-mers, so a spot re-fold of a 103-mer sample would close that gap
if needed. Full detail in manuscript `sec:production`.

## C1 redo at 103 (native 103) — done 2026-07-11

Why: John noted that in DP3, 104-mers truncated to 103 gave generally weaker binding signal. It is
probably assay run-to-run variation, but enough to avoid carrying a single-residue truncation. Since C1
was 104 and C5/C6 derive from it, three of the four antibody components rode on that truncation. So
instead of trimming after the fact, I regenerated C1 natively at 103.

How (`episcaf_pipeline/build_whole_epitope_designs.py`): take Lawson's exact whole-epitope contigs from
`dp2` (each `N-N/A…/spacer/A…/C-C`, summing to 104) and drop one scaffold residue — the larger terminal
flank, or the largest interior spacer when the islands are flush at both termini. Every epitope residue
and Lawson's inter-island spacing are preserved; only the length changes, so C1 stays a direct DP3
comparator. This is why I edit his contigs rather than resample fresh ones as C2 does. Native 103 also
dissolves the `3ux9_1P` case: its spacer is shortened and both islands are kept, so there is no drop and
no rank-21 substitution.

Ledger: `results/whole_epitope_designs.csv` — 2,206 contigs, 56 mAbs, all 103, giving
2,206 × 8 RFD3 × 8 MPNN = 141,184 designs / AF3 structures. Verified: no island edits, exactly one
scaffold residue dropped per contig.

Run (Gemini): `bash scripts/run_whole_epitope_rfd3.sh` (init→stage01→stage02 + the chunked RFD3 sbatch),
then after RFD3 finishes `bash scripts/run_whole_epitope_mpnn_af3.sh runs/whole_epitope_rfd3` (MPNN wave
→ AF3 wave). Then re-run stage05 metrics + `stage06_select` for the new C1, re-case-encode, and rebuild
C5/C6 on the new pool. (8VDL is run separately; see its section above.)

## Reproduce (exact commands)

Every file is regenerable by a named script, and the reported counts and coverage numbers are printed by
these scripts, not computed by hand. Metric CSVs live in the local sibling data dirs
(`$D = /Users/bneff/Desktop/projects/episcaf`, see `filesystem-map`); C2 and case-encoding run on Gemini.

```bash
# C1 redo at 103 (local: build the ledger; Gemini: run RFD3->MPNN->AF3)
python scripts/build_whole_epitope_designs.py --drop-targets 2h32,4xwo,7a3t \
  --out results/whole_epitope_designs.parquet      # -> 2,206 contigs, all 103
bash scripts/run_whole_epitope_rfd3.sh             # Gemini: init->stage01->stage02 + RFD3 sbatch
bash scripts/run_whole_epitope_mpnn_af3.sh runs/whole_epitope_rfd3   # Gemini: after RFD3 done

# C1 selection (local; on the 103 metrics)
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

# Case-encode C1 + C5 (Gemini) -> results/dp4_C{1,5}_scaffoldEPITOPE.csv
sbatch scripts/case_encode_selected.sbatch

# C6 (local; seeded, reproducible)
python episcaf_pipeline/scaffolded_epitope_controls/build_c6_mutants.py \
  --input results/dp4_C1_scaffoldEPITOPE.csv \
  --id-col token --target-col target --seq-col scaffoldEPITOPE \
  --drop-targets 2h32,4xwo,7a3t --out results/dp4_C6_controls.csv

# 8VDL arm (Gemini: RFD3->MPNN->AF3; then consolidate top-10 per definition)
python dp4_8vdl/scripts/07_consolidate.py --out results/dp4_8vdl_top10.csv

# Assemble the library (local) -> data/libraries/dp4_library.csv, 15,324 constructs
python scripts/stage06_assemble.py --depth 20   # C1/C2 top-20; C3 top-10 (--c3-depth default)

# Export the oligo-encoder input -> data/libraries/dp4_named_peptides.csv
python scripts/stage07_named_peptides.py \
  --library data/libraries/dp4_library.csv --out data/libraries/dp4_named_peptides.csv
```

Scorer weights and transforms are config, not magic numbers (`episcaf_analysis/presets.py`). C5 and C6
are deterministic (FPS is seed-free deterministic; C6 seeds its RNG).

## Budget and depth

DP4 is a 36k library that includes all minibinders, which leaves roughly 10–15k slots for Episcaf designs
(`memory: dp4-budget`). Depth is set at top-20 for C1/C2 — the most the ranked files hold, and the depth
C6 was built at, so no C6 rebuild is needed. **C3 is top-10** (2026-07-14): John flagged that the spare
capacity (~2k slots to reach 36k) is best spent maximizing polyclonal hits, given C3's weak clash
distribution. That brings the library to **15,324**. Per component: C1 1,120, C2 1,660, C3 4,390,
C4 2,034, C5 3,000, C6 3,100, 8VDL 20.

C3 depth is the one elastic dial left: top-3 = 1,317, top-5 = 2,195, top-10 = 4,390 (shipped). Set it
with `stage06_assemble.py --c3-depth <n>` if the final minibinder count moves the headroom.

## Pending

0. C1 redo at 103 — done 2026-07-11. RFD3→MPNN→AF3 on the 2,206-contig ledger, metrics (140,716 designs,
   all `status==ok`), C1 re-selected to top-20 (1,120), re-case-encoded (`case_encode_whole_epitope.py`),
   C5 (3,000) and C6 (3,100) rebuilt on the new pool. C1/C5/C6 are native 103; the 104→103 trim is a no-op.
1. Assembly (`06_library`) — done 2026-07-13. `scripts/stage06_assemble.py --depth 20` concatenated the
   seven components into `data/libraries/dp4_library.csv` (15,324 constructs), applying the 56-exclusion
   (C1/C2), the top-20 (C1/C2) and top-10 (C3) depth cuts, and global numbering. `library_member` and
   `design_ID` are unique and every sequence is 103 residues. Ships the 8 annotation columns + the 5
   scoring columns, with `designedSequence` in EPITOPEscaffold casing on every row.
2. Oligo encoding — in progress. Encoder input is exported and validated:
   `scripts/stage07_named_peptides.py` → `data/libraries/dp4_named_peptides.csv` (`name,seq`, no header,
   all 103-mers; regenerate after any library change). Next: run the LadnerLab encoder on Gemini (`main` step 1 → `encoding_with_nn.py`
   step 2, DP3 recipe + `codon_weights_updated.csv`; `episcaf_pipeline/oligo_encoding/`), then the order-file
   and Twist-adapter step (confirm the exact recipe with Erin).
