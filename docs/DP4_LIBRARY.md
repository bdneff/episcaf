# DP4 library — components, selection, and status

Reference for the DP4 PepSeq library: what each component is, how designs were selected, and where
each file lives. Manuscript counterpart: `manuscript/sections/dp4_library.tex` (`sec:dp4`). Related:
`docs/CASE_ENCODING.md`, `docs/CYLINDER_PARAMS.md`.

Status (2026-07-16): **built, encoded, and synthesis-ready.** The seven components (C1–C6 plus the 8VDL
arm) are selected, built, case-encoded, and concatenated into `data/libraries/dp4_library.csv` — 15,324
constructs, each a 103-mer with a unique `library_member` and `design_ID`. The shipping depth is top-20
per group for C1/C2 and top-10 per window for C3; C4/C5/C6 and 8VDL are fixed size. Selection ran under
the soft-gate scorer (`antibody_softgate`). The library has since been oligo-encoded and gated into
`data/libraries/dp4_order_file.csv` — 37,083 oligos, all verified — which is the file that goes to Twist.

> **Minibinder arm added (2026-07-20).** `dp4_library.csv` now also carries the **21,759** filter-passing
> **LX PfEMP1/EPCR minibinders** (`category=minibinder`), so the file is a single view of the whole DP4
> library — **37,083 rows** (15,324 episcaf + 21,759 minibinder). These are a separate de-novo binder arm
> (same PfEMP1 project as 8VDL), not scaffolded or scored by episcaf, so their five episcaf-metric columns
> are blank — but **every native LatentX column is carried as `lx_<name>`** (plddt, pae, rmsd, ipae,
> iptm, plddt_binder, hotspots, uuid, …) for post-hoc analysis. The library also carries the full metric
> + scoring set (not condensed to the lean 5): the PAE decomposition, ptm, `composite`, `rank_in_group`,
> `is_global_pass`, `island_index`. **33 columns** total (episcaf rows blank in `lx_`). The **15,324** count elsewhere in this doc refers to the episcaf-scaffolded portion — what was
> selected and case-encoded. **The oligo order file covers the WHOLE library** (confirmed 2026-07-20):
> all **37,083** 103-mers are oligo-encoded together into one PepSeq assay, minibinders included (they
> carry no oligos of their own). This is how the lab runs PepSeq — several projects combined into a
> single assay, tracked in one file. `--exclude-category minibinder` would revert to episcaf-only.
> Added by `dp4_8vdl/scripts/08_add_minibinders.py` (idempotent; run after assembly). The current
> committed `dp4_order_file.csv` is stale (15,324, pre-minibinder) — pending re-encode at 37,083.

## Paths used in this doc

Three shell variables appear throughout. Define them before running anything here:

```bash
# On the cluster (Gemini):
export REPO=/scratch/bneff/episcaf                                       # the git checkout (DISPOSABLE
                                                                         # -- it is a clone; re-clone it)
export WS=/tgen_labs/altin/alphafold3/workspace/episcaf_v2_bneff         # the DURABLE workspace

# Locally (laptop): $D is the parent of the repo, holding the non-git sibling data dirs
export D=/Users/bneff/Desktop/projects/episcaf                           # has known_antigen/, 12mer_tiling/
```

The rule that matters: **`/scratch` is ephemeral (swept on ~30 days); `/tgen_labs` is persistent.**
Anything long-lived belongs under `$WS`. Never `rsync --delete` toward `/tgen_labs`, and never
`git init` inside a data directory. Full map: memory `filesystem-map`.

## Output files

The deliverable and the input to the next step:

- **`data/libraries/dp4_library.csv`** — the library. All seven components (+ minibinders) merged into
  one file, 37,083 rows, in the 8-column PepSeq annotated format **plus the full metric + scoring set and
  the `lx_` minibinder columns** — 33 columns, schema below. This is the
  file to hand off. `designedSequence` is the full 103-mer in EPITOPEscaffold casing (epitope uppercase /
  scaffold lowercase) for every row; `sequence` is the plain uppercase 103-mer that gets synthesized.
- **`data/libraries/dp4_order_file.csv`** — the Twist synthesis order file. 37,083 oligos, two columns
  (`Seq ID`, `nucleotide_encoding_with_twist_adapters`). Every row verified by
  `scripts/stage07_order_file.py`: 349 nt, the 20-mer adapters on both ends, and a core that translates
  back to exactly its own peptide. **This is what goes to Twist.**
- **`data/libraries/dp4_superset.csv.gz`** — the all-designs superset (committed gzipped, 357,789 rows
  across all arms; see below). Not a deliverable — it exists for looking at the distributions the library
  was drawn from, and it is a true superset (the shipped library is a strict subset of it).
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
| | | | | **15,324 episcaf** | `data/libraries/dp4_library.csv` (episcaf portion) |
| LX | PfEMP1/EPCR minibinders | de-novo binders (not scaffolded) | LX-filter passers | 21,759 | `dp4_8vdl/scripts/08_add_minibinders.py` (source gitignored) |
| | | | | **37,083 total** | `data/libraries/dp4_library.csv` (whole file) |

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

Designs are ranked, not gated: there is no hard pass/fail cutoff, so no target is ever dropped for want
of a good enough design. Each design gets one composite score and I take the top *n* per group. Tooling:
`scripts/stage06_select.py`; weights and transforms in `episcaf_analysis/presets.py`; scorer in
`episcaf_analysis/score.py`.

**C1/C2 use the `antibody_softgate` preset** (adopted 2026-07-16). Each metric is squashed by its own
sigmoid and the four are weighted-summed:

> **In flux (2026-07-20):** the shipped `dp4_library.csv` was selected under this preset with the
> **epitope-PAE midpoint at 5** (its C1 members reproduce that ranking exactly). The midpoint has since
> been retuned to **2.5** (data-driven — see below), which swaps ~107 of the 1120 C1 selections (~10%),
> counts and components unchanged. Re-selection under 2.5 (`stage06_select` → `stage06_assemble`) is
> pending; nothing is synthesized yet, so the shift is free. Until it runs, the code scorer (2.5) and
> the shipped file (5) differ by design.

```
composite = 0.45 · sigma(af3_n_clash_res;    midpoint 6,   k 0.5)   accessibility -- ranked
          + 0.25 · sigma(epitope_chunk_rmsd; midpoint 1,   k 4.0)   fidelity      -- soft gate
          + 0.20 · sigma(overall_rmsd;       midpoint 2,   k 4.0)   fold          -- soft gate
          + 0.10 · sigma(epitope_pae;        midpoint 2.5, k 1.2)   epitope rigidity -- rank nudge

sigma(x) = 1 / (1 + exp(k * (x - midpoint)))          all four metrics are lower-is-better
```

The point of the steepnesses: broad on clash (k 0.5) so accessibility is genuinely *ranked* across its
range, and steep on the two fold metrics (k 4) so they act as **gates** — a design far the wrong side of
the threshold scores near zero on that term. As k grows the sigmoid approaches a step, and in the limit
this IS Lawson's hard four-filter. Keeping k finite is the whole trick: nothing is ever fully zeroed, so
a target whose designs are all mediocre still contributes its best ones rather than vanishing from the
library. That is what "soft gate" means here.

The `epitope_pae` term is gentler (k 1.2) — a rank *nudge*, not a gate — and its midpoint is **2.5**,
set from the data rather than borrowed from the global `mean_pae < 5` threshold. `epitope_pae` is the
intra-epitope PAE block (short-range pairs), so it runs well below the whole-matrix `mean_pae`: a good
epitope sits near ~2 A (four-filter passers median 1.85; ~1.98 even with the epitope-RMSD filter
removed, so it is not merely a conditioning artifact of the r≈0.8 epitope-PAE/RMSD correlation), while a
marginal one sits ~3.6. So 2.5 is the half-credit point; the old 5.0 sat in the tail and barely
discriminated. This term encodes a **design assumption** — that a scaffold presents its epitope best
when the epitope is *rigid* — which may be wrong; the spanning subset is what tests it experimentally.
Provenance and the measurement live in `episcaf_analysis/presets.py`.

**Global-pass promotion** (`pass_bonus`) sits on top: `composite += 2.0 * P`, where `P` is the product of
four steep sigmoids (k 12) at the four-filter thresholds (`epitope_chunk_rmsd` 1, `overall_rmsd` 2,
`mean_pae` 5, `af3_n_clash_res` 0.5). `P` is ~1 only if a design clears *all four* — a product, not a sum,
so it is a soft AND rather than a count of how many filters passed. A gain of 2 exceeds the composite's
own range, so every four-filter passer floats above every non-passer, and the composite then breaks ties
within each band. This is John's rule ("all global-passing designs ranked above any non-global-passing")
implemented without a hard gate.

- Accessibility is the real AF3 clash (`af3_n_clash_res`) for C1/C2, where the antibody is known, and the
  native-aware cylinder surrogate for C3/C5, where it is not.
- **C3 uses the `twelvemer` preset, not the soft-gate** — there is no antibody, so there is no real clash
  to gate on, and accessibility rests entirely on the cylinder.
- Groups for "top *n*": per mAb / `id` (C1); per `(id, island_index)` (C2); per `(antigen, id)` (C3).

The weights are a hand-set prior from the DP3 binding data, where accessibility and epitope RMSD were the
strongest within-antibody predictors of enrichment and overall RMSD and PAE carried little signal — and
that data is itself a set in which every design had already passed the filters. C5 is built to span the
metric space so these dials can be re-fit on the real DP4 binding data (manuscript Q2).

**Selection is budget-bound, not scorer-bound** (measured on the superset, 2026-07-16). Of C1's 727
four-filter passers only 184 shipped, which sounds like the scorer failing until you group it: the
passers sit in just 15 of the 56 targets, and 7 of those hold more passers than their 20 slots (`7ox3_0P`
alone has 360). `sum over targets of min(passers, 20)` is exactly 184 — the number that shipped — and in
zero targets does a non-passing design outrank a passing one. The promotion does what it says; the
per-target budget is what binds.

## The all-designs superset (`$WS/dp4_superset.csv`)

`dp4_library.csv` holds only the designs that shipped, which makes it impossible to ask the obvious
question: what did they beat? The superset answers that — and it is a **true superset**, so the shipped
library is a strict subset of it. It is every candidate design across **every arm** — **357,789** rows
(C1 140,716 + C2 111,322 + C3 82,712 + 8VDL 1,280 + 21,759 passing LX minibinders) — in the **union of
the library's and the superset's columns** (36): the scoring internals (`selected`, `library_member`,
`is_global_pass`, `composite`, `rank_in_group`, full PAE decomposition) plus the library's synthesis and
minibinder columns (`model`, `designedSequenceLength`, `design_ID`, the 13 `lx_*`). Episcaf/8VDL rows are
blank in `lx_*`; minibinder rows blank in the episcaf metric/scoring columns (never scored on our axes).
So a design that shipped and a design that lost sit in the same table, in the same shape. Of these,
**28,949 shipped** (C1 1,120 + C2 1,660 + C3 4,390 + 8VDL 20 + LX 21,759) — the library's 37,083 minus
the C4/C5/C6 controls, which aren't candidate-pool designs — and the four-filter passers are C1 727 +
C2 407 (C3/minibinders have no antibody-based global pass defined).

It is ranked under the same preset that picked the library (`antibody_softgate` for C1/C2, `twelvemer`
for C3), which makes the ranks reconcile rather than merely resemble: every selected design is exactly
the top-*n* of its group, 56 targets × ranks 1–20 for C1 and 439 windows × ranks 1–10 for C3, with no
gaps. That equality is worth treating as a test — if a future scorer change breaks it, the superset and
the shipped library have drifted apart.

`sequence` is filled for the selected designs (copied verbatim from `dp4_library.csv`, so the two agree
by construction) and for the global-passing ones (read from each design's AF3 chain A), and left blank
for the rest of the episcaf pool. Filling every one would mean reading every design PDB, and the
distributions this file exists for live in the metrics, not the sequences. C3 is the exception and comes
out fully sequenced (82,712), because its metrics already carry `design_seq` — free, so we take it; the
minibinders likewise arrive fully sequenced from the LX file. `designedSequence` is
selected-only, since case-encoding was only ever run on the selections.

Build it in **one cluster pass**: `sbatch scripts/build_superset.sbatch` does everything —
(a) the C1/C2/C3 episcaf pool into an `*_episcaf.csv` intermediate (the sequence pass reads
`runs/*/04_af3/outputs`, and because the C1/C2 metrics record the **absolute `/scratch` paths** those runs
were built at, the job passes `--af3-remap /scratch/bneff/episcaf/runs:$WS/runs` to point them at the
durable copies — without it the AF3 lookups resolve nothing once `/scratch` is swept, and it hard-fails
rather than writing a blank `sequence` column); (b) `extend_superset.py` folds in the 8VDL candidates and
the passing minibinders and unions the columns → the full **357,789-row / 36-column** `$WS/dp4_superset.csv`;
(c) gzips it to `data/libraries/dp4_superset.csv.gz` for you to commit. So `$WS/dp4_superset.csv` and the
committed `.gz` are **the same file** — no more partial-vs-full confusion. The extend step needs the LX
minibinder source on the cluster (`dp4_8vdl/data/LX_*.csv`, ~168 MB, gitignored — `rsync` it up once); the
job fails loudly if it is missing. The raw `.csv` (~116 MB) is gitignored; the gzipped copy (~34 MB) is
committed. Everything is idempotent. `scripts/extend_superset.py` can also be run standalone if the
episcaf intermediate already exists.

Verified on the 2026-07-16 build: every shipped member matched a design in its pool (1120/1120 C1,
1660/1660 C2, 4390/4390 C3), and every passing design's sequence was readable (543/543 C1, 233/233 C2 —
the rest of the passers had already shipped and got theirs from the library).

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

Scoring is by `07_consolidate.py`, which aligns each predicted epitope onto the native chain-C frame,
counts the real H/L clash there, and (2026-07-20) ranks under the **same `antibody_softgate` scorer as
C1/C2** — real Fab clash, no cylinder surrogate. The two definitions separate cleanly, and the hotspot
one fails informatively:
- **epitope** (contiguous C652–673): top-10 epitope RMSD 0.97–1.36 Å, **1–2** clashing residues — the
  soft-gate found accessible designs.
- **hotspots** (F655/F656/E666): top-10 epitope RMSD 0.04–0.11 Å but **39–71** clashing residues. This
  is **not** a selection artifact — across all **640** hotspot designs the clash floor is **14** and
  none fall below 10, so an accessible minimal-hotspot graft was never generated (generation-limited,
  the same phenomenon as C1). The three discontiguous residues appear to force the scaffold into the
  Fab footprint. We ship the top-10 anyway as a **documented negative result** — a limitation of the
  minimal-graft approach, not a scorer failure. (Per-design metrics: `results/dp4_8vdl_top10_allmetrics.csv`.)

Pipeline: `dp4_8vdl/scripts/` (`01_generate_contigs` → `02_emit_rfd3_inputs` → `03_rfd3_array.sbatch` →
`04_make_fixed_pdbs`, then the shared episcaf MPNN/AF3 via the `*_fixed_dldesign_*` naming contract, then
`07_consolidate.py`). Its 20 designs are merged into the library at assembly (`category=scaffolded8VDL`).

## Assembly format (8-column annotated format)

The synthesis file uses the 8-column annotated format (the column schema of the earlier PepSeq library
annotation; reference `episcaf_pipeline/scaffolded_epitope_controls/reference_dp3/DP3_annot.csv`),
validated on C4 (`dp4_tiled30mers_fasta.csv` is the reference instance):

### Column dictionary (33 columns)

Every column, what it means, and **which rows are blank and why** — nothing is ever imputed, so a blank
is always "not measured for this kind of design," never "measured as zero." The fill pattern is measured
from the file itself (`● filled, · blank`; components: C1 whole-epitope, C2 single-island, C3 polyclonal
12-mer, C4 linear tiles, C5 titration, C6 mutant controls, 8V = 8VDL, LX = minibinders):

| column | meaning | C1 C2 C3 C4 C5 C6 8V LX | why blank where it is |
|---|---|---|---|
| `library_member` | global id `DP4_<N>` | ● ● ● ● ● ● ● ● | never blank |
| `sequence` | the 103-mer synthesized (plain uppercase) | ● ● ● ● ● ● ● ● | never blank — this is what's ordered |
| `category` | component type (`scaffoldedAbEpitope`, `minibinder`, …) | ● ● ● ● ● ● ● ● | never blank |
| `model` | design method: `RFD` (RFdiffusion), `LX` (LatentX minibinder), `(none)` (C4 linear) | ● ● ● ● ● ● ● ● | never blank |
| `designedSequence` | the 103-mer in **EPITOPEscaffold** casing (epitope UPPER, scaffold lower); for minibinders = the plain sequence (no epitope to case) | ● ● ● ● ● ● ● ● | never blank |
| `designedSequenceLength` | `len(designedSequence)` — 103 everywhere | ● ● ● ● ● ● ● ● | never blank |
| `design_ID` | per-design id (globally unique) | ● ● ● ● ● ● ● ● | never blank |
| `target` | antigen / mAb id (minibinders: `fold_pfemp1_epcr_*`) | ● ● ● ● ● ● ● ● | never blank |
| `epitope_rmsd` | epitope-backbone RMSD to native | ● ● ● · ● · ● · | **C4/C6 were never folded** (linear tiles / sequence-only mutants); **LX** not scored on our axes |
| `overall_rmsd` | whole-construct backbone RMSD | ● ● ● · ● · ● · | same as above |
| `epitope_pae` | intra-epitope PAE block | ● ● ● · ● · ● · | same as above |
| `scaffold_pae` | scaffold×scaffold PAE block | ● ● ● · ● · ● · | same as above (8VDL: added 2026-07-20) |
| `mean_pae` | whole-matrix PAE (the four-filter's "overall PAE") | ● ● ● · ● · ● · | same as above |
| `ptm` | AF3 predicted TM-score | ● ● ● · ● · ● · | same as above |
| `af3_clashes` | **real** antibody (H/L) clash-residue count | ● ● · · ● · ● · | **C3 has no antibody** (polyclonal) → uses the cylinder instead; C4/C6/LX unfolded/unscored |
| `cylinder_clashes` | **cylinder-surrogate** clash (accessibility proxy) | ● ● ● · ● · · · | **8VDL has the real Fab** so no surrogate was computed; C4/C6/LX unfolded/unscored. C1/C2/C5 carry BOTH clash flavors (the whole-epitope pipeline computes both) |
| `composite` | the `antibody_softgate` (or `twelvemer` for C3) composite score | ● ● ● · · · ● · | **C5 is a farthest-point metric-space sample, not composite-ranked**; C4/C6/LX not ranked |
| `rank_in_group` | rank within the design's selection group (1 = best) | ● ● ● · · · ● · | same as `composite` — only the composite-ranked arms |
| `is_global_pass` | clears all four Lawson filters (soft-AND > 0.5) | ● ● · · · · ● · | **C3 has no antibody**, so no antibody-based global pass is defined; C4/C5/C6/LX not applicable |
| `island_index` | which island (0-based) the design scaffolds | · ● · · · · · · | **only C2 is per-island**; C1 scaffolds the whole epitope (no split), everything else is single-target |
| `lx_*` (13 cols) | every native LatentX column (`lx_plddt`, `lx_pae`, `lx_rmsd`, `lx_ipae`, `lx_iptm`, `lx_plddt_binder`, `lx_hotspots`, `lx_uuid`, …) | · · · · · · · ● | **only the minibinders have a LatentX record**; all episcaf rows are blank here (`lx_iptm` is also blank on many LX rows — LatentX left it empty at the source) |

Two non-obvious blanks worth remembering: **C5 has no `composite`/`rank`** because it is sampled to
*span* the metric space, not to top-rank it (that is the whole point of the titration arm); and **C3 has
no `is_global_pass` or `af3_clashes`** because it is the polyclonal / no-antibody arm, where accessibility
is the cylinder surrogate and there is no single antibody to define a clash or a four-filter pass against.

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
these scripts, not computed by hand. `$D`, `$REPO`, and `$WS` are defined under **Paths used in this doc**
above — export them first. Metric CSVs live in the local sibling data dirs under `$D`; C2 and
case-encoding run on Gemini.

```bash
# C1 redo at 103 (local: build the ledger; Gemini: run RFD3->MPNN->AF3)
python scripts/build_whole_epitope_designs.py --drop-targets 2h32,4xwo,7a3t \
  --out results/whole_epitope_designs.parquet      # -> 2,206 contigs, all 103
bash scripts/run_whole_epitope_rfd3.sh             # Gemini: init->stage01->stage02 + RFD3 sbatch
bash scripts/run_whole_epitope_mpnn_af3.sh runs/whole_epitope_rfd3   # Gemini: after RFD3 done

# C1 selection (on the native-103 metrics; --preset antibody_softgate is the adopted scorer)
python scripts/stage06_select.py --preset antibody_softgate \
  --metrics-csv $D/known_antigen/analysis/data/metrics_whole_epitope_103.csv \
  --group id --topk 20 --out results/dp4_C1_whole_epitope_ranked.csv

# C2 (Gemini, at the dual-island run's metrics; soft-gate is where it matters -- clash 6->2, 6cyf 14.5->3)
python scripts/stage06_select.py --preset antibody_softgate \
  --metrics-csv runs/dual_island_rfd3/05_analysis/metrics_dual_island.parquet \
  --group id,island_index --topk 20 --out results/dp4_C2_single_island_ranked.csv

# C3 (local; polyclonal/no-antibody -> twelvemer preset, NOT soft-gate; cylinder accessibility)
python scripts/stage06_select.py --preset twelvemer \
  --metrics-csv $D/12mer_tiling/analysis/data/metrics_12mer.csv \
  --group antigen,id --topk 20 --out results/dp4_C3_12mer_ranked.csv

# C4 (local; defaults -> data/libraries/dp4_tiled30mers_fasta.csv)
python -m episcaf_pipeline.build_dp4_tiled30mers_fasta

# C5 (local; deterministic FPS -- scorer-independent, samples the metric space of the 103 pool)
python scripts/stage06_sample_c5.py \
  --metrics-csv $D/known_antigen/analysis/data/metrics_whole_epitope_103.csv \
  --total 3000 --out results/dp4_C5_titration.csv

# Case-encode the selections (Gemini). NOTE: the native-103 C1 uses the CONTIG-POSITION encoder
# (case_encode_whole_epitope.py) -- NOT case_encode_selected.sbatch, which is the old token->dp2 path
# and errors with KeyError 'token' on the 103 ranked file.
python scripts/case_encode_whole_epitope.py \
  --selected results/dp4_C1_whole_epitope_ranked.top20.csv \
  --ledger   results/whole_epitope_designs.csv \
  --out      results/dp4_C1_scaffoldEPITOPE.csv       # C1 (native 103)
sbatch scripts/case_encode_c2.sbatch                  # C2 (single-island)
sbatch scripts/case_encode_selected.sbatch            # C5 titration (token->dp2); C3: case_encode_c3.py

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

# Oligo-encode the library (Gemini). Run BOTH steps in the same working dir; step 1 is the long pole
# (C++ sampler, hours). ADAPTER is pinned to the 20-mers inside encode_step2_select.sbatch -- do NOT
# rely on the encoder's own --adapter default, which is the 19-mer form (see oligo-adapter-trap).
cd $REPO/runs/dp4_encoding_full     # with dp4_named_peptides.csv + codon_weights_updated.csv present
# TOOL_DIR is REQUIRED (no default, by design -- a bad default once ran an unexecutable binary). Point it
# at the built LadnerLab encoder; the env var must come BEFORE `sbatch`, not after (else it is read as the
# script path). Same TOOL_DIR for both steps. Step 2 pins the 20-mer ADAPTER internally -- do not pass it.
TOOL=/home/bneff/Library-Design/oligo_encoding
JID1=$(TOOL_DIR=$TOOL INPUT=dp4_named_peptides.csv sbatch --parsable --time=12:00:00 \
  $REPO/episcaf_pipeline/oligo_encoding/encode_step1_generate.sbatch)
TOOL_DIR=$TOOL sbatch --dependency=afterok:$JID1 \
  $REPO/episcaf_pipeline/oligo_encoding/encode_step2_select.sbatch

# Emit + VERIFY the Twist order file (Gemini). Checks every row: 20-mer adapters, 349 nt, and that each
# core translates back to its own peptide. Writes nothing if any row fails.
python $REPO/scripts/stage07_order_file.py \
  --best-encodings $REPO/runs/dp4_encoding_full/DP4_best_encodings \
  --peptides       $REPO/runs/dp4_encoding_full/dp4_named_peptides.csv \
  --out            $REPO/data/libraries/dp4_order_file.csv     # -> 37,083 oligos, all verified

# ALL-DESIGNS SUPERSET (John's ask -- every candidate design, not just the selected ones). ONE cluster
# pass: build C1/C2/C3 -> extend with 8VDL + passing minibinders -> gzip. Needs the LX source on-cluster
# (rsync dp4_8vdl/data/LX_*.csv up once, ~168MB, gitignored). Emits $WS/dp4_superset.csv (== the .gz).
sbatch scripts/build_superset.sbatch
git add data/libraries/dp4_superset.csv.gz && git commit -m "superset: rebuild" && git push
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
   `design_ID` are unique and every sequence is 103 residues. Ships the full **33-column schema** (see the
   column dictionary): 8 annotation columns + the metric/scoring set + the 13 `lx_` minibinder columns,
   with `designedSequence` in EPITOPEscaffold casing on every scaffolded row.
2. Oligo encoding — **whole library encoded 2026-07-20** (37,083 = 15,324 episcaf + 21,759 minibinders;
   supersedes the 2026-07-16 episcaf-only 15,324 encode after the 2.5 re-selection + soft-gate 8VDL +
   minibinder addition). Encoder input: `scripts/stage07_named_peptides.py` →
   `data/libraries/dp4_named_peptides.csv` (`name,seq`, no header, all 103-mers; regenerate after any
   library change). Ran on Gemini in `runs/dp4_encoding_full/` with the LadnerLab encoder (step 1 sampler
   → step 2 NN selector, DP3 recipe + `codon_weights_updated.csv`; `episcaf_pipeline/oligo_encoding/`, see
   its README and manuscript `sec:oligo`) → `DP4_best_encodings`, all 37,083, nothing dropped. The order
   file is **not** a further encoding step — step 2 already emits the adapter-flanked oligo, so the order
   file is a two-column slice, emitted + verified by `scripts/stage07_order_file.py` →
   `data/libraries/dp4_order_file.csv`. All 37,083 verified: 349 nt, 20-mer adapters both ends, every
   core translates back to its own peptide, one encoding per peptide, none missing. **Twist-ready.**
   **Adapter length resolved 2026-07-16: John confirmed the 20-mers** (`ACCTATACTTCCAAGGCGCA` /
   `GGTGACTCTCTGTCTTGGCT` → 349 nt), the same ones DP3's order file carried. This **supersedes Erin's
   interim 2026-07-14 "19"**; whether or not the length caused a DP3 issue, 20 is the spec, and 349 sits
   one base under Twist's next price tier at 350. Pinned explicitly (`ADAPTER=` in
   `encode_step2_select.sbatch`); the earlier 19-mer smoke test is therefore the wrong length and was
   superseded by the full run. See memory `oligo-adapter-trap`.
