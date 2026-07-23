# DP4 library — components, selection, and status

Reference for the DP4 PepSeq library: what each component is, how designs were selected, and where
each file lives. Manuscript counterpart: `manuscript/sections/dp4_library.tex` (`sec:dp4`). Related:
`docs/CASE_ENCODING.md`, `docs/CYLINDER_PARAMS.md`.

Status (2026-07-23): **built, encoded, gated, and synthesis-ready.** `data/libraries/dp4_library.csv` is
the combined library: **36,000 rows** = **14,241 episcaf** constructs (the seven components C1–C6 + 8VDL)
plus **21,759 LX PfEMP1/EPCR minibinders** (`category=minibinder`, a separate de-novo binder arm folded in
so the whole DP4 library is one file). Each row is a 103-mer with a unique `library_member` and
`design_ID`. The episcaf shipping depth is top-20 per group for C1/C2 and top-10 per window for C3;
C4/C5/C6 are fixed size and 8VDL is top-29 per run. Selection ran under the soft-gate scorer
(`antibody_softgate`, epitope-PAE midpoint 2.5). The library was oligo-encoded and gated into
`data/libraries/dp4_order_file.csv` — **36,000 oligos, all verified** — the file that goes to Twist, with
`data/libraries/dp4_quote_file.csv` as the 2-column vendor-quote view of the same molecules. (The library
was assembled at 37,083, culled to 35,962 on 2026-07-21, then topped up to 36,000 on 2026-07-23 — see
**Cull to 35,962** and **Top-up to 36,000** below.)

> **Two things to know reading the numbers.** (1) **14,241 vs 36,000:** the episcaf-scaffolded portion is
> 14,241 — what the scoring and selection below concern; the encoding and the order file cover the whole
> **36,000** (episcaf + minibinders encoded together into one PepSeq assay — the lab runs one pooled oligo
> library per assay, so two projects on the same antigen share one synthesis order).
> (2) **Columns:** the library carries the **full 33-column** set (not the lean 5): 8 identity + 8 metrics
> (RMSDs, the PAE decomposition epitope/scaffold/mean, ptm, both clash flavors) + 4 scoring (`composite`,
> `rank_in_group`, `is_global_pass`, `island_index`) + 13 `lx_*` minibinder metrics. Minibinder rows are
> blank in the 12 episcaf metric/scoring columns (never scored on our axes) and carry the `lx_*`; episcaf
> rows are blank in `lx_*`. See the column dictionary below for which columns are blank where, and why.
> Minibinders added by `dp4_8vdl/scripts/08_add_minibinders.py` (idempotent; run after assembly);
> `stage07_named_peptides.py --exclude-category minibinder` would revert the encode to episcaf-only.

## Cull to 35,962 (2026-07-21)

After John's QC of the 37,083 library, four changes brought it to **35,962** (target ~36k). They are
applied in `scripts/stage06_assemble.py` (except the C6 depth, below), so the shipped file is regenerable:

1. **C6 controls → top-15 per antibody** (3,100 → 2,325). C6 was the single biggest block for what are
   controls; scaling the controls to the top-15 ranked designs per antibody keeps all three flavors
   (island1→A, island2→A, scaffold-6mer) and recovers ~775 slots. **How, and why it matters for
   reproducibility:** we do *not* re-run `build_c6_mutants` at top-15 — the scaffold-disruption control
   places its hexamers with a global RNG, so re-running on a different input set **re-randomizes** those
   sequences. That would have been scientifically fine (any scaffold disruption is a valid control), but
   the new sequences would then differ from the ones the oligo encoder already encoded, forcing a full
   re-encode. Instead we **filter the already-built top-20 C6 down to the top-15 base designs**
   (`scripts/filter_c6_depth.py`: strip each C6 `design_ID` to its base C1 `predID`, join `rank_in_group`,
   keep ≤15), preserving the exact encoded sequences. The encoded top-20 C6 is kept as
   `results/dp4_C6_controls.full.csv` so the filter is reproducible from committed files:
   `python scripts/filter_c6_depth.py --in results/dp4_C6_controls.full.csv --depth 15 --out results/dp4_C6_controls.csv`.
2. **Dedup (286).** 286 peptides were picked by more than one component — a C1 top-ranked design also
   drawn by C5's farthest-point metric-space sample of the same pool (285), plus 2 tiled-30mer collisions.
   The assembler drops duplicate `sequence`s keeping the first (C1 over C5); the C5 point is still assayed,
   just under the C1 barcode. So each peptide is ordered once.
3. **Cull the 60 no-accessibility designs.** Three C2 single-island targets (`5eu7_1P`, `5fhx_1P`,
   `6ztr_0P`) have a **2-residue island**, which gives too few epitope-Cα pairs (<3) to define the
   design→native rigid-body fit, so *both* the real Fab clash and the cylinder surrogate are uncomputable
   (`stage05` `too_few_epitope_pairs`). Their top-20 (60 designs) shipped with zero accessibility credit;
   the assembler now drops any design with `af3_clash_status == too_few_epitope_pairs`, so every shipped
   C1/C2 design carries accessibility.
4. **Fixed the `design_ID` naming.** The predID doubled the RFdiffusion backbone name (an AF3-output-dir
   filename artifact: `<name>_<name>_<d>_model…`), for both the C1/C2/C5 and 8VDL patterns. The assembler
   now collapses the repeated block; cosmetic, sequences/metrics untouched. The same collapse rule is kept
   in `stage06_superset.py`/`extend_superset.py` so the superset's `selected` join still matches.

**Order file:** regenerated at 35,962 from the *existing* encodings — no re-encode. Because the cull only
removes peptides and renumbers IDs, every surviving sequence was already encoded, so
`stage07_order_file.py --by-sequence` matches each oligo to the renumbered library by translated sequence
and drops the 1,121 no-longer-shipped encodings. Gate re-passed (349 nt, 20-mer adapters, round-trip).

**Composition after the cull (35,962):** C1 1,120 · C2 1,600 · C3 4,390 · C4 2,033 · C5 2,715 · C6 2,325 ·
8VDL 20 · minibinder 21,759. (Superseded by the top-up below.)

## Top-up to 36,000 (2026-07-23)

The cull left the library 38 short of the 36,000 target, so those slots went to the 8VDL arm (John,
2026-07-23: *"if you have a little leftover space, I'd suggest adding the next n best designs"*). Both
8VDL runs went from **top-10 to top-29 per run** — **+19 `8VDL_epitope` and +19 `8VDL_hotspots`**, 20 → 58
designs — landing the library exactly on **36,000**.

**The depth 29 is arithmetic, not a threshold.** 36,000 − 35,962 = 38 spare slots; split evenly over the
two 8VDL runs that is 19 each; 10 + 19 = 29. Nothing in the score distribution justifies a cut at 29 over
28 or 30 — it is simply how far down the existing ranking the leftover space reached. This is worth being
explicit about because 29 otherwise looks like a tuned parameter, and it is not one: the ranking, scorer,
and every other component's depth are untouched.

**No cluster job was needed to select them.** `07_consolidate.py --metrics-out` had already dumped every
8VDL design's scored metrics (`results/dp4_8vdl_top10_allmetrics.csv`, all 1,280), and a deeper top-k is a
deterministic re-slice of that dump: `rank_in_group <= topk` is exactly `nlargest(topk, composite)`. The
new `--from-metrics` mode does that re-rank with no AF3 re-read, reproducing a fresh `--topk 29` run:

```bash
python dp4_8vdl/scripts/07_consolidate.py --runs epitope,hotspots --topk 29 \
    --from-metrics results/dp4_8vdl_top10_allmetrics.csv --out results/dp4_8vdl_top29.csv
python scripts/stage06_assemble.py --depth 20 --c3-depth 10      # --vdl defaults to the top-29 file
python dp4_8vdl/scripts/08_add_minibinders.py --lx dp4_8vdl/data/LX_20260626.csv
```

Verified: top-29 is a strict superset of the shipped top-10 (nothing already ordered was displaced), no
misfolds entered (the epitope arm's rank-36 design — epitope RMSD 0.76 Å but *overall* RMSD 19.3 Å — sits
below the cut), and the assembled file has unique sequences and `design_ID`s with contiguous
`library_member` 1..36,000. Quality of the additions: the epitope arm spans epitope RMSD 0.97–2.16 Å at
clash 0–3; the hotspots arm has near-perfect epitope geometry (0.03–0.15 Å) but clashes 39–72, which is
the known generation-limited behaviour of that arm, not a selection artifact.

**Order file:** only the 38 new peptides needed encoding — `scripts/stage07_new_peptides.py` emits exactly
those, they were encoded on the cluster, merged with the whole-library encodings, and
`stage07_order_file.py --by-sequence` rebuilt the order file at **36,000** (37,121 encodings in, 1,121
dropped as no-longer-in-library). Matching is by *sequence* because a top-up renumbers `library_member`
for everything after the insertion point — adding 38 8VDL rows shifted all 21,759 minibinders.

> **Encoder input must be LF, not CRLF.** The first attempt at these 38 failed silently: the encoder
> rejected every line (`Processed 0 lines`) yet still printed `done` and wrote a 0-byte result, because the
> input had `\r\n` endings (`csv.writer`'s default dialect). `stage07_new_peptides.py` now writes bare
> `\n` and hard-fails on any CR. Always check step 1 reports `Processed <N> lines` before running step 2.

**Shipped composition (36,000):** C1 1,120 · C2 1,600 · C3 4,390 · C4 2,033 · C5 2,715 · C6 2,325 ·
8VDL 58 · minibinder 21,759.

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
  one file, 36,000 rows, in the 8-column PepSeq annotated format **plus the full metric + scoring set and
  the `lx_` minibinder columns** — 33 columns, schema below. This is the
  file to hand off. `designedSequence` is the full 103-mer in EPITOPEscaffold casing (epitope uppercase /
  scaffold lowercase) for every row; `sequence` is the plain uppercase 103-mer that gets synthesized.
- **`data/libraries/dp4_order_file.csv`** — the Twist synthesis order file. 36,000 oligos, two columns
  (`Seq ID`, `nucleotide_encoding_with_twist_adapters`). Every row verified by
  `scripts/stage07_order_file.py`: 349 nt, the 20-mer adapters on both ends, and a core that translates
  back to exactly its own peptide. **This is what goes to Twist.**
- **`data/libraries/dp4_quote_file.csv`** — the vendor **quote** view of that same order: 36,000 rows of
  `name,sequence` with flat zero-padded names (`DP4_00001`…`DP4_36000`) and the identical 349-mers.
  Vendors quote off a plain two-column list and the sequences need not be final, so this can go out while
  checking continues. Emitted by `scripts/stage07_quote_file.py`, which re-shapes the order file and never
  re-derives a sequence — quote and order are the same molecules by construction.
- **`data/libraries/dp4_superset.csv.gz`** — the all-designs superset (committed gzipped, 357,789 rows
  across all arms; see below). Not a deliverable — it exists for looking at the distributions the library
  was drawn from, and it is a true superset (the shipped library is a strict subset of it).
  **Stale as of the 2026-07-23 top-up:** its `selected` flag still marks 28,889 and needs a rebuild to
  reach 28,927. The row set is unaffected — the 38 added designs were already candidate rows in it.
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
| C2 | known-Ab, single island | best per (mAb, island); 87 island contigs | ranked, top-20 per island | 1,600 | `results/dp4_C2_single_island_ranked.top20.csv` |
| C3 | polyclonal 12-mer tiles | best per window, no antibody | ranked, top-10 per window † | 4,390 | `results/dp4_C3_12mer_ranked.top20.csv` |
| C4 | linear 30-mer controls | bare tiled peptides, no scaffold | exhaustive tiling | 2,033 | `data/libraries/dp4_tiled30mers_fasta.csv` |
| C5 | metric-space titration | designs spread across metrics | farthest-point sample | 2,715 § | `results/dp4_C5_titration.csv` |
| C6 | scaffolded-epitope controls | island→Ala + scaffold disruption | C1 top-15 base × flavors | 2,325 ‡§ | `results/dp4_C6_controls.csv` |
| 8VDL | PfEMP1 conserved epitope | two epitope definitions | ranked, top-29 each ¶ | 58 | `results/dp4_8vdl_top29.csv` |
| | | | | **14,241 episcaf** | `data/libraries/dp4_library.csv` (episcaf portion) |
| LX | PfEMP1/EPCR minibinders | de-novo binders (not scaffolded) | LX-filter passers | 21,759 | `dp4_8vdl/scripts/08_add_minibinders.py` (source gitignored) |
| | | | | **36,000 total** | `data/libraries/dp4_library.csv` (whole file) |

*(Shipped counts, after the 2026-07-21 cull and the 2026-07-23 top-up. C2 is 1,600 not 1,660 — 60
no-accessibility 2-residue-island designs culled; C4 2,033 and C5 2,715 net the dedup; C6 2,325 is
top-15. See **Cull to 35,962** and **Top-up to 36,000** above.)*

¶ **Why 29, specifically?** It is a budget number, not a quality threshold. The cull left the library at
35,962, i.e. 38 slots under the 36,000 target, and those spare slots were given to the 8VDL arm; split
evenly across its two runs that is +19 each, taking both from the original top-10 to **top-29**. Nothing
about the ranking changes at 29 — read it as "top-10 plus however many the leftover space allowed,"
which is why the depth is not a round number like the other components'.

† C3 is shipped deep (top-10). Its windows are **12-mers stepping by 2 residues**, so neighbouring tiles
are highly redundant (adjacent windows share 10 of 12 residues) — which argues for a shallow cut. We
ship top-10 anyway because that redundancy does not protect against the failure mode here: C3 has the
weakest clash distribution of any component, so the per-design success probability is low and
overlapping tiles can fail together. It is also the arm with the most to gain (a few hits in the
no-antibody setting would be the first evidence the approach works without a known antibody), so the
spare budget capacity is spent here (John, 2026-07-14). *(Not to be confused with C4, which tiles
**30-mers at step 6**.)*

‡ C6 is derived from the C1 top-15 base (island1→Ala + island2→Ala for dual-island epitopes + scaffold
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

> **Midpoint history (2026-07-20):** the epitope-PAE midpoint was originally 5 (borrowed from the global
> `mean_pae < 5` threshold), then **retuned to 2.5**, set from the data (see below). That swapped ~107 of
> the 1120 C1 selections (~10%), counts and components unchanged. The library was **re-selected,
> re-assembled, and re-encoded at 2.5**, so the shipped `dp4_library.csv` and order file already reflect
> 2.5 — the code scorer and the shipped file agree.

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
**28,927 ship** (C1 1,120 + C2 1,600 + C3 4,390 + 8VDL 58 + LX 21,759) — the library's 36,000 minus the
C4/C5/C6 controls, which aren't candidate-pool designs — and the four-filter passers are C1 727 + C2 407
(C3/minibinders have no antibody-based global pass defined).

> **The committed `.gz` is one rebuild behind.** It still flags **28,889** selected, from before the
> 2026-07-23 top-up. Only the flag is stale: all 357,789 rows are correct, and the 38 added 8VDL designs
> were already candidate rows in it — the top-up promoted them, it did not introduce them. Rerun
> `sbatch scripts/build_superset.sbatch` to refresh (analysis artifact, not on the deliverable path).

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

Two epitope definitions, top-29 each (58 total; originally top-10/20 — the extra 19 each came from the
2026-07-23 top-up that spent the leftover order slots, see **Top-up to 36,000**):
- `epitope` — the contiguous contact window C652–C673 (22 residues, covering all 13 residues with a
  heavy atom within 4 Å of the Fab). This is the strong constraint, the analog of C1.
- `hotspots` — only F655, F656, and E666, fixed at their native crystal coordinates, with the design
  building around them (a minimal hotspot graft).

Scoring is by `07_consolidate.py`, which aligns each predicted epitope onto the native chain-C frame,
counts the real H/L clash there, and (2026-07-20) ranks under the **same `antibody_softgate` scorer as
C1/C2** — real Fab clash, no cylinder surrogate. The two definitions separate cleanly, and the hotspot
one fails informatively:
- **epitope** (contiguous C652–673): shipped top-29 spans epitope RMSD 0.97–2.16 Å at **0–3** clashing
  residues (the top-10 was 0.97–1.36 Å at 1–2) — the soft-gate found accessible designs, and the extra
  19 stay accessible, only loosening on epitope RMSD.
- **hotspots** (F655/F656/E666): shipped top-29 spans epitope RMSD 0.03–0.15 Å but **39–72** clashing
  residues (top-10 was 0.04–0.11 Å at 39–71). This is **not** a selection artifact — across all **640**
  hotspot designs the clash floor is **14** and none fall below 10, so an accessible minimal-hotspot graft
  was never generated (generation-limited, the same phenomenon as C1). The three discontiguous residues
  appear to force the scaffold into the Fab footprint. We ship these anyway as a **documented negative
  result** — a limitation of the minimal-graft approach, not a scorer failure. Going deeper to 29 does not
  change that conclusion; it adds more of the same clashing designs.
  (Per-design metrics for all 1,280: `results/dp4_8vdl_top10_allmetrics.csv`.)

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

**Which columns are populated for which component** is given in full by the 33-column dictionary above
(measured fill matrix + the reason for every blank). In short: C1/C2/C5 carry both clash flavors; C3 has
no `af3_clashes` (no antibody) but has the cylinder; 8VDL has the real clash but no cylinder; C4 (linear
tiles) and C6 (mutants) never went through AF3, so all fold metrics are blank; only C2 has `island_index`;
C5 has no `composite`/`rank` (it is a metric-space sample). Blank cells are honest gaps, not missing work.

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
python episcaf_pipeline/build_whole_epitope_designs.py --drop-targets 2h32,4xwo,7a3t \
  --out results/whole_epitope_designs.parquet   # -> 2,206 contigs; also writes the .csv sibling the
                                                #    case-encoder reads (2,206 x 8 x 8 = 141,184 designs
                                                #    generated; 140,716 land status==ok after AF3)
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

# C6 (local; seeded). Build at top-20 -> keep as the full/encoded reference, then FILTER to top-15
# (do NOT re-run build_c6 at 15 -- it re-randomizes the scaffold-disruption controls; see Cull to 35,962).
python episcaf_pipeline/scaffolded_epitope_controls/build_c6_mutants.py \
  --input results/dp4_C1_scaffoldEPITOPE.csv \
  --id-col token --target-col target --seq-col scaffoldEPITOPE \
  --drop-targets 2h32,4xwo,7a3t --out results/dp4_C6_controls.full.csv
python scripts/filter_c6_depth.py --in results/dp4_C6_controls.full.csv --depth 15 \
  --out results/dp4_C6_controls.csv                         # top-15 (2,325), encoded sequences preserved

# 8VDL arm (Gemini: RFD3->MPNN->AF3; then consolidate per definition under antibody_softgate). The
# --metrics-out dump this writes is what makes later depth changes free -- see the re-rank below.
python dp4_8vdl/scripts/07_consolidate.py --runs epitope,hotspots --topk 10 --out results/dp4_8vdl_top10.csv

# Shipped 8VDL depth is top-29/run (the 2026-07-23 top-up). This is a LOCAL deterministic re-slice of the
# metrics dump above -- no AF3 re-read, no cluster job (rank_in_group <= topk == nlargest(topk, composite)).
python dp4_8vdl/scripts/07_consolidate.py --runs epitope,hotspots --topk 29 \
  --from-metrics results/dp4_8vdl_top10_allmetrics.csv --out results/dp4_8vdl_top29.csv

# Assemble the episcaf library (local) -> data/libraries/dp4_library.csv, 14,241 episcaf constructs
python scripts/stage06_assemble.py --depth 20   # C1/C2 top-20; C3 top-10; --vdl defaults to the top-29 file

# Fold in the passing LX minibinders -> 36,000 rows (needs the LX source; see the minibinder arm above)
python dp4_8vdl/scripts/08_add_minibinders.py --lx dp4_8vdl/data/LX_20260626.csv

# Export the oligo-encoder input (whole 36,000-row library) -> data/libraries/dp4_named_peptides.csv
python scripts/stage07_named_peptides.py \
  --library data/libraries/dp4_library.csv --out data/libraries/dp4_named_peptides.csv
  # add `--exclude-category minibinder` to encode only the 14,241 episcaf portion instead

# Oligo-encode the library (Gemini). Run BOTH steps in the same working dir; step 1 is the long pole
# (C++ sampler, hours). ADAPTER is pinned to the 20-mers inside encode_step2_select.sbatch -- do NOT
# rely on the encoder's own --adapter default, which is the 19-mer form (see oligo-adapter-trap).
# Run in a FRESH rundir under the DURABLE $WS, not $REPO (/scratch sweeps at ~30 days, and reusing a
# rundir name is how the whole-library encodings ended up shadowed by an older file -- see the
# "Which DP4_best_encodings?" box below). Copy the input + weights in first.
mkdir -p $WS/runs/dp4_encoding_<tag> && cd $WS/runs/dp4_encoding_<tag>
cp $REPO/data/libraries/dp4_named_peptides.csv $REPO/episcaf_pipeline/oligo_encoding/codon_weights_updated.csv .
# TOOL_DIR is REQUIRED (no default, by design -- a bad default once ran an unexecutable binary). Point it
# at the built LadnerLab encoder; the env var must come BEFORE `sbatch`, not after (else it is read as the
# script path). Same TOOL_DIR for both steps. Step 2 pins the 20-mer ADAPTER internally -- do not pass it.
TOOL=/home/bneff/Library-Design/oligo_encoding
JID1=$(TOOL_DIR=$TOOL INPUT=dp4_named_peptides.csv sbatch --parsable --time=12:00:00 \
  $REPO/episcaf_pipeline/oligo_encoding/encode_step1_generate.sbatch)
TOOL_DIR=$TOOL sbatch --dependency=afterok:$JID1 \
  $REPO/episcaf_pipeline/oligo_encoding/encode_step2_select.sbatch

# TOP-UP ENCODE (only when the library GREW). Encode just the new peptides, never all 36,000:
python scripts/stage07_new_peptides.py \
  --peptides data/libraries/dp4_named_peptides.csv \
  --encoded  $WS/runs/dp4_encoding_full/DP4_best_encodings.wholelib_37083 \
  --out      data/libraries/dp4_named_peptides.new38.csv        # -> only the peptides lacking an encoding
# ...copy that file into a fresh rundir and run the same two steps. The input MUST be LF, never CRLF:
# the encoder rejects every CRLF line, prints "Processed 0 lines", then still says "done" and writes a
# 0-BYTE result. ALWAYS confirm step 1 before step 2:  grep Processed slurm-*.out  (want N, not 0).
# Then merge the new encodings under the whole-library ones (same header, strip the duplicate):
MERGED=$WS/DP4_best_encodings_36000
head -1 $WS/runs/dp4_encoding_full/DP4_best_encodings.wholelib_37083  > $MERGED
tail -n +2 $WS/runs/dp4_encoding_full/DP4_best_encodings.wholelib_37083 >> $MERGED
tail -n +2 $WS/runs/dp4_encoding_new38/DP4_best_encodings               >> $MERGED

# Emit + VERIFY the Twist order file (Gemini). Checks every row: 20-mer adapters, 349 nt, and that each
# core translates back to its own peptide. Writes nothing if any row fails.
python $REPO/scripts/stage07_order_file.py --by-sequence \
  --best-encodings $MERGED \
  --peptides       $REPO/data/libraries/dp4_named_peptides.csv \
  --out            $REPO/data/libraries/dp4_order_file.csv     # -> 36,000 oligos, all verified
  # --by-sequence: the cull and the top-up both renumber library_member (adding 38 8VDL rows shifted all
  # 21,759 minibinders), so oligos are matched to the new names by TRANSLATED SEQUENCE, not by name.
  # Encodings for peptides no longer in the library are simply dropped -- NO re-encode of survivors.
  # NB: the whole-library encodings are `DP4_best_encodings.wholelib_37083`. The plain
  # `DP4_best_encodings` in that same directory is the OLDER episcaf-only 15,324 file -- see below.

# Emit the 2-column vendor QUOTE file (local) -- same molecules as the order file, flat DP4_00001 names.
python scripts/stage07_quote_file.py \
  --order-file data/libraries/dp4_order_file.csv \
  --out        data/libraries/dp4_quote_file.csv               # -> 36,000 x 349mer, adapters verified

# ALL-DESIGNS SUPERSET (John's ask -- every candidate design, not just the selected ones). ONE cluster
# pass: build C1/C2/C3 -> extend with 8VDL + passing minibinders -> gzip. extend reads the minibinders
# from the committed dp4_library.csv, so NO external file is needed. Emits $WS/dp4_superset.csv (== the .gz).
sbatch scripts/build_superset.sbatch
git add data/libraries/dp4_superset.csv.gz && git commit -m "superset: rebuild" && git push
```

Scorer weights and transforms are config, not magic numbers (`episcaf_analysis/presets.py`). C5 and C6
are deterministic (FPS is seed-free deterministic; C6 seeds its RNG).

> ### Which `DP4_best_encodings`? (two files, same name — read this before using one)
> There are **two different encodings files with the same name, in directories with the same name.**
> Check size or line count, never the path alone:
>
> | file | size | what it is |
> |---|---|---|
> | `$WS/runs/dp4_encoding_full/DP4_best_encodings.wholelib_37083` | **31 MB** (37,084 lines) | **whole library, 37,083** — what every shipped oligo came from |
> | `$WS/runs/dp4_encoding_full/DP4_best_encodings` | 12.9 MB (15,325 lines) | the earlier **episcaf-only 15,324** encode, superseded |
> | `/scratch/.../runs/dp4_encoding_test/DP4_best_encodings` | 42 KB | 50-peptide smoke test, **wrong (19-mer) adapters** |
>
> Records are ~840 bytes, so **size alone identifies the file**: 15,324 → ~12.9 MB, 37,083 → ~31 MB.
>
> **How this happened (2026-07-23):** the whole-library encode ran on `/scratch` on 2026-07-20, *after*
> the migration had already copied `dp4_encoding_full` to `$WS` on 07-16. The durable copy was therefore
> the stale episcaf-only one while the real file sat on a filesystem that sweeps at ~30 days. Earlier
> versions of this doc claimed the whole-library encode lived in `$WS/runs/dp4_encoding_full/` — it did
> not. The file has since been rescued to `$WS` under the distinct `.wholelib_37083` name so the two can
> no longer be confused. Lesson: a rundir name is not an identity; re-running a step into a same-named
> directory can leave the durable copy pointing at the older artifact.
>
> **If encodings are ever lost:** `data/libraries/dp4_order_file.csv` is committed to git and holds every
> finished adapter-flanked oligo, so a library can be rebuilt from the order file plus encodings for only
> the *new* peptides. A full re-encode is never required.

## Budget and depth

DP4 targets a **36,000-row** library including all minibinders, which leaves roughly 10–15k slots for
Episcaf designs (`memory: dp4-budget`). Depth is set at top-20 for C1/C2 — the most the ranked files hold,
and the depth C6 was built at, so no C6 rebuild is needed. **C3 is top-10** (2026-07-14): John flagged that
the spare capacity (~2k slots to reach 36k) is best spent maximizing polyclonal hits, given C3's weak clash
distribution. As built this was **15,324** episcaf; the 2026-07-21 cull took it to 14,203 and the
2026-07-23 top-up returned it to **14,241** (C1 1,120, C2 1,600, C3 4,390, C4 2,033, C5 2,715, C6 2,325,
8VDL 58) + 21,759 minibinders = **36,000** total, exactly on target. See **Cull to 35,962** and
**Top-up to 36,000**.

The budget is now **fully spent** — the library sits exactly on 36,000, so any future addition has to
displace something. The 8VDL depth (top-29) is the marginal dial: it absorbed the last 38 slots and is
the natural place to give them back.

C3 depth is the one elastic dial left: top-3 = 1,317, top-5 = 2,195, top-10 = 4,390 (shipped). Set it
with `stage06_assemble.py --c3-depth <n>` if the final minibinder count moves the headroom.

## Build log (all steps complete)

*(This section is the completed build history, not open work — every item below is done. The whole
library, **36,000 rows**, is assembled, encoded, gated, and shipped as of 2026-07-23: assembled at 37,083,
culled to 35,962 on 07-21, topped up to 36,000 on 07-23.)*

0. C1 redo at 103 — done 2026-07-11. RFD3→MPNN→AF3 on the 2,206-contig ledger, metrics (140,716 designs,
   all `status==ok`), C1 re-selected to top-20 (1,120), re-case-encoded (`case_encode_whole_epitope.py`),
   C5 (3,000) and C6 (3,100) rebuilt on the new pool. C1/C5/C6 are native 103; the 104→103 trim is a no-op.
1. Assembly (`06_library`) — done 2026-07-13. `scripts/stage06_assemble.py --depth 20` concatenated the
   seven components into `data/libraries/dp4_library.csv` (15,324 episcaf constructs; the 21,759
   minibinders are folded in next by `08_add_minibinders.py` for 37,083), applying the 56-exclusion
   (C1/C2), the top-20 (C1/C2) and top-10 (C3) depth cuts, and global numbering. `library_member` and
   `design_ID` are unique and every sequence is 103 residues. Ships the full **33-column schema** (see the
   column dictionary): 8 annotation columns + the metric/scoring set + the 13 `lx_` minibinder columns,
   with `designedSequence` in EPITOPEscaffold casing on every scaffolded row.
2. Oligo encoding — **whole library encoded 2026-07-20** (37,083 = 15,324 episcaf + 21,759 minibinders;
   supersedes the 2026-07-16 episcaf-only 15,324 encode after the 2.5 re-selection + soft-gate 8VDL +
   minibinder addition). Encoder input: `scripts/stage07_named_peptides.py` →
   `data/libraries/dp4_named_peptides.csv` (`name,seq`, no header, all 103-mers; regenerate after any
   library change). Ran on Gemini in `/scratch/bneff/episcaf/runs/dp4_encoding_full/` — note **`/scratch`**,
   not `$WS`: this ran after the migration, so it was NOT captured by it, and the same-named `$WS` dir holds
   the earlier episcaf-only encode. Rescued to `$WS/.../DP4_best_encodings.wholelib_37083` on 2026-07-23
   (see the "Which `DP4_best_encodings`?" box). Encoder: the LadnerLab tool (step 1 sampler
   → step 2 NN selector, DP3 recipe + `codon_weights_updated.csv`; `episcaf_pipeline/oligo_encoding/`, see
   its README and manuscript `sec:oligo`) → `DP4_best_encodings`, all 37,083, nothing dropped. The order
   file is **not** a further encoding step — step 2 already emits the adapter-flanked oligo, so the order
   file is a two-column slice, emitted + verified by `scripts/stage07_order_file.py` →
   `data/libraries/dp4_order_file.csv`. All 37,083 verified: 349 nt, 20-mer adapters both ends, every
   core translates back to its own peptide, one encoding per peptide, none missing. *(Step 3 below then
   culls this to the shipped **35,962** by dropping the culled rows from the order file by sequence — no
   re-encode; the 37,083 encodings still exist, the order file is the 35,962 subset.)* **Twist-ready.**
   **Adapter length resolved 2026-07-16: John confirmed the 20-mers** (`ACCTATACTTCCAAGGCGCA` /
   `GGTGACTCTCTGTCTTGGCT` → 349 nt), the same ones DP3's order file carried. This **supersedes Erin's
   interim 2026-07-14 "19"**; whether or not the length caused a DP3 issue, 20 is the spec, and 349 sits
   one base under Twist's next price tier at 350. Pinned explicitly (`ADAPTER=` in
   `encode_step2_select.sbatch`); the earlier 19-mer smoke test is therefore the wrong length and was
   superseded by the full run. See memory `oligo-adapter-trap`.
3. Cull to 35,962 — **2026-07-21** (after John's QC). C6 → top-15 (filter, not rebuild), dedup the 286
   picked-twice, drop the 60 no-accessibility 2-residue-island C2 designs, and collapse the doubled
   `design_ID`. Order file regenerated at 35,962 by `stage07_order_file.py --by-sequence` from the
   existing encodings (no re-encode). Full detail in **Cull to 35,962** at the top of this doc.
4. Top-up to 36,000 — **2026-07-23**. The cull left 38 slots under target; they went to the 8VDL arm,
   top-10 → **top-29 per run** (+19 epitope, +19 hotspots; 20 → 58 designs). Selection was a local
   deterministic re-slice of the existing metrics dump (`07_consolidate.py --from-metrics`) — no cluster
   job. Only the 38 new peptides were encoded (`stage07_new_peptides.py` → a fresh rundir); merged with
   the whole-library encodings and the order file rebuilt at **36,000**, gate re-passed (349 nt, 20-mer
   adapters, every core round-trips). Added `data/libraries/dp4_quote_file.csv` for the vendor quote.
   Two traps surfaced and are written up above: the **CRLF encoder input** (first attempt silently
   encoded nothing) and the **two same-named `DP4_best_encodings` files** (the whole-library encode was
   on `/scratch`, not `$WS`, and has now been rescued). Full detail in **Top-up to 36,000**.
   *Outstanding:* the superset `.gz` still flags 28,889 selected and needs a rebuild to 28,927.
