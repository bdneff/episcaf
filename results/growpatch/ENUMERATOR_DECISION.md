# growPatch enumerator: benchmark + decision

**Settled 2026-11-13.** What the no-antibody epitope enumerator should be, decided on data.

## The question

Given an antigen with no known antibody, enumerate candidate epitope patches so that whatever the
real epitope is, one candidate covers it — then scaffold the candidates. Linear sequence tiling can
only propose contiguous stretches, so it is blind to **conformational** epitopes (residues close in
3D, scattered in sequence), which are most real antibody epitopes. The enumeration has to work in 3D.

### This is not epitope *prediction*

growPatch answers a different question than DiscoTope, BepiPred, SEMA, or any B-cell epitope
predictor. Those rank the surface by how *likely* each residue is to be an epitope and return the
most probable one(s) — they optimize **precision**, and they are wrong often enough that betting one
scaffolding campaign on their top call is risky. We want the opposite: an **exhaustive, unbiased
enumeration of every surface region that could plausibly be an antibody footprint**, so we can
scaffold them all and let the experiment decide. The objective is **recall and coverage of the
surface**, not a probability ranking. A learned predictor is therefore not a competitor to the
enumerator; at most it is an optional *ranking layer on top* of the exhaustive list, for when compute
forces us to prioritize — and even then the enumeration stays complete underneath. This is why the
benchmark metric below is coverage/recall, and why the shipped default unions two strategies rather
than pruning to a favorite.

Two strategies, both ported here:

- **growpatch** — seed-and-grow along the residue contact graph (John Altin's `segmentation/growPatch.py`;
  `episcaf_analysis/growpatch/core.py` + `read.py`). From every surface residue, accrete the nearest
  spatially-connected surface residue to ~`patch_size`, then gap-fill / island-trim / NMS.
- **eyeview** — antibody's-eye-view direction caps (`episcaf_analysis/growpatch/eyeview.py`). Sample
  approach directions over the sphere; each direction's candidate is the surface cap an antibody from
  there would contact. Same cleanup as growpatch, so only the enumeration geometry differs.

## Benchmark

`scripts/growpatch_benchmark.py` on `reference/sabdb_epitopes_4A_sample200_chain1.csv` (SAbDab,
single-chain epitopes; vendored here as `sample200_chain1.csv`). Kept the **conformational** ones —
gap-10 island count in {2,3}: 2+ islands = discontinuous in sequence, ≤3 = inside the enumerators'
≤3-island footprint design. The 1-island (continuous) epitopes are a separate control. **117
conformational + 19 continuous** benchmarked (1 skipped: `9shd`, chain-id mismatch).

Metric per epitope, per method: **best_single_coverage** = max over candidates of
|candidate ∩ true epitope| / |true epitope| (true = nonGapFilled 4 Å contacts present with a complete
backbone; both methods get it as `protect`, John's convention — this asks "is a good candidate in the
pool"). Also candidate count (fewer is better at equal coverage) and best-patch size (coverage must
be read against size — a bigger patch covers more trivially).

Reproduce:

    python scripts/growpatch_benchmark.py \
        --structures ~/Desktop/projects/epitope-universe/data/structures \
        --radii 11,12,13,14 --islands 2,3 --control-islands 1 \
        --out results/growpatch/benchmark_coverage.csv

(structures = the gitignored local SAbDab download, `<pdb>.cif.gz`.)

### Conformational epitopes (n=117, gap-10 islands 2–3)

| method | mean cov | median | recall ≥0.8 | mean #cand | mean best-patch size |
|---|---|---|---|---|---|
| growpatch | 0.899 | 1.000 | 0.795 | 55.7 | 31.7 |
| eyeview r=11 | 0.723 | 0.727 | 0.359 | 30.7 | 21.5 |
| eyeview r=12 | 0.793 | 0.786 | 0.496 | 37.5 | 26.1 |
| eyeview r=13 | 0.858 | 0.909 | 0.718 | 39.8 | 30.8 |
| eyeview r=14 | 0.895 | 0.923 | 0.821 | 38.2 | 35.7 |
| **union (growpatch ∪ eyeview r=14)** | **0.940** | **1.000** | **0.923** | ~94 | — |

### Continuous control (n=19, 1 island)

| method | mean cov | recall ≥0.8 | mean #cand |
|---|---|---|---|
| growpatch | 0.975 | 0.947 | 45.7 |
| eyeview r=14 | 0.933 | 0.947 | 27.3 |

## What the numbers say

1. **Head-to-head is a tie.** At matched patch size (~31 res: growpatch 31.7 vs eyeview r=13 30.8),
   growpatch wins (0.899 vs 0.858 coverage, 0.795 vs 0.718 recall@0.8). eyeview only matches
   growpatch's coverage by growing bigger patches (r=14, size 35.7), but does it with ~30% fewer
   candidates and marginally higher recall. Neither dominates per residue.
2. **They are complementary — the real result.** At ≥0.8 coverage on the 117 conformational epitopes:
   both methods get 81, **growpatch-only 12, eyeview-only 15**, neither 9. Each recovers ~12–15
   epitopes the other misses. Unioning the candidate sets lifts recall@0.8 from ~0.80 (either alone)
   to **0.923** and mean coverage to 0.94, for ~94 candidates/antigen.
3. The 9 "neither" epitopes (far-apart multi-island) are a real ceiling of a ≤3-island footprint.

## Decision

**Ship both, unioned.** `read.enumerate_candidates(cif, chains, methods=("growpatch","eyeview"), ...)`
runs both and cross-dedups (patches both find are tagged as co-found — a confidence signal). growpatch
is the validated base (byte-parity with John's code); eyeview is the complementary second strategy.
For DP5's comprehensive surface search, run the union — it is the ~92%-recall candidate set that
realizes the "don't miss the real epitope" goal. This follows the project's discipline: we did not
swap on opinion; we measured, and the data says *both*.

## User knobs (all parameters of `enumerate_candidates`)

| knob | controls | default |
|---|---|---|
| `sasa_cut` | surface cutoff (relSASA ≥) | 0.2 |
| `sasa_scope` / `context` | per-chain vs antigen-complex SASA; extra buried chains | antigen |
| `patch_size`, `grow_link` | how far a growpatch seed grows | 16, 5.0 Å |
| `footprint_radius`, `depth` | eyeview cap size | 14 Å, = radius |
| `gap` | gap-fill size (bridge sequence gaps ≤) | 10 |
| `min_res`, `max_islands`, `min_island_len` | footprint-plausibility filter | 12, 3, 3 |
| `nms_jaccard` | **exhaustiveness** — lower = less dedup = more, overlapping candidates | 0.8 |
| `top_k` | hard cap on candidates returned | none |
| `n_dirs` | eyeview approach directions sampled | 150 |

## Caveats

- Coverage is read against patch size; the benchmark protects the known epitope (measures pool
  quality, the right question for enumeration). Margins would shift on multi-chain epitopes (not yet
  benchmarked) and on our own antigens.
- Patches run large (~31 res from gap-fill); tightening the footprint is a separate knob
  (`gap`, `min_res`), orthogonal to the enumerator choice.
- eyeview is new (not independently validated like growpatch); the benchmark above is its validation.
