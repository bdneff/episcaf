# DP4 — 8VDL PfEMP1 conserved-epitope scaffolding

Self-contained arm that scaffolds the conserved EPCR-binding epitope of **PfEMP1 CIDRα1.4**
(structure **8VDL**) and folds the top designs into the DP4 PepSeq library. It is kept as its own
subdirectory with its own Gemini job because it belongs to the sibling *minibinder / PfEMP1* project
(`~/Desktop/projects/minibinder`); it rides into DP4 only for synthesis, and is **consolidated from
the repo root at the end** (score under `antibody_softgate` → top 10 per run → the full metric/scoring
column set → `stage06_assemble`).

## The epitope region — the conserved EPCR-binding surface

```
Target : 8VDL  (HB3VAR03 CIDRα1.4 + C7 Fab; Reyes et al., Nature 636:182–189, 2024)
Chain  : C   (the CIDR antigen; chains H/L are the C7 antibody)
Epitope: residues 652–673  (22 aa, one contiguous island) — the scaffolded window (see "On the window" below)
Hotspots (paper): F655, F656, E666
```

PfEMP1 is the hypervariable *P. falciparum* surface antigen that drives severe malaria and escapes
immunity by antigenic variation. But the residues its CIDRα1 domain uses to grip the host receptor
**EPCR cannot vary freely** — so a binder or antibody aimed at that conserved contact surface could
recognize **many PfEMP1 variants at once**. In 8VDL that conserved, EPCR-binding contact surface is the
single contiguous stretch **chain C 652–673**; the 2024 paper singles out **F655/F656/E666** as the
functional hotspots within it. We take this **22-residue contiguous window as the epitope region**
because it is the conserved, functionally constrained EPCR-binding surface the whole broadly-acting
strategy is built on — and because it is contiguous and fully resolved in the crystal, so it scaffolds
as one clean island with no gap-bridging guesswork. (An earlier **651–670** framing was a 20-mer guess
that missed the contact residue K673; 652–673 is the corrected window — see "On the window" below.)
A Rosetta predecessor
(6SNY) scaffolded this epitope and recapitulated EPCR-binding geometry but elicited weaker
adhesion-inhibitory antibodies than the native protein; RFD3 is not constrained to a particular fold,
so it explores a broader scaffold space for presenting it.

Because chains H/L are the C7 antibody, this is a **known-antibody run** — the real AF3 clash filter
(`af3_n_clash_res`) applies, exactly as for DP4 components C1/C2.

## Two runs: full epitope vs. minimal hotspots

We generate two independent contig sets and score each on its own, taking the **top 10** from each:

| run | fixed motif | what it asks |
|---|---|---|
| **`epitope`** | the contiguous window **C652–673** (22 aa, 1 island) | present the *entire* conserved epitope in native geometry — strongest constraint, the direct analog of C1 |
| **`hotspots`** | only **F655, F656, E666** (2 islands) | present just the *functionally critical* residues, letting RFD3 design everything else — minimal constraint / hotspot graft |

Both are scaffolded into a constant **103-mer** (the PepSeq maximum); the fixed atoms carry the residues'
native crystal coordinates, so each run preserves the motif's 3-D arrangement. Comparing them in the
assay asks whether the minimal hotspots suffice, or the full epitope is needed.

**On the window 652–673.** The AbDb/IEDB 4 Å contact epitope — the residues within a heavy atom of the
C7 Fab (`contact_epitope.py`: 652,653,655–657,659–661,666,667,669,670,673) — all fall inside the
contiguous stretch 652–673. Rather than scaffold that footprint as a fragmented six-island motif (hard
to hold rigidly), we take the single contiguous 652–673 window, which presents every contact residue
(including **K673**) and scaffolds cleanly as one island. This window supersedes both the earlier
651–670 guess (which missed K673) and the six-island contact target — they are one and the same epitope.

## Pipeline (self-contained; runs on Gemini)

Same RFD3→ProteinMPNN→AlphaFold3 shape as the rest of episcaf, but chain-C / real-residue aware
(adapted from `minibinder/pfemp1/run_pfemp1_scaffold`). All commands from `dp4_8vdl/`.

```bash
# 1. contigs — one CSV per run (--n-contigs is John's "n designs" knob; 10 each → 10×8×8=640 designs)
python scripts/01_generate_contigs.py --target epitope  --n-contigs 10 --out 01_contigs/epitope.csv
python scripts/01_generate_contigs.py --target hotspots --n-contigs 10 --out 01_contigs/hotspots.csv

# 2. RFD3 input JSONs
python scripts/02_emit_rfd3_inputs.py --contigs_csv 01_contigs/epitope.csv --out_dir 02_rfd3/epitope/inputs
python scripts/02_emit_rfd3_inputs.py --contigs_csv 01_contigs/hotspots.csv  --out_dir 02_rfd3/hotspots/inputs

# 3. RFD3 (Gemini GPU) — one task per contig, 8 backbones each
sbatch --array=1-10%200 scripts/03_rfd3_array.sbatch 02_rfd3/epitope
sbatch --array=1-10%200 scripts/03_rfd3_array.sbatch 02_rfd3/hotspots

# 4. FIXED backbone PDBs — the ONE 8VDL-specific step (multi-island FIXED positions), then
#    everything downstream reuses episcaf's proven stage03/stage04 (run per-run, from repo root):
python dp4_8vdl/scripts/04_make_fixed_pdbs.py \
    --contigs_csv dp4_8vdl/01_contigs/epitope.csv \
    --rfd3_outputs_dir dp4_8vdl/02_rfd3/epitope/outputs \
    --out_dir dp4_8vdl/runs/epitope/03_mpnn/fixed_pdbs

# 5. MPNN (episcaf, proven) — 8 seqs/backbone
python scripts/stage03_mpnn_submit.py --fixed_pdb_dir dp4_8vdl/runs/epitope/03_mpnn/fixed_pdbs \
    --outdir dp4_8vdl/runs/epitope/03_mpnn/mpnn_pdbs --batch_size 300 --tag 8vdl_ep20

# 6. AF3 (episcaf, proven)
python scripts/stage04_af3_emit_jsons.py --mpnn_pdb_dir dp4_8vdl/runs/epitope/03_mpnn/mpnn_pdbs \
    --out_dir dp4_8vdl/runs/epitope/04_af3/inputs --seed 1
sbatch --array=1-N scripts/stage04_af3_array.sbatch dp4_8vdl/runs/epitope   # N = ceil(#jsons/100)
#   (repeat 4–6 for hotspots)

# 7. consolidate — rank under antibody_softgate (real H/L clash), top 10 per run, full metric/scoring set
python dp4_8vdl/scripts/07_consolidate.py --runs epitope,hotspots --topk 10 \
    --out results/dp4_8vdl_top10.csv   # also writes *_allmetrics.csv (all 1,280 candidates)
```

`data/8VDL.pdb` is the crystal (chains C/H/L), committed here so the arm is self-contained. The final
`results/dp4_8vdl_top10.csv` is what `stage06_assemble` folds into the library (as its own target).

## Notes
- The RFD3/AF3 `sbatch` are adapted from the minibinder project; the stale hardcoded `--chdir` has
  been removed (they inherit the submit dir), but dry-run once on the current cluster before the full array.
- `--n-contigs` × 8 RFD3 × 8 MPNN designs per run; top 10 come out, so 10 contigs is ample. Note
  `epitope` (one island, 2 scaffold segments) caps at 64 unique contigs; `hotspots` (3 segments)
  has no practical cap. RFD3 still makes 8 stochastic backbones per contig, so 10 contigs = 640 designs.
- MPNN/AF3 are **not ported** — the only 8VDL-specific step is `04_make_fixed_pdbs.py` (computes the
  multi-island FIXED positions and reuses episcaf's `cif_to_fixed_pdb`); steps 5–6 are episcaf's
  just-run `stage03_mpnn_submit.py` / `stage04_*`. The prior art's `04/05/06` are superseded and not
  carried in.
- **`07_consolidate.py` (done):** scores the AF3 outputs under `antibody_softgate` (real H/L clash, not
  the cylinder), ranks top 10 per run, and emits the full metric/scoring set — epitope/scaffold/mean PAE,
  ptm, composite, rank_in_group, is_global_pass — for `stage06_assemble`, plus a `*_allmetrics.csv` with
  every candidate. The `epitope` arm is antibody-accessible (top-10 clash 1–2); the `hotspots` arm is
  **generation-limited** (top-10 clash 39–71; the whole 640-design pool clashes, floor 14) — shipped as a
  documented negative result on the minimal-hotspot-graft idea (manuscript `sec:dp4_8vdl`).

## LX minibinder arm (the PfEMP1/EPCR de-novo binders)

The same PfEMP1/EPCR project also produced a large set of **de-novo minibinders** (the "LX" designs) —
separate from the 8VDL *scaffolds* above, but folded into the shared DP4 deliverable.

- **Source:** `dp4_8vdl/data/LX_YYYYMMDD.csv` — ~484k generations with the LX filter flag `passes`;
  ~21.8k pass. 103-mer sequences, targets `fold_pfemp1_epcr_*`. **Gitignored** (~168 MB); an external
  artifact, provenance recorded here.
- **`scripts/08_add_minibinders.py`** — folds the *passing* LX designs into `data/libraries/dp4_library.csv`
  (`category=minibinder`), mapping sequence/uuid/target into the library shape, leaving the 12 episcaf
  metric/scoring columns blank (never scored on our axes), and carrying every native LX column as `lx_*`
  (plddt, pae, rmsd, iptm, hotspots, …) for post-hoc analysis. Idempotent; run after
  `scripts/stage06_assemble.py`. Takes the library from 14,241 episcaf → 36,000 rows.
- These are **not** episcaf-scaffolded (no episcaf metrics), but they ARE oligo-encoded **along with**
  the episcaf designs — one order file over the whole 36,000-row library (confirmed 2026-07-20; they
  carry no oligos of their own). This is the lab's PepSeq model: several projects combined into a single
  assay. Encode via `stage07_named_peptides.py` on the full library; `--exclude-category minibinder`
  would revert to episcaf-only.
