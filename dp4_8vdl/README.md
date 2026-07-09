# DP4 — 8VDL PfEMP1 conserved-epitope scaffolding

Self-contained arm that scaffolds the conserved EPCR-binding epitope of **PfEMP1 CIDRα1.4**
(structure **8VDL**) and folds the top designs into the DP4 PepSeq library. It is kept as its own
subdirectory with its own Gemini job because it belongs to the sibling *minibinder / PfEMP1* project
(`~/Desktop/projects/minibinder`); it rides into DP4 only for synthesis, and is **consolidated from
the repo root at the end** (score → top 10 per run → 8-column annotated format → `stage06_assemble`).

## The epitope region — and why this 20-mer

```
Target : 8VDL  (HB3VAR03 CIDRα1.4 + C7 Fab; Reyes et al., Nature 636:182–189, 2024)
Chain  : C   (the CIDR antigen; chains H/L are the C7 antibody)
Epitope: residues 651–670  =  FDSFFFQVIYKFNEGEAKWN   (20 aa, one contiguous island)
Hotspots (paper): F655, F656, E666
```

PfEMP1 is the hypervariable *P. falciparum* surface antigen that drives severe malaria and escapes
immunity by antigenic variation. But the residues its CIDRα1 domain uses to grip the host receptor
**EPCR cannot vary freely** — so a binder or antibody aimed at that conserved contact surface could
recognize **many PfEMP1 variants at once**. In 8VDL that conserved, EPCR-binding contact surface is
the single contiguous stretch **chain C 651–670**; the 2024 paper singles out **F655/F656/E666** as
the functional hotspots within it. We take this **20-residue contiguous chain as the epitope region**
because it is the conserved, functionally constrained EPCR-binding surface the whole broadly-acting
strategy is built on — and because it is contiguous and fully resolved in the crystal (all 20 Cα
present), so it scaffolds as one clean island with no gap-bridging guesswork. A Rosetta predecessor
(6SNY) scaffolded this epitope and recapitulated EPCR-binding geometry but elicited weaker
adhesion-inhibitory antibodies than the native protein; RFD3 is not constrained to a particular fold,
so it explores a broader scaffold space for presenting it.

Because chains H/L are the C7 antibody, this is a **known-antibody run** — the real AF3 clash filter
(`af3_n_clash_res`) applies, exactly as for DP4 components C1/C2.

## Three runs: full epitope vs. minimal hotspots vs. contact footprint

We generate three independent contig sets and score each on its own, taking the **top 10** from each:

| run | fixed motif | what it asks |
|---|---|---|
| **`epitope20`** | the whole island **C651–670** (20 aa, 1 island) | present the *entire* conserved epitope in native geometry — strongest constraint, the direct analog of C1 |
| **`hotspots`** | only **F655, F656, E666** (2 islands) | present just the *functionally critical* residues, letting RFD3 design everything else — minimal constraint / hotspot graft |
| **`contact`** | the **4 Å contact epitope**: 13 residues (652,653,655–657,659–661,666,667,669,670,673) in **6 islands** | present exactly the residues that touch the C7 Fab (the AbDb/IEDB standard definition) — the true antibody footprint, most faithful but the hardest scaffold |

All are scaffolded into a constant **103-mer** (the PepSeq maximum); the fixed atoms carry the residues'
native crystal coordinates, so each run preserves the motif's 3-D arrangement. The `contact` islands are
natively close-packed, so it is generated with **`--native-gaps`** (hold the native inter-island spacing
1/1/4/1/2, randomize only the flanks) rather than the large random gaps the others allow — otherwise RFD3
would have to bridge natively-adjacent islands with strained loops. Comparing the three in the assay asks
whether the minimal hotspots suffice, the full window is needed, or the exact contact footprint is best.

## Pipeline (self-contained; runs on Gemini)

Same RFD3→ProteinMPNN→AlphaFold3 shape as the rest of episcaf, but chain-C / real-residue aware
(adapted from `minibinder/pfemp1/run_pfemp1_scaffold`). All commands from `dp4_8vdl/`.

```bash
# 1. contigs — one CSV per run (--n-contigs is John's "n designs" knob; 10 each → 10×8×8=640 designs)
python scripts/01_generate_contigs.py --target epitope20 --n-contigs 10 --out 01_contigs/epitope20.csv
python scripts/01_generate_contigs.py --target hotspots  --n-contigs 10 --out 01_contigs/hotspots.csv
python scripts/01_generate_contigs.py --target contact --native-gaps --n-contigs 10 --out 01_contigs/contact.csv

# 2. RFD3 input JSONs
python scripts/02_emit_rfd3_inputs.py --contigs_csv 01_contigs/epitope20.csv --out_dir 02_rfd3/epitope20/inputs
python scripts/02_emit_rfd3_inputs.py --contigs_csv 01_contigs/hotspots.csv  --out_dir 02_rfd3/hotspots/inputs

# 3. RFD3 (Gemini GPU) — one task per contig, 8 backbones each
sbatch --array=1-10%200 scripts/03_rfd3_array.sbatch 02_rfd3/epitope20
sbatch --array=1-10%200 scripts/03_rfd3_array.sbatch 02_rfd3/hotspots

# 4. FIXED backbone PDBs — the ONE 8VDL-specific step (multi-island FIXED positions), then
#    everything downstream reuses episcaf's proven stage03/stage04 (run per-run, from repo root):
python dp4_8vdl/scripts/04_make_fixed_pdbs.py \
    --contigs_csv dp4_8vdl/01_contigs/epitope20.csv \
    --rfd3_outputs_dir dp4_8vdl/02_rfd3/epitope20/outputs \
    --out_dir dp4_8vdl/runs/epitope20/03_mpnn/fixed_pdbs

# 5. MPNN (episcaf, proven) — 8 seqs/backbone
python scripts/stage03_mpnn_submit.py --fixed_pdb_dir dp4_8vdl/runs/epitope20/03_mpnn/fixed_pdbs \
    --outdir dp4_8vdl/runs/epitope20/03_mpnn/mpnn_pdbs --batch_size 300 --tag 8vdl_ep20

# 6. AF3 (episcaf, proven)
python scripts/stage04_af3_emit_jsons.py --mpnn_pdb_dir dp4_8vdl/runs/epitope20/03_mpnn/mpnn_pdbs \
    --out_dir dp4_8vdl/runs/epitope20/04_af3/inputs --seed 1
sbatch --array=1-N scripts/stage04_af3_array.sbatch dp4_8vdl/runs/epitope20   # N = ceil(#jsons/100)
#   (repeat 4–6 for hotspots)

# 7. consolidate — score (real H/L clash), take top 10 per run, emit 8-column for DP4  [TODO: 07]
python dp4_8vdl/scripts/07_consolidate.py --runs epitope20,hotspots --topk 10 \
    --out results/dp4_8vdl_top10.csv
```

`data/8VDL.pdb` is the crystal (chains C/H/L), committed here so the arm is self-contained. The final
`results/dp4_8vdl_top10.csv` is what `stage06_assemble` folds into the library (as its own target).

## Notes
- The RFD3/AF3 `sbatch` are adapted from the minibinder project; the stale hardcoded `--chdir` has
  been removed (they inherit the submit dir), but dry-run once on the current cluster before the full array.
- `--n-contigs` × 8 RFD3 × 8 MPNN designs per run; top 10 come out, so 10 contigs is ample. Note
  `epitope20` (one island, 2 scaffold segments) caps at 64 unique contigs; `hotspots` (3 segments)
  has no practical cap. RFD3 still makes 8 stochastic backbones per contig, so 10 contigs = 640 designs.
- MPNN/AF3 are **not ported** — the only 8VDL-specific step is `04_make_fixed_pdbs.py` (computes the
  multi-island FIXED positions and reuses episcaf's `cif_to_fixed_pdb`); steps 5–6 are episcaf's
  just-run `stage03_mpnn_submit.py` / `stage04_*`. The prior art's `04/05/06` are superseded and not
  carried in.
- **Still TODO:** `07_consolidate.py` — score the AF3 outputs (antibody preset, real H/L clash),
  take top 10 per run, emit the 8-column rows for `stage06_assemble`. I'll write it against the real
  AF3 outputs once they land.
