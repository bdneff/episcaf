# Stage 07 — peptide → DNA oligo encoding

The final synthesis-prep step. Takes the assembled DP4 **peptide** library (the `06_library`
output, Section `sec:assembly` of the manuscript) and encodes each 103-mer into a DNA oligo, then
adds Twist adapters to produce the synthesis order file. External LadnerLab tooling
(`oligo_encoding`); our repo holds the config + recipe + the handoff, the run happens on Gemini.

## How it works
Codon degeneracy means each peptide has many possible DNA encodings. The tool:
1. **`main`** (compiled C++) samples many candidate encodings per peptide, weighted by the codon
   table and targeting a GC ratio.
2. **`encoding_with_nn.py`** scores the candidates with a pretrained H2O deep-learning model and
   keeps the best per peptide, then appends Twist adapters → order file.

## Files here
- `codon_weights_updated.csv` — the **updated codon-weight table** (same as DP3; reuse for DP4).
  Format, no header: `aa,codon,weight,index` (e.g. `K,AAA,1.596,0`).
- `run_pepseq_design_step1.sh` / `run_pepseq_design_step2.sh` — the **exact DP3 recipe** (ekelley),
  kept as the reference to mirror for DP4. Pinned DP3 parameters:
  - step 1 (`main`): `-t 10000` trials, `-g 0.55` GC target, `-n 300` candidates/peptide,
    `-c 2` cores, `-p codon_weights_updated.csv`, `-i <named_peptides.csv>`,
    outputs `-s out_seqs -r output_ratio`.
  - step 2 (`encoding_with_nn.py`): `-m DeepLearning_model_R_1539970074840_1_20181019`,
    `--subsample 300`, `--read_per_loop 10`, `-n 1` (one best encoding per peptide),
    `-o <best_encodings>`.
- `examples/DP3_named_peptides.sample.csv` — the **`-i` input format**: no header, `name,seq`
  (e.g. `DP3_0001,DWTQVALAN...`), max line length 128. (6000 peptides in the full DP3 file.)
- `examples/DP3_order_file.sample.csv` — the **final output format**: header
  `Seq ID,nucleotide_encoding_with_twist_adapters`; `Seq ID` is `<pep>_<enc##>`. Twist adapters
  flank each oligo (5′ `…CCTATACTTCCAAGGCGCA`, 3′ `GGTGACTCTCTGTCTTGGC…`).
- `examples/DP3_named_peptides.csv`, `examples/DP3_order_file.csv` — full DP3 input/output, kept
  locally as fixtures (gitignored via `*.csv`; not committed).

## What is NOT here yet (the encoder tool)
The encoder, NN selector, and model are ekelley's cluster install, referenced by the run scripts:
`/home/ekelley/bin/Library-Design/oligo_encoding/` (`main`, `encoding_with_nn.py`, and
`DeepLearning_model_R_1539970074840_1_20181019`). The current GitHub master
(`github.com/LadnerLab/Library-Design/tree/master/oligo_encoding`) has diverged from this install
(it ships `oligo_encoding.py`, not `encoding_with_nn.py`, and a `oligo_encoding` make target, not
`main`), so to reproduce DP3 exactly we mirror **ekelley's** version, not master. [decide: vendor a
pinned copy of that dir here, or reference the cluster install]

## Running DP4
The only thing that changes from DP3 is the **`-i` input**: point it at the DP4 named-peptides file
that `06_library` exports (`library_member,sequence`, no header). Everything else — weights, model,
GC, trials, adapters — stays as the DP3 recipe above.

Two parameterized SLURM scripts here drive the run (better-documented rewrites of ekelley's
`run_pepseq_design_step{1,2}.sh`, which are kept as the pinned DP3 reference):

```bash
# on Gemini, in a working dir holding codon_weights_updated.csv + DP4_named_peptides.csv:
sbatch episcaf_pipeline/oligo_encoding/encode_step1_generate.sbatch   # candidates -> out_seqs, output_ratio
sbatch episcaf_pipeline/oligo_encoding/encode_step2_select.sbatch     # NN-select  -> DP4_best_encodings
```

Both default `TOOL_DIR` to ekelley's install and carry the DP3 parameters as defaults; override any
by exporting it first (e.g. `INPUT=...`, `MODEL=...`, `GC=...`). The encoder runs on Gemini — the
compiled `main`, the `encoding_with_nn.py` selector, and the H2O model all live there; we only
supply the input file and reuse the DP3 config.

## The order-file gap (open item)
Step 2 writes `*_best_encodings`. The DP3 **order file** that actually goes to synthesis
(`examples/DP3_order_file.sample.csv`, columns `Seq ID,nucleotide_encoding_with_twist_adapters`,
Twist adapters flanking each oligo) is a **further reformat + adapter step that is not in ekelley's
two scripts**. It is most likely either (a) `encoding_with_nn.py`'s `--adapter` default already
flanking the oligo so the order file is just a column rename of `*_best_encodings`, or (b) a small
unshared script. **Confirm with Erin** before treating the two-step run as the whole pipeline; once
known, it becomes a tiny `stage07c` (or a flag on step 2) and gets documented here.
