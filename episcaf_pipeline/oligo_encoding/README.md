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

## Encoder setup (one-time, on Gemini)
Per Erin (2026-07-13): **build the encoder from the GitHub master repo, using her conda env** — the
env carries everything needed to both compile and run (g++/OpenMP, pandas, and h2o 3.20.0.8 for the
model). So the tool is not vendored here; it is built once on the cluster.

```bash
# 1. create the env from Erin's yml (has g++/OpenMP/pandas/h2o):
conda env create -f /tgen_labs/Immunology/ekelley/DP3_PepSeq_Library_Design/new_codon_weights/pepseq_encoding.yml
conda activate pepseq_encoding      # confirm the name with `conda env list`
conda install -c conda-forge "openjdk=8" -y   # env ships no JRE; H2O (step 2) needs Java 8 -- see (d)

# 2. build the encoder from master (recipe verified 2026-07-13 -- see the three gotchas below):
git clone https://github.com/LadnerLab/Library-Design.git
cd Library-Design/oligo_encoding
sed -i 's/oligo_ecoding/oligo_encoding/g' Makefile      # (a) fix upstream Makefile typo
module load GCC/9.3.0                                   # (b) env has no g++; conda gxx won't solve
make optimized CC="g++ -static-libstdc++ -static-libgcc"  # (c) static C++ runtime -> ./main runs under the env
export TOOL_DIR=$(pwd)               # this dir holds: main, oligo_encoding.py, the model
```

Four gotchas hit on the first run, all handled above (the encode step is VERIFIED end-to-end on a
50-peptide smoke test, 2026-07-14):

- **(a) Makefile typo.** Master misspells its build target as `oligo_ecoding` (missing the `n`) in the
  `optimized`/`profile`/`debug` rules, so `make optimized` fails with `No rule to make target
  'oligo_ecoding'`. The `sed` fixes it. The compiled executable is named **`main`** (the real target
  runs `... -o main`), which happens to match ekelley's binary name.
- **(b) No compiler in the env.** `pepseq_encoding` ships no `g++` (only `/usr/bin/gcc` exists on the
  node), and `conda install gxx_linux-64` cannot solve against the env's pinned old libs
  (libstdcxx-ng 8.2.0 / libgomp 11.2.0). So load a cluster `GCC` module for the build instead.
- **(c) libstdc++ skew.** The sbatch jobs activate the pinned env, whose `libstdc++` is older than the
  build compiler's, so a normally-linked `main` would fail at runtime with a `GLIBCXX` error. Building
  with `-static-libstdc++ -static-libgcc` bakes the C++ runtime into `main`; the only shared lib it then
  needs is `libgomp` (OpenMP), which the env has and is backward-compatible. Verified: `./main` runs.
- **(d) No Java → step 2 dies.** The NN selector is H2O, which is a *Java* platform driven from Python,
  and the env ships no JRE (`H2OStartupError: Cannot find Java`). Install **Java 8** into the env --
  h2o 3.20.0.8 is a 2018 release supporting Java 7--10, so a modern JDK is the wrong choice:
  ```bash
  conda install -c conda-forge "openjdk=8" -y     # into pepseq_encoding
  ```
  Verified: H2O 3.20.0.8 starts on OpenJDK 1.8.0_472 and scores the candidates.

The selector script is `oligo_encoding.py` (master) rather than ekelley's `encoding_with_nn.py`. The
two sbatch scripts here **auto-detect both** binary (`main`/`oligo_encoding`) and selector
(`encoding_with_nn.py`/`oligo_encoding.py`); override with `BIN`/`SEL`. The DP3 model
`DeepLearning_model_R_1539970074840_1_20181019` is in the master repo. The sbatch scripts default
`CONDA_ENV=pepseq_encoding` and activate it themselves, so pointing `TOOL_DIR` at the built dir is all
that's needed.

## Running DP4
The only thing that changes from DP3 is the **`-i` input**: point it at the DP4 named-peptides file
(`library_member,sequence`, no header). Everything else — weights, model, GC, trials, adapters —
stays as the DP3 recipe above.

**The input is generated:** `scripts/stage07_named_peptides.py` slices `library_member,sequence` out
of the assembled library and validates the encoder's format (no header, unique names, standard
residues, every line ≤ 128) → **`data/libraries/dp4_named_peptides.csv`** (12,251 lines, all 103-mers,
line length 109–113). Regenerate with:

```bash
python scripts/stage07_named_peptides.py \
  --library data/libraries/dp4_library.csv --out data/libraries/dp4_named_peptides.csv
```

Two parameterized SLURM scripts here drive the run (better-documented rewrites of ekelley's
`run_pepseq_design_step{1,2}.sh`, which are kept as the pinned DP3 reference). Validate on the
50-peptide smoke-test input first (`data/libraries/dp4_named_peptides.test50.csv`, from
`stage07_named_peptides.py --sample 50`), then run the full file:

```bash
export TOOL_DIR=/path/to/Library-Design/oligo_encoding   # the built master dir (see setup above)

# smoke test (50 peptides), in a rundir holding the weights + test input:
cp episcaf_pipeline/oligo_encoding/codon_weights_updated.csv data/libraries/dp4_named_peptides.test50.csv <rundir>/
cd <rundir>
INPUT=dp4_named_peptides.test50.csv sbatch <repo>/episcaf_pipeline/oligo_encoding/encode_step1_generate.sbatch
sbatch <repo>/episcaf_pipeline/oligo_encoding/encode_step2_select.sbatch   # after step 1 finishes

# full run (12,251 peptides): same two commands with INPUT=dp4_named_peptides.csv
```

The scripts carry the DP3 parameters as defaults; override any by exporting it first (`INPUT=...`,
`TOOL_DIR=...`, `MODEL=...`, `GC=...`, `CONDA_ENV=...`). Step 1 writes `out_seqs` + `output_ratio`;
step 2 writes `DP4_best_encodings`. The whole thing runs on Gemini in the `pepseq_encoding` env.

## The order-file gap (open item)
Step 2 writes `*_best_encodings`. The DP3 **order file** that actually goes to synthesis
(`examples/DP3_order_file.sample.csv`, columns `Seq ID,nucleotide_encoding_with_twist_adapters`,
Twist adapters flanking each oligo) is a **further reformat + adapter step that is not in ekelley's
two scripts**. It is most likely either (a) `encoding_with_nn.py`'s `--adapter` default already
flanking the oligo so the order file is just a column rename of `*_best_encodings`, or (b) a small
unshared script. **Confirm with Erin** before treating the two-step run as the whole pipeline; once
known, it becomes a tiny `stage07c` (or a flag on step 2) and gets documented here.
