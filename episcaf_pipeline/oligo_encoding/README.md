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
residues, every line ≤ 128) → **`data/libraries/dp4_named_peptides.csv`** (15,324 lines, all 103-mers,
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

# full run (15,324 peptides): same two commands with INPUT=dp4_named_peptides.csv
```

The scripts carry the DP3 parameters as defaults; override any by exporting it first (`INPUT=...`,
`TOOL_DIR=...`, `MODEL=...`, `GC=...`, `CONDA_ENV=...`). Step 1 writes `out_seqs` + `output_ratio`;
step 2 writes `DP4_best_encodings`. The whole thing runs on Gemini in the `pepseq_encoding` env.

## The order file, and the adapter length (RESOLVED: build with 19)

There is **no missing script**. Step 2 already flanks each oligo with the Twist adapters and emits it
as the column `Nucleotide Encoding w/ Adapters`, so the synthesis order file is a two-column slice of
`*_best_encodings` — `Seq ID` and that column, renamed to `nucleotide_encoding_with_twist_adapters`.
Emit it with `scripts/stage07_order_file.py`.

**Use the 19-mer adapters (347-nt oligos).** `--adapter` is a *defaulted* flag, and the two encoder
installs ship **different defaults**:

| install | 5′ adapter | 3′ adapter | oligo (103-mer) |
|---|---|---|---|
| **standard / GitHub master `oligo_encoding.py`** (**use this**) | `CCTATACTTCCAAGGCGCA` (19) | `GGTGACTCTCTGTCTTGGC` (19) | **347 nt** |
| ekelley's install, as it happened to run for DP3 | `ACCTATACTTCCAAGGCGCA` (20) | `GGTGACTCTCTGTCTTGGCT` (20) | 349 nt |

Erin confirmed (2026-07-14): *"It is usually 19. It would matter, so I'll try to figure out what the
deal was. So let's build with 19."* So **DP4 uses the 19-mers**, and `encode_step2_select.sbatch`
now passes `--adapter` **explicitly** (`ADAPTER=`, defaulting to the 19-mers) rather than inheriting a
tool default. The 20-mers the DP3 order file carried appear to be a DP3-specific anomaly Erin is still
tracking down, not the intended length — so we do **not** copy them forward.

**How this came up.** Neither ekelley's DP3 script nor our first sbatch passed `--adapter`, so both
silently took whatever their local copy defaulted to. `oligo_encoding.py` *is* `encoding_with_nn.py`
renamed upstream (`9d3ef58`, 2025-03-27) — same script, not a rewrite — and upstream has carried the
19-mers since at least 2022 (`e9c7de2`). So a fresh clone of master gives 19; ekelley's install, via a
local unpushed edit, gave 20 for DP3. The gap was invisible: no error, nothing in the log, just a
2-nt-longer oligo. Note Erin's *"it would matter"* — do not treat the length as free; the earlier
guess here that the extra bases were probably inert (they sit outside the primer footprint) was **not**
confirmed, and she believes length is load-bearing. Pin it; don't inherit it.

**Reference-file provenance.** The local `examples/DP3_order_file.csv` (the 20-mer file, used only to
build and self-test the checker) is byte-for-byte
`/tgen_labs/Immunology/ekelley/DP3_PepSeq_Library_Design/new_codon_weights/DP3_order_file.csv`
(md5 `bb537db52035593264b6ac2da3644d65`, verified 2026-07-14) — the updated-codon-table run DP4 reuses.

**What the adapters are.** The adapter is a primer *binding site*, not a primer. It is the constant
landing pad that lets one universal primer pair amplify the whole pool, whatever each member carries in
the middle; the primer is a separate ordered reagent, the adapter's **reverse complement**. Confirmed
against the 10x/GEM primers Heather Mead circulated (Slack, 2026-07-10):
`rc(GCCAAGACAGAGAGTCACC) = GGTGACTCTCTGTCTTGGC`, exactly the 19-mer 3′ adapter. (An earlier note here
speculated the 3′ adapter was *translated* as a `GDSLSW` linker — that was wrong; its function is
priming.)

### Emitting and checking the order file
`scripts/stage07_order_file.py` writes the order file and checks **every** row: adapters present and
correct, total length as expected, and the core translating back to exactly the peptide it claims to
encode. An oligo order is expensive and irreversible, so this runs before anything is sent.

```bash
python scripts/stage07_order_file.py \
  --best-encodings <rundir>/DP4_best_encodings \
  --peptides data/libraries/dp4_named_peptides.csv \
  --out data/libraries/dp4_order_file.csv
```

The checker defaults to the standard **19-mer** adapters (what DP4 ships). The DP3 example file carries
the 20-mers, so it is verified by passing them explicitly — a self-test that the checker both accepts
the real DP3 order (all 6,000, cores translating) and, on the 19-mer default, correctly *rejects* it:

```bash
# passes: DP3 file checked against its own 20-mers
python scripts/stage07_order_file.py --verify \
  --order-file episcaf_pipeline/oligo_encoding/examples/DP3_order_file.csv \
  --peptides episcaf_pipeline/oligo_encoding/examples/DP3_named_peptides.csv \
  --prefix ACCTATACTTCCAAGGCGCA --suffix GGTGACTCTCTGTCTTGGCT

# fails (as it should): DP3 file against the 19-mer default
python scripts/stage07_order_file.py --verify \
  --order-file episcaf_pipeline/oligo_encoding/examples/DP3_order_file.csv \
  --peptides episcaf_pipeline/oligo_encoding/examples/DP3_named_peptides.csv
```
