# Gemini bring-up: RFD4-Proteina + AlphaFold2 (ColabFold)

Phase 0 of the RFD4 engine swap (see `~/.claude/plans/jolly-questing-perlis.md`). This stands up the
**two compute tools** the new engine shells out to, and proves each one runs on Gemini with a smoke
test, **before** any orchestration is ported. Nothing here touches the live pilot.

Division of labor: Brandon runs everything on Gemini (cluster boundary). Claude wrote this runbook and
the two smoke jobs; Claude never ssh's to Gemini or submits these.

The engine (built next) only needs, at the end of this:
- a way to run `rfd4-generate ...` on a GPU node (a SIF **or** an activated venv), and
- a way to run `colabfold_batch ...` on a GPU node,
plus the paths to both, which go into `~/.episcaf-pilot/worker.env` and the new `handoff_rfd4.sbatch`.

Everything installs from **public** PyPI + the **public** HF weights (`nvidia/RFD4-Proteina`, ~1.7 GB) —
no HF token, no credentials, no antigen data leaves Gemini.

---

## 0. Confirm the GPU + CUDA version first

On a GPU node (`salloc`/`srun` onto `gpu-a100`, or inside the smoke sbatch):

```bash
nvidia-smi    # note "CUDA Version: 12.x" or "13.x" in the top-right
```

RFD4-Proteina ships `cu13` (default), `cu12`, and `cpu` builds. **Pick the extra that matches the
driver:** CUDA 13 → `all-cuda13` / `ACCELERATOR=cu13`; CUDA 12 → `all-cuda12` / `cu12`. Python must be
**3.12**. Fill the chosen value in once here and reuse it below:

```bash
ACCEL=cu13            # or cu12  — set to match nvidia-smi
```

---

## 1. Install RFD4-Proteina — pick ONE route

Both end at the same CLI (`rfd4-generate`). Route A (native venv) is usually simplest on a cluster
where you can't build containers; Route B (Apptainer) matches how we run gcloud and isolates deps.
**Try Route A first unless the lab standard is containers.**

### Route A — native venv with uv (simplest if Gemini has Python 3.12 + the CUDA driver)

```bash
cd "$HOME"
git clone --filter=blob:none https://github.com/RosettaCommons/RFD4-Proteina.git
cd RFD4-Proteina
# pip 26.2.2 first (Oct-2026 CA fix), then uv, then the package for your CUDA major:
pip install "pip @ https://files.pythonhosted.org/packages/3f/06/3e4420302363fbb71be3c7f35ce042f54bfca330e3e66b3a2e3c4cca579f/pip-26.2.2-py3-none-any.whl#sha256=fd7f4dfa8190e014013b9f08f3e973f4e4c26ffd179a204895533a7254813a70"
pip install "uv>=0.13.4"
uv venv --python 3.12 && source .venv/bin/activate
# CUDA 13: --torch-backend cu132 "rfd4-proteina[all-cuda13]"
# CUDA 12: --torch-backend cu126 "rfd4-proteina[all-cuda12]"
uv pip install --torch-backend cu132 "rfd4-proteina[all-cuda13]"
```

The weights cache on first run into `HF_HOME` (default `~/.cache/huggingface`). To pre-cache, run the
smoke in §3. The engine will activate this venv in `handoff_rfd4.sbatch`.

**Gemini specifics to confirm:** is Python 3.12 available (module load / conda)? Does a build node vs
a GPU node matter for the install? Pick a durable location for the checkout + venv on `$HOME` (NOT
`/scratch`, which is ephemeral).

### Route B — Apptainer/Singularity SIF (if you prefer a container)

```bash
cd "$HOME/RFD4-Proteina"      # same checkout as above
# Build the SIF from the shipped definition (needs fakeroot on the cluster):
apptainer build --fakeroot "$HOME/.episcaf-pilot/rfd4-proteina-$ACCEL.sif" \
    --build-arg ACCELERATOR=$ACCEL apptainer/rfd4-proteina.def
```

If `--fakeroot` isn't enabled for your account, ask the Gemini admins to enable it, or build on a host
where it is and copy the `.sif` over. Run form (used by the engine):

```bash
apptainer run --nv \
  -v "${HF_HOME:-$HOME/.cache/huggingface}:/opt/rfd4/cache/huggingface" \
  "$HOME/.episcaf-pilot/rfd4-proteina-$ACCEL.sif" rfd4-generate <args...>
```

---

## 2. Install ColabFold (AlphaFold2) — the oracle

John's composite is calibrated to **AF2**, so we use ColabFold (not AF3). The standard cluster install
is localcolabfold:

```bash
# On a durable $HOME path. The installer sets its target to $(pwd)/localcolabfold, so run it from
# $HOME (not from inside the clone, or you get a nested localcolabfold/localcolabfold).
cd "$HOME"
git clone https://github.com/YoshitakaMo/localcolabfold.git
bash "$HOME/localcolabfold/install_colabbatch_linux.sh"   # note: colabBATCH; downloads Miniforge + AF2 params
# adds colabfold_batch under $HOME/localcolabfold/colabfold-conda/bin
export PATH="$HOME/localcolabfold/colabfold-conda/bin:$PATH"
```

We run per design with: 5 models, 3 recycles, no relax, MSA on, no templates (matches John's AF2
settings). The MSA step uses the public MMseqs2 API by default; if Gemini blocks outbound traffic from
GPU nodes, we'll need a local MSA/database or the data-mover node for the MSA step — **flag if the
smoke in §4 can't reach the MSA server.**

---

## 3. Smoke test RFD4-Proteina (GPU)

Edit the two marked lines in `rfd4_smoke.sbatch` (partition/account and the Route A venv path **or**
Route B SIF path), then:

```bash
sbatch episcaf_service/gcloud/rfd4_smoke.sbatch
```

It runs the quickstart one-design generation (`inference_experiment=quickstart`, which bundles its own
example input, turns metrics off and skips the compile warm-up — no `inputs=` needed) — success means
RFD4-Proteina generates on Gemini. Expect a PDB/CIF under `runs/smoke/`. Paste the tail of the log if
it errors.

Note (Gemini QOS): an srun/sbatch of ≤2 h lands on the `axshort` QOS. A full RFD4 campaign may need a
longer wall time / different QOS — size `handoff_rfd4.sbatch` accordingly in Phase 2.

## 4. Smoke test ColabFold (GPU)

```bash
sbatch episcaf_service/gcloud/af2_smoke.sbatch
```

It folds one short sequence and writes a `*_rank_001_*.pdb` + the per-model scores (PAE). Success means
AF2 runs and we can read `Mean PAE`. If it hangs at the MSA step, that's the outbound-traffic caveat in §2.

---

## 5. Record the paths for the engine

Once both smoke tests pass, note these for the next phase (they go into `worker.env` / `handoff_rfd4.sbatch`):

| Setting | Route A (venv) | Route B (SIF) |
|---|---|---|
| `EPISCAF_RFD4_RUN` | `source <checkout>/.venv/bin/activate; rfd4-generate` | `apptainer run --nv -v …:… <sif> rfd4-generate` |
| `EPISCAF_COLABFOLD` | `<…>/colabfold-conda/bin/colabfold_batch` | same |
| `HF_HOME` | the weights cache dir | the weights cache dir (mounted) |
| CUDA `ACCEL` | `cu13`/`cu12` | `cu13`/`cu12` |

That's Phase 0. The engine (`episcaf_service/rfd4/`) is being built in parallel and will call these two
tools; nothing in the live pilot changes until Phase 2 wires it in.
