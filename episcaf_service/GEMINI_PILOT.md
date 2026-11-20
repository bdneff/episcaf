# Advisor pilot on Gemini

John can upload a PDB or mmCIF, preview it, choose its antigen chain, select growPatch regions and submit a private request. The request follows him across sign-ins under **Results → My runs**. For the first automatic GPU run, the candidate page provides the checked-in hen egg lysozyme PDB from 3HFM. Select chain Y. Other structures enter **Awaiting lab review**. The browser shows which path applies before submission.

## Start the private service

Brandon runs these commands. Use the usual `ssh gemini` alias. In a fresh checkout on the submission host, keep the service and dispatcher on that same host and code revision. The GPU node must see the checkout and the run directory at the same absolute paths. The repository's earlier Gemini service used a shared home directory. Verify shared access, available quota, GPU partition, model environments and the site's policy for long-lived login-host services before the first live run. The login host may lack the compute nodes' `/tgen_labs` mount.

```bash
ssh gemini
umask 077
git clone --branch pilot/private-simulated-runs --single-branch https://github.com/bdneff/episcaf.git "$HOME/episcaf-advisor"
cd "$HOME/episcaf-advisor"
mkdir -p "$HOME/.episcaf-advisor" "$HOME/episcaf_runs/advisor-jobs"
chmod 700 "$HOME/.episcaf-advisor" "$HOME/episcaf_runs/advisor-jobs"
```

Save the configuration once on that host. The example uses the private home directory for SQLite. If home is network-backed, confirm reliable SQLite file locking with your administrators or choose a private persistent host-local volume for the database. Only the web service and dispatcher access SQLite. Store uploads, final models and logs outside both the checkout and the legacy workspace's `runs/` tree.

```bash
cat > "$HOME/.episcaf-advisor/session.env" <<'ENV'
umask 077
cd "$HOME/episcaf-advisor"
module load Mamba/25.3.0-3
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$HOME/rfd3/env/rfd3_py312"
export EPISCAF_PILOT_MODE=gemini-lysozyme
export EPISCAF_PILOT_DB="$HOME/.episcaf-advisor/runs.sqlite"
export EPISCAF_GEMINI_RUNS="$HOME/episcaf_runs/advisor-jobs"
ENV
printf 'export EPISCAF_GEMINI_HOST=%q\n' "$(hostname -f)" >> "$HOME/.episcaf-advisor/session.env"
source "$HOME/.episcaf-advisor/session.env"
python -m pip install -r episcaf_service/requirements.txt
python -c 'import fastapi, uvicorn, gemmi, pyarrow, freesasa, Bio, scipy'
command -v sbatch sacct squeue rfd3
python -m episcaf_service.pilot issue john > "$HOME/.episcaf-advisor/john.token"
```

Issue John's token once and give it to him privately. Issuing again for `john` rotates his token and preserves his runs. Revoke with `python -m episcaf_service.pilot revoke john`.

Use two persistent terminal sessions on the configured host, for example two tmux windows. In the first, run the web service.

```bash
source "$HOME/.episcaf-advisor/session.env"
python -m uvicorn episcaf_service.app:app --host 127.0.0.1 --port 8810
```

In the second, start the dispatcher when ready to accept lysozyme jobs under your Slurm account.

```bash
source "$HOME/.episcaf-advisor/session.env"
python -m episcaf_service.gemini
```

To exercise one dispatch or status refresh manually, run `python -m episcaf_service.gemini --once`. The continuous dispatcher checks for new requests every ten seconds and refreshes Slurm status about every thirty minutes. It stores the numeric job ID, retrieves the pending start estimate from the queue and records the actual completion timestamp from accounting. A lost receipt triggers a lookup by the unique Slurm job name. For a prolonged **submitting** state, inspect accounting and queue records for `ep-<run ID>` before taking action. A failed job's stage, exit status and job ID appear in My runs. The private run directory holds the detailed Slurm log.

## What to send John

If John can SSH to the same submission host with his own account, he can run the following on his laptop. His alias must resolve to the host recorded in `EPISCAF_GEMINI_HOST`. Otherwise arrange a lab-approved private HTTPS address for the web service. The public Cloud Run URL serves staging.

```bash
ssh -N -L 8810:127.0.0.1:8810 gemini
```

With that tunnel open, visit [Find candidate regions](http://127.0.0.1:8810/growpatch_tool.html). For the first GPU check, download the linked reference, upload it, choose Y, turn off DiscoTope ranking, generate candidates and select one. Sign in beneath the candidates. **Prepare design inputs** shows eight backbone attempts. **Start Gemini run** opens My runs. The dispatcher continues after the browser closes. Sign back in to download the scores CSV and the ZIP containing a predicted structure for each design.

## Compute and score interpretation

The worker uses the repository's 3HFM chain Y coordinates. It maps Y to A internally for the pipeline and reports selected residues in their original Y numbering. Each candidate gets eight backbone attempts across its contig variants, followed by eight ProteinMPNN sequences per backbone. Three candidates therefore request 24 backbones and 192 sequence predictions. The job checks per-backbone sequence counts and motif preservation. The AF3 batch count comes from the actual inputs. Final results require complete metrics for every submitted design.

The first allocation requests one A100, eight CPUs, 64 GB and twelve hours. Walltime and throughput are unverified for this installation. The compute node needs the existing RFd3 environment and weights, ProteinMPNN/PyRosetta installation, Singularity, AF3 container and licensed model files. The wrapper uses the repository's established AF3 paths under `/tgen_labs` and `/ref_genomes`. Verify quota in the private shared run directory for intermediates and final structures. Use one region for the first check.

The final table ranks all predictions within the run using the existing no-antibody percentile composite. The raw occupancy metric is `s_p1_free`. The viewer's four-filter pass needs an actual reference-antibody clash measurement. The reference job applies no bare-backbone cutoff. Neither a composite rank nor a native control measures experimental binding.

## Handle other uploaded structures

In the private service environment, list pending requests and export an individual submission for inspection.

```bash
python -m episcaf_service.reviewed list
python -m episcaf_service.reviewed export "$RUN_ID"
```

The private directory contains the original upload, its SHA-256 digest and the frozen request with owner, chain and selected regions. Review and disposition happen under the operator's account. Return a reason to the researcher with `python -m episcaf_service.reviewed decline "$RUN_ID" "Reason for the researcher"`. To deliver results from an independently completed lab job, supply its real numeric Slurm receipt, a CSV with `run_id,design_id,epitope,sequence,composite` and a directory containing one `<design_id>.cif` per row. The epitope must match a selected region in the frozen request. The operator is responsible for checking scientific provenance.

```bash
python -m episcaf_service.reviewed publish "$RUN_ID" --job-id "$JOB_ID" --ranked "$SCORES_CSV" --models "$MODEL_DIR"
```

## Offline checks

Run `python tests/test_gemini_pilot.py` in the service environment with Gemmi, PyArrow and the FastAPI test-client dependencies. It uses real lysozyme coordinates, mocked Slurm receipts and a synthetic confidence matrix for the native control. It exercises ownership, queue reconciliation, trial allocation, failure handling, metric extraction, scoring and downloads. The browser acceptance test in the service README covers closing and reopening the simulator. A GPU run under Brandon's account is still required to verify the installed model stack and access from John's machine.
