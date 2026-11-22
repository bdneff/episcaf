"""Gemini acceptance job for the checked-in hen egg lysozyme reference (3HFM chain Y).

The frozen selection identifies residue ranges on that reference. All model input
coordinates come from the repository. This module never submits a scheduler job.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import subprocess
import sys
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
import pandas as pd

from episcaf_analysis.presets import PRESETS
from episcaf_analysis.score import score
from episcaf_service.staging import candidate_plan

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "episcaf_v3/energetics/md/3hfm_apo/structures/3hfm.pdb"
SCRIPTS = ROOT / "scripts"


def reference_pdb(chain="Y"):
    atoms = [line[:21] + chain + line[22:] for line in REFERENCE.read_text().splitlines()
             if line.startswith("ATOM  ") and line[21] == "Y"]
    if not atoms:
        raise ValueError("Missing 3HFM lysozyme coordinates")
    return "\n".join(atoms) + "\nTER\nEND\n"


def validate_regions(regions):
    from episcaf_pipeline.build_growpatch_designs import parse_epitope
    for region in regions:
        islands = parse_epitope(region)
        if not islands or any(c != "Y" or not 1 <= a <= b <= 129 for c, a, b in islands):
            raise ValueError("Select regions on 3HFM lysozyme chain Y (residues 1–129).")


def prepare(directory, regions, plan=None):
    from episcaf_pipeline.build_growpatch_designs import parse_epitope
    validate_regions(regions)
    # The shipped fixed-PDB and stage05 extractors use chain A and positional
    # residue indices. 3HFM Y has contiguous native positions 1..129.
    if plan is None:
        plan = candidate_plan("lysozyme", "A", [r.replace("Y:", "A:") for r in regions])
    native = directory / "native"
    native.mkdir(mode=0o700)
    coordinates = reference_pdb("A")
    ledger = directory / "designs.csv"
    with ledger.open("x", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(plan["contigs"][0]))
        writer.writeheader()
        writer.writerows(plan["contigs"])
    inputs = directory / "inputs"
    inputs.mkdir()
    for row in plan["contigs"]:
        pdb = native / (row["id"] + ".pdb")
        pdb.write_text(coordinates)
        key = design_key(row)
        positions = {f"A{i}": "BKBN" for _, a, b in parse_epitope(row["growpatch_epitope"])
                     for i in range(a, b + 1)}
        spec = dict(input=str(pdb), length=row["contig_length"], contig=row["contig_string"].replace("/", ","),
                    select_fixed_atoms=positions, select_unfixed_sequence=False)
        (inputs / (key + ".json")).write_text(json.dumps({key: spec}, indent=2))
    return plan


def design_key(row):
    return f'{row["id"]}__contig{row["contig_id"]}__seed0__rep0'


def stage(directory, name):
    temp = directory / "stage.tmp"
    temp.write_text(name)
    temp.replace(directory / "stage.txt")


def collect(directory, expected):
    from episcaf_analysis.aomclash import score as accessibility
    from scripts.stage04_af3_emit_jsons import extract_sequence_from_pdb

    output = directory / "04_af3/outputs"
    analysis = directory / "05_analysis"
    analysis.mkdir()
    args = ["--af3_out_dir", str(output), "--mpnn_pdb_dir", str(directory / "03_mpnn"),
            "--ledger", str(directory / "designs.csv"), "--native_dir", str(directory / "native"),
            "--out", str(analysis / "metrics.parquet"), "--workers", "4"]
    subprocess.run([sys.executable, str(SCRIPTS / "stage05_extract_metrics.py"), *args], check=True)
    rows = pd.read_csv(analysis / "metrics.csv")
    ids = rows["predID"]
    if set(ids) != set(expected) or ids.duplicated().any() or not rows["status"].eq("ok").all():
        raise ValueError("Final metrics do not cover every submitted design.")
    # Keep every design and report the supplied no-antibody preset on its native
    # scale. Raw occupancy is an input to the percentile rank, never a pass gate.
    rows["antigen"] = "3HFM lysozyme"
    rows["design_id"] = ids
    ledger = pd.read_csv(directory / "designs.csv").set_index(["id", "contig_id"])
    targets = {}
    models = {}
    for row in rows.itertuples():
        contig = ledger.loc[(row.id, row.contig_id)]
        region = contig["growpatch_epitope"]
        if region not in targets:
            native = directory / "native" / (row.id + ".pdb")
            if native.read_text() != reference_pdb("A"):
                raise ValueError("Scoring requires the frozen 3HFM lysozyme reference.")
            targets[region] = accessibility.build_target(str(native), region, exclude_pdb="3hfm")
        if targets[region] is None:
            raise ValueError(f"Could not score lysozyme region {region}.")
        model, = (output / row.predID.lower()).glob("*_model.cif")
        positions = accessibility.design_epi_positions(contig["contig_string"], 1)
        rows.at[row.Index, "s_p1_free"] = accessibility.aomclash_score(targets[region], str(model), positions)
        rows.at[row.Index, "epitope"] = region.replace("A:", "Y:")
        models[row.predID] = model
    rows["source"] = "gemini"
    rows["sequence"] = [extract_sequence_from_pdb(expected[i]) for i in ids]
    preset = PRESETS["twelvemer"]
    metrics = dict(preset["metrics"])
    metrics["s_p1_free"] = metrics.pop("aomclash_pred")
    if not np.isfinite(rows[list(metrics)].to_numpy(dtype=float)).all():
        raise ValueError("One or more final designs lack a metric required for scoring.")
    ranked = score(rows, {**preset, "metrics": metrics, "select": {}}).drop(columns=["af3_dir"])
    ranked.to_csv(analysis / "ranked.csv", index=False)
    temp = directory / "designs.zip.tmp"
    with ZipFile(temp, "x", compression=ZIP_DEFLATED) as archive:
        archive.write(analysis / "ranked.csv", "ranked.csv")
        for design, model in models.items():
            archive.write(model, f"models/{design}.cif")
    temp.replace(directory / "designs.zip")


def execute(directory):
    from collections import Counter
    from scripts.stage03_mpnn_fixed_pdbs import cif_to_fixed_pdb, fixed_from_contig
    from scripts.stage04_af3_emit_jsons import extract_sequence_from_pdb
    from episcaf_pipeline.build_growpatch_designs import parse_epitope
    from scripts.stage03_mpnn_submit import MPNN_ENV, MPNN_SCRIPT, SEQS_PER_STRUCT
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run the lysozyme job inside a Gemini GPU allocation.")
    if "\n" in str(directory) or "'" in str(directory):
        raise ValueError("The AF3 wrapper requires a run directory without quotes or newlines.")
    plan = prepare(directory, json.loads((directory / "selection.json").read_text())["regions"])
    stage(directory, "Generating backbones")
    fixed = directory / "fixed"
    fixed.mkdir()
    motifs = {}
    for row in plan["contigs"]:
        key = design_key(row)
        out = directory / "02_rfd3" / key
        subprocess.run(["rfd3", "design", f"++inputs={directory / 'inputs' / (key + '.json')}",
                        f"++out_dir={out}", "n_batches=1", f"diffusion_batch_size={row['backbone_trials']}"],
                       cwd=directory, check=True)
        models = sorted(out.rglob("*.cif.gz"))
        if len(models) != row["backbone_trials"]:
            raise ValueError(f"Unexpected RFd3 model count for {key}: {len(models)}")
        positions = fixed_from_contig(row["contig_string"])
        native = extract_sequence_from_pdb(directory / "native" / (row["id"] + ".pdb"))
        motif = "".join(native[a-1:b] for _, a, b in parse_epitope(row["growpatch_epitope"]))
        for i, model in enumerate(models):
            path = fixed / f"{key}_model_{i}_fixed.pdb"
            if cif_to_fixed_pdb(model, positions, path) != "ok":
                raise ValueError(f"Unreadable backbone for {key}")
            motifs[path.stem] = (positions, motif)
    stage(directory, "Designing sequences")
    sequences = directory / "03_mpnn"
    sequences.mkdir()
    subprocess.run([str(Path(MPNN_ENV) / "bin/python"), MPNN_SCRIPT, "-pdbdir", str(fixed),
                    "-outpdbdir", str(sequences), "-seqs_per_struct", str(SEQS_PER_STRUCT),
                    "-temperature", "0.1", "-relax_cycles", "0"], check=True)
    expected = {}
    counts = Counter()
    for pdb in sequences.rglob("*_fixed_dldesign_*.pdb"):
        backbone = pdb.stem.rsplit("_dldesign_", 1)[0]
        if backbone not in motifs or pdb.stem in expected:
            raise ValueError("Unexpected ProteinMPNN design ID.")
        positions, motif = motifs[backbone]
        sequence = extract_sequence_from_pdb(pdb)
        if len(sequence) != 103 or "".join(sequence[i] for i in positions) != motif:
            raise ValueError(f"ProteinMPNN design {pdb.stem} does not preserve the selected lysozyme region.")
        expected[pdb.stem] = pdb
        counts[backbone] += 1
    if any(counts[key] != SEQS_PER_STRUCT for key in motifs):
        raise ValueError("ProteinMPNN output count does not match the frozen run.")
    stage(directory, "Predicting structures")
    inputs = directory / "04_af3/inputs"
    subprocess.run([sys.executable, str(SCRIPTS / "stage04_af3_emit_jsons.py"), "--mpnn_pdb_dir",
                    str(sequences), "--out_dir", str(inputs)], check=True)
    if {p.stem for p in inputs.glob("*.json")} != set(expected):
        raise ValueError("AF3 inputs do not match the designed sequences.")
    for start in range(0, len(expected), 100):  # existing stage04 wrapper packs 100 sequences per task
        subprocess.run(["bash", str(SCRIPTS / "stage04_af3_array.sbatch"), str(directory)],
                       env={**os.environ, "SLURM_ARRAY_TASK_ID": str(start // 100 + 1)}, check=True)
    stage(directory, "Scoring designs")
    collect(directory, expected)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run_directory", type=Path)
    execute(parser.parse_args().run_directory.resolve())
