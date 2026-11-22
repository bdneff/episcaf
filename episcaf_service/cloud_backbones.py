"""Bounded Google Batch job for the two 3HFM lysozyme reference routes."""
import argparse
import csv
from functools import lru_cache
import json
import math
from pathlib import Path
import subprocess
import tempfile
import time
from zipfile import ZipFile, ZIP_DEFLATED

from episcaf_analysis.growpatch.read import detect_complex
from episcaf_pipeline.build_growpatch_designs import parse_epitope
from episcaf_service import lysozyme_job as lysozyme
from scripts.stage03_mpnn_fixed_pdbs import cif_to_fixed_pdb, fixed_from_contig
from scripts.aomclash_prefilter import keep_mask

ATTEMPTS = 16
LIMIT = 8


@lru_cache(maxsize=1)
def contact_region():
    info = detect_complex(str(lysozyme.REFERENCE), "Y")
    if info["partner_chains"] != ["H", "L"] or not info["epitope"]:
        raise ValueError("The 3HFM reference must contain the HyHEL-10 heavy and light chains.")
    return info["epitope"]


def region_for(route, region):
    if route == "antibody":
        if region != contact_region():
            raise ValueError("Use the contact epitope from chains H and L of the 3HFM reference.")
    elif route != "growpatch":
        raise ValueError("Unknown epitope route.")
    lysozyme.validate_regions([region])
    return region


def generate(directory, route, region):
    from episcaf_analysis.aomclash import score as raw
    from scripts.stage04_af3_emit_jsons import extract_sequence_from_pdb

    plan = lysozyme.prepare(directory, [region_for(route, region)])
    rows = plan["contigs"]
    scores = []
    models = directory / "models"
    models.mkdir()
    native = directory / "native" / (rows[0]["id"] + ".pdb")
    target = raw.build_target(str(native), region.replace("Y:", "A:"), exclude_pdb="3hfm")
    if target is None:
        raise ValueError("The selected lysozyme region could not be scored.")
    sequence = extract_sequence_from_pdb(native)
    motif = "".join(sequence[a-1:b] for _, a, b in parse_epitope(region))
    for index, row in enumerate(rows):
        key = lysozyme.design_key(row)
        attempts = ATTEMPTS // len(rows) + int(index < ATTEMPTS % len(rows))
        out = directory / key
        subprocess.run(["rfd3", "design", f"++inputs={directory / 'inputs' / (key+'.json')}",
                        f"++out_dir={out}", "n_batches=1", f"diffusion_batch_size={attempts}"],
                       cwd=directory, check=True)
        structures = sorted(out.rglob("*.cif.gz"))
        if len(structures) != attempts:
            raise ValueError("RFdiffusion3 output count did not match the fixed trial budget.")
        fixed = fixed_from_contig(row["contig_string"])
        for i, structure in enumerate(structures):
            name = f"{key}_model_{i}"
            pdb = models / (name + ".pdb")
            value = None
            if cif_to_fixed_pdb(structure, fixed, pdb) == "ok":
                seq = extract_sequence_from_pdb(pdb)
                if len(seq) == 103 and "".join(seq[j] for j in fixed) == motif:
                    value = raw.aomclash_score(target, str(pdb), [j+1 for j in fixed])
            scores.append(dict(source="cloud-backbone", design_id=name, epitope=region, route=route,
                               s_p1_free=value if value is not None and math.isfinite(value) else None))
    keep_mask(scores, "epitope", keep_frac=0.5, cutoff=None)
    kept = sorted((r for r in scores if r["kept"]), key=lambda r:r["s_p1_free"])
    if len(scores) != ATTEMPTS or len(kept) > LIMIT:
        raise ValueError("The run exceeded its backbone budget.")
    def write_csv(path, records):
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(scores[0]))
            w.writeheader()
            w.writerows(records)
    write_csv(directory / "backbone_scores.csv", scores)
    write_csv(directory / "retained.csv", kept)
    with ZipFile(directory / "designs.zip", "x", compression=ZIP_DEFLATED) as archive:
        for table in ("retained.csv", "backbone_scores.csv"):
            archive.write(directory / table, table)
        for r in kept:
            name = r["design_id"] + ".pdb"
            archive.write(models / name, "backbones/" + name)
    return len(kept)


def main():
    from google.cloud import storage
    import re
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--route", required=True, choices=("antibody", "growpatch"))
    parser.add_argument("--region", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--bucket", required=True)
    args = parser.parse_args()
    if not re.fullmatch("[0-9a-f]{32}", args.run):
        parser.error("Invalid run ID.")
    bucket = storage.Client().bucket(args.bucket)
    with tempfile.TemporaryDirectory() as tmp:
        root = Path(tmp)
        count = generate(root, args.route, args.region)
        for name in ("retained.csv", "designs.zip"):
            bucket.blob(f"{args.run}/{name}").upload_from_filename(str(root/name), if_generation_match=0)
        bucket.blob(f"{args.run}/done.json").upload_from_string(
            json.dumps({"attempted": ATTEMPTS, "retained": count, "completed": time.time()}), if_generation_match=0)


if __name__ == "__main__":
    main()
