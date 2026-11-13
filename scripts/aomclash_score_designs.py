#!/usr/bin/env python3
"""Compute aomclash_pred on a run's designs -- the end-of-pipeline accessibility metric.

The counterpart to the backbone pre-filter: there AOMclash prunes RFd3 backbones; here it scores the
final AF3 designs so `aomclash_pred` can be the composite's accessibility term (replacing the
cylinder) and the Results tab can rank by it. Same scorer (episcaf_analysis/aomclash/score.py), just
pointed at the AF3 output (CIF) instead of the backbones (PDB).

Cluster-side (AF3 designs live on $WS); build/test local, run there.

Inputs:
  --designs DIR   AF3 (or MPNN) design structures, named <design_id>*.{cif,pdb}
                  (design_id = <id>__contig<c>__seed<s>__rep<r>...)
  --contigs CSV   the designs table (build_growpatch_designs / stage01): id, contig_id,
                  contig_string, and the epitope column (growpatch_epitope or epitope)
  --antigen PDB   the native antigen (epitope frame + ref)
  --map PKL       occupancy map (default: vendored aom_map.pkl)
Output:
  --out CSV       design_id, id, epitope, aomclash_pred
  --merge CSV     optional: left-join aomclash_pred into this metrics CSV (match on design_id) and
                  write the merged table to --out instead
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
import warnings

warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from episcaf_analysis.aomclash import score as S   # noqa: E402

DEFAULT_MAP = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                           "episcaf_analysis", "aomclash", "aom_map.pkl")
_DESIGN_ID = re.compile(r"^(?P<id>.+?)__contig(?P<cid>\d+)__seed(?P<seed>\d+)__rep(?P<rep>\d+)")


def score_designs(designs_dir, contigs_csv, antigen, map_path, epitope_col="",
                  resid_offset=1, exclude_pdb=None):
    """-> {design_id: aomclash_pred} plus a list of per-design rows (design_id,id,epitope,aomclash_pred)."""
    contigs = list(csv.DictReader(open(contigs_csv)))
    epi_col = epitope_col or ("growpatch_epitope" if "growpatch_epitope" in contigs[0]
                              else ("epitope" if "epitope" in contigs[0] else None))
    if not epi_col:
        raise SystemExit("no epitope column in --contigs (expected growpatch_epitope or epitope)")
    by_id_cid = {(r["id"], str(r["contig_id"])): r for r in contigs}

    targets, rows, scores = {}, [], {}
    files = sorted(glob.glob(os.path.join(designs_dir, "*.cif"))
                   + glob.glob(os.path.join(designs_dir, "*.pdb")))
    for path in files:
        m = _DESIGN_ID.match(os.path.basename(path))
        if not m:
            continue
        row = by_id_cid.get((m.group("id"), m.group("cid")))
        if row is None:
            continue
        epitope = row[epi_col]
        if epitope not in targets:
            targets[epitope] = S.build_target(antigen, epitope, map_path, exclude_pdb=exclude_pdb)
        tgt = targets[epitope]
        epi_pos = S.design_epi_positions(row["contig_string"], resid_offset)
        val = S.aomclash_score(tgt, path, epi_pos) if tgt is not None else None
        did = re.sub(r"\.(cif|pdb)$", "", os.path.basename(path))
        scores[did] = round(val, 3) if val is not None else ""
        rows.append(dict(design_id=did, id=m.group("id"), epitope=epitope,
                         aomclash_pred=scores[did]))
    return scores, rows


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--designs", required=True)
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--antigen", required=True)
    ap.add_argument("--map", default=DEFAULT_MAP)
    ap.add_argument("--epitope-col", default="")
    ap.add_argument("--resid-offset", type=int, default=1)
    ap.add_argument("--exclude-pdb", default=None)
    ap.add_argument("--merge", default=None, help="metrics CSV to join aomclash_pred into (by design_id)")
    ap.add_argument("--out", default="aomclash_pred.csv")
    a = ap.parse_args(argv)

    scores, rows = score_designs(a.designs, a.contigs, a.antigen, a.map, a.epitope_col,
                                 a.resid_offset, a.exclude_pdb)
    n_ok = sum(r["aomclash_pred"] != "" for r in rows)

    if a.merge:
        src = list(csv.DictReader(open(a.merge)))
        if not src:
            raise SystemExit(f"--merge {a.merge} is empty")
        idcol = "design_id" if "design_id" in src[0] else ("id" if "id" in src[0] else None)
        if not idcol:
            raise SystemExit("--merge CSV has no design_id/id column to join on")
        cols = list(src[0].keys()) + (["aomclash_pred"] if "aomclash_pred" not in src[0] else [])
        for r in src:
            r["aomclash_pred"] = scores.get(r[idcol], r.get("aomclash_pred", ""))
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols)
            w.writeheader()
            w.writerows(src)
        print(f"scored {n_ok}/{len(rows)} designs; merged aomclash_pred into {len(src)} rows -> {a.out}")
    else:
        with open(a.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["design_id", "id", "epitope", "aomclash_pred"])
            w.writeheader()
            w.writerows(rows)
        print(f"scored {n_ok}/{len(rows)} designs -> {a.out}")


if __name__ == "__main__":
    main()
