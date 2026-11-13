#!/usr/bin/env python3
"""AOMclash backbone pre-filter: drop clashy RFd3 backbones before MPNN + AF3.

Runs AFTER RFd3, BEFORE MPNN. Every RFd3 backbone is scored for antibody-clash propensity with
AOMclash (`s_p1_free`, via episcaf_analysis/aomclash/score.py), and only the least-clashy backbones
per epitope are passed to MPNN -- the rest never cost an MPNN or AF3 run. This is John's biggest
compute win: the cheap learned score prunes the expensive pipeline up front. (`aom_score.py` backbone
mode is the reference; John's suggested backbone cutoff is aom_pred <= 7.)

This is a cluster-side stage (RFd3 backbones live on $WS); build and unit-test locally, run there.

Inputs:
  --backbones DIR   RFd3 output PDBs, named <design_id>*.pdb (design_id = <id>__contig<c>__seed<s>__rep<r>)
  --contigs CSV     the designs table (build_growpatch_designs / stage01): id, contig_id,
                    contig_string, and the epitope column (growpatch_epitope or epitope)
  --antigen PDB     the native antigen (the uploaded structure; provides the epitope frame + ref)
  --map PKL         the occupancy map (default: the vendored aom_map.pkl)
Keep rule (one of):
  --keep-frac F     keep the best fraction F of backbones per group (default 0.5)
  --cutoff C        keep backbones with s_p1_free <= C
  --group COL       rank within this column (default 'id' = per epitope)

Output: a scores CSV (design_id, epitope, s_p1_free, kept) and a kept-manifest of design_ids for MPNN.
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


def design_epi_positions(contig_string, resid_offset=1):
    """Residue ids of the epitope in an RFd3 output backbone, from the contig (chain-general).

    Walks the slash-delimited contig: a free span `n-n` advances the output position by n; a motif
    segment `<chain><a>-<b>` occupies the next (b-a+1) output positions (that is where the fixed
    epitope lands in the scaffold). Returns the output residue ids (chain A, 1-based by default), in
    the island order the contig lists them -- the same order `score.build_target` builds `ref` in.
    """
    pos, out = 0, []
    for tok in str(contig_string).split("/"):
        tok = tok.strip()
        if not tok:
            continue
        if re.match(r"^\d+-\d+$", tok):                      # free span "n-n"
            pos += int(tok.split("-")[0])
        else:                                                # motif "<chain><a>-<b>" (or "<chain><a>")
            m = re.match(r"^[A-Za-z][A-Za-z0-9]*?(-?\d+)-(-?\d+)$", tok) or \
                re.match(r"^[A-Za-z][A-Za-z0-9]*?(-?\d+)$", tok)
            if not m:
                raise ValueError(f"unparseable contig token {tok!r} in {contig_string!r}")
            a = int(m.group(1))
            b = int(m.group(2)) if m.lastindex == 2 else a
            span = b - a + 1
            out.extend(range(pos + resid_offset, pos + resid_offset + span))
            pos += span
    return out


def keep_mask(scored, group_col, keep_frac, cutoff):
    """Mark kept=True on the least-clashy backbones per group (by s_p1_free; lower is better)."""
    from collections import defaultdict
    groups = defaultdict(list)
    for r in scored:
        if r["s_p1_free"] is not None:
            groups[r.get(group_col, "")].append(r)
    for r in scored:
        r["kept"] = False
    for rows in groups.values():
        rows.sort(key=lambda r: r["s_p1_free"])
        if cutoff is not None:
            for r in rows:
                r["kept"] = r["s_p1_free"] <= cutoff
        else:
            n = max(1, round(len(rows) * keep_frac))
            for r in rows[:n]:
                r["kept"] = True


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--backbones", required=True)
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--antigen", required=True)
    ap.add_argument("--map", default=DEFAULT_MAP)
    ap.add_argument("--epitope-col", default="", help="column with the epitope spec (auto: growpatch_epitope|epitope)")
    ap.add_argument("--group", default="id")
    ap.add_argument("--keep-frac", type=float, default=0.5)
    ap.add_argument("--cutoff", type=float, default=None)
    ap.add_argument("--exclude-pdb", default=None, help="drop this antigen's own SAbDab complexes from the map")
    ap.add_argument("--resid-offset", type=int, default=1, help="first RFd3 output residue id")
    ap.add_argument("--out", default="aomclash_prefilter.csv")
    ap.add_argument("--kept-out", default="aomclash_kept_designs.txt")
    a = ap.parse_args(argv)

    contigs = list(csv.DictReader(open(a.contigs)))
    epi_col = a.epitope_col or ("growpatch_epitope" if "growpatch_epitope" in contigs[0]
                                else ("epitope" if "epitope" in contigs[0] else None))
    if not epi_col:
        sys.exit("no epitope column in --contigs (expected growpatch_epitope or epitope)")
    by_id_cid = {(r["id"], str(r["contig_id"])): r for r in contigs}

    targets = {}   # epitope spec -> Target (built once)
    scored, missing = [], 0
    for pdb in sorted(glob.glob(os.path.join(a.backbones, "*.pdb"))):
        m = _DESIGN_ID.match(os.path.basename(pdb))
        if not m:
            continue
        row = by_id_cid.get((m.group("id"), m.group("cid")))
        if row is None:
            missing += 1
            continue
        epitope = row[epi_col]
        if epitope not in targets:
            targets[epitope] = S.build_target(a.antigen, epitope, a.map, exclude_pdb=a.exclude_pdb)
        tgt = targets[epitope]
        epi_pos = design_epi_positions(row["contig_string"], a.resid_offset)
        val = S.aomclash_score(tgt, pdb, epi_pos) if tgt is not None else None
        scored.append(dict(design_id=os.path.basename(pdb).split(".pdb")[0], id=m.group("id"),
                           epitope=epitope, s_p1_free=(round(val, 3) if val is not None else None)))

    keep_mask(scored, a.group, a.keep_frac, a.cutoff)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["design_id", "id", "epitope", "s_p1_free", "kept"])
        w.writeheader()
        w.writerows(scored)
    kept = [r["design_id"] for r in scored if r["kept"]]
    with open(a.kept_out, "w") as fh:
        fh.write("\n".join(kept) + ("\n" if kept else ""))
    n_sc = sum(r["s_p1_free"] is not None for r in scored)
    print(f"scored {n_sc}/{len(scored)} backbones ({missing} had no contig row); kept {len(kept)} "
          f"({'cutoff ' + str(a.cutoff) if a.cutoff is not None else f'best {a.keep_frac:.0%}/' + a.group})")
    print(f"wrote {a.out} and {a.kept_out}")


if __name__ == "__main__":
    main()
