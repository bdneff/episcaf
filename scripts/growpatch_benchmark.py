#!/usr/bin/env python3
"""Head-to-head benchmark of growPatch enumerators on known conformational epitopes.

Question: given an antigen with no known antibody, which surface-enumeration strategy best recovers
the true epitope among its candidates? We compare, on the same epitopes with the same cleanup and the
same coverage metric:

  * growpatch   -- seed-and-grow along the residue-adjacency graph (read.enumerate_patches)
  * eyeview-R   -- antibody's-eye-view direction-sampled facing caps (eyeview.enumerate_eyeview),
                   one column per footprint radius R

Epitopes come from `reference/sabdb_epitopes_4A_sample200_chain1.csv` (SAbDab, single-chain). We keep
the CONFORMATIONAL ones -- gap-10 island count in the kept range (default 2..3): 2+ islands means
discontinuous in sequence (the case linear tiling cannot do), <=3 keeps them inside the enumerators'
<=3-island footprint design. The 1-island (continuous) epitopes are run separately as a control.

For each epitope the target is its true nonGapFilled 4 A contacts that are present in the structure
with a complete backbone. Both methods get that target as `protect` (shields it from short-island
trimming), exactly as John's best_patch benchmark protects the reference epitope -- so the comparison
is about enumeration geometry, not cleanup. Metric per epitope, per method:

  best_single_coverage = max over candidates of |candidate ∩ target| / |target|
  n_candidates         = how many candidates the method emitted (fewer is better at equal coverage)

Structures are the gitignored local SAbDab download (`--structures`, <pdb>.cif.gz). Writes a per-epitope
CSV and prints a summary (mean/median best-single coverage, recall at >=0.8, mean candidate count).

Usage:
  python scripts/growpatch_benchmark.py \
      --table results/growpatch/sample200_chain1.csv \
      --structures ~/Desktop/projects/epitope-universe/data/structures \
      --islands 2,3 --radii 12,14 --out results/growpatch/benchmark_coverage.csv
"""
from __future__ import annotations

import argparse
import csv
import gzip
import os
import sys
import tempfile
import warnings

import numpy as np

warnings.simplefilter("ignore")
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from episcaf_analysis.growpatch import read as R, eyeview as EV   # noqa: E402


def _gunzip_to_tmp(gz_path, tmpdir):
    out = os.path.join(tmpdir, os.path.basename(gz_path)[:-3])
    with gzip.open(gz_path, "rt") as fh, open(out, "w") as o:
        o.write(fh.read())
    return out


def best_coverage(patches, target):
    """(best_single_coverage, best patch's n_res, best patch's n_islands) for a method's output."""
    if not patches or not target:
        return 0.0, None, None
    best = max(patches, key=lambda p: len(set(p["keys"]) & target))
    return len(set(best["keys"]) & target) / len(target), best["n_res"], best["n_islands"]


def run_one(cif, chain, true_spec, radii, sasa_cut, n_dirs):
    """All methods on one epitope -> dict of per-method (coverage, n_candidates, best size/islands)."""
    ag = R.read_antigen(cif, chain, relsasa=True)
    target = ag.parse_epitope(true_spec) & set(ag.by_key)
    if not target:
        return None, 0
    out = {"n_target": len(target), "n_target_islands": len(ag_islands(ag, target))}

    gp, _ = R.enumerate_patches(cif, chain, antigen=ag, protect_epitope=target)
    cov, nres, nisl = best_coverage(gp, target)
    out.update(gp_cov=round(cov, 3), gp_ncand=len(gp), gp_best_nres=nres, gp_best_nisl=nisl)

    for r in radii:
        ev = EV.enumerate_eyeview(ag, sasa_cut=sasa_cut, n_dirs=n_dirs,
                                  footprint_radius=r, depth=r, protect=target)
        cov, nres, nisl = best_coverage(ev, target)
        out.update({f"ev{r}_cov": round(cov, 3), f"ev{r}_ncand": len(ev),
                    f"ev{r}_best_nres": nres, f"ev{r}_best_nisl": nisl})
    return out, len(target)


def ag_islands(ag, keys):
    from episcaf_analysis.growpatch import core
    return core.islands(keys)


def summarize(rows, radii, label):
    print(f"\n=== {label}  (n={len(rows)} epitopes) ===")
    methods = ["gp"] + [f"ev{r}" for r in radii]
    names = {"gp": "growpatch"}
    names.update({f"ev{r}": f"eyeview r={r}" for r in radii})
    print(f"{'method':<16} {'mean_cov':>9} {'median_cov':>11} {'recall>=.8':>11} "
          f"{'mean_ncand':>11} {'mean_bestsz':>12}")
    for m in methods:
        cov = np.array([r[f"{m}_cov"] for r in rows], dtype=float)
        nc = np.array([r[f"{m}_ncand"] for r in rows], dtype=float)
        bsz = np.array([r[f"{m}_best_nres"] for r in rows if r[f"{m}_best_nres"] is not None], dtype=float)
        print(f"{names[m]:<16} {cov.mean():>9.3f} {np.median(cov):>11.3f} "
              f"{(cov >= 0.8).mean():>11.3f} {nc.mean():>11.1f} "
              f"{(bsz.mean() if len(bsz) else float('nan')):>12.1f}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--table", default="results/growpatch/sample200_chain1.csv")
    ap.add_argument("--structures", required=True, help="dir of <pdb>.cif.gz")
    ap.add_argument("--islands", default="2,3", help="gap-10 island counts to keep (conformational)")
    ap.add_argument("--radii", default="12,14", help="eyeview footprint radii to test")
    ap.add_argument("--sasa-cut", type=float, default=0.2)
    ap.add_argument("--n-dirs", type=int, default=150)
    ap.add_argument("--control-islands", default="1", help="island counts run as a separate control")
    ap.add_argument("--out", default="results/growpatch/benchmark_coverage.csv")
    ap.add_argument("--limit", type=int, default=0, help="cap epitopes (debug)")
    a = ap.parse_args(argv)

    radii = [int(x) for x in a.radii.split(",")]
    keep = set(a.islands.split(","))
    ctrl = set(a.control_islands.split(","))
    table = list(csv.DictReader(open(a.table)))

    done, skipped = [], []
    n = 0
    for row in table:
        icount = row["epitope_island_count_GapFilled"]
        if icount not in keep | ctrl:
            continue
        pdb = row["pdbId"].strip()
        chain = row["epitope_chains_nonGapFilled"].strip().split(",")[0].split(";")[0]
        gz = os.path.join(os.path.expanduser(a.structures), f"{pdb.lower()}.cif.gz")
        if not os.path.exists(gz):
            skipped.append((pdb, "no structure"))
            continue
        with tempfile.TemporaryDirectory() as td:
            try:
                cif = _gunzip_to_tmp(gz, td)
                rec, ntgt = run_one(cif, chain, row["epitope_island_ranges_nonGapFilled"],
                                    radii, a.sasa_cut, a.n_dirs)
            except Exception as e:
                skipped.append((pdb, f"{type(e).__name__}: {e}"))
                continue
        if rec is None:
            skipped.append((pdb, "no target residues present"))
            continue
        rec.update(pdbId=pdb, chain=chain, island_count=icount,
                   true_epitope=row["epitope_island_ranges_nonGapFilled"], is_control=(icount in ctrl))
        done.append(rec)
        n += 1
        if n % 20 == 0:
            print(f"  ... {n} epitopes done")
        if a.limit and n >= a.limit:
            break

    if not done:
        sys.exit("no epitopes benchmarked")
    cols = ["pdbId", "chain", "island_count", "is_control", "n_target", "n_target_islands",
            "true_epitope", "gp_cov", "gp_ncand", "gp_best_nres", "gp_best_nisl"]
    for r in radii:
        cols += [f"ev{r}_cov", f"ev{r}_ncand", f"ev{r}_best_nres", f"ev{r}_best_nisl"]
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows({k: r.get(k) for k in cols} for r in done)

    conf = [r for r in done if not r["is_control"]]
    cont = [r for r in done if r["is_control"]]
    if conf:
        summarize(conf, radii, f"CONFORMATIONAL epitopes (islands {','.join(sorted(keep))})")
    if cont:
        summarize(cont, radii, f"CONTINUOUS control (islands {','.join(sorted(ctrl))})")
    print(f"\nwrote {a.out}  ({len(done)} epitopes; skipped {len(skipped)})")
    if skipped:
        for pdb, why in skipped[:12]:
            print(f"  skip {pdb}: {why}")


if __name__ == "__main__":
    main()
