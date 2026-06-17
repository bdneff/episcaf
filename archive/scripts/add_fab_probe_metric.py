#!/usr/bin/env python3
"""
add_fab_probe_metric.py

A faithful clash proxy: instead of a geometric solid, dock a panel of REAL
generic antibodies over the epitope and count how many scaffold residues they'd
hit. For each design and each probe, the probe's antibody is placed in its own
native binding pose (real standoff + approach direction from the probe complex),
re-oriented onto this epitope's outward normal, spun about the approach axis, and
clashed against the scaffold. The metric is the MINIMUM clashing-residue count
over probes and spins -- "can any real antibody shape fit here."

This reuses the same clash definition as the ground-truth filter (scaffold heavy
atoms within --clash_cutoff of antibody heavy atoms), so units match
af3_n_clash_res. A design's OWN target complex is excluded as a probe (no leakage).

Adds columns:
    fab_min_clash   : best (least-clashing) placement over all probes/spins
    fab_mean_clash  : mean over probes (each at its best spin)
    fab_n_probes    : how many probes were usable for this design

Example
-------
    python scripts/add_fab_probe_metric.py \
        --metrics_csv runs/run_rfd3_mpnn/04_filter/metrics_cone.csv \
        --dp2_parquet datasets/dp2.parquet \
        --probe_dir   /tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned \
        --out_csv     runs/run_rfd3_mpnn/04_filter/metrics_fab.csv \
        --limit 1000
"""

import argparse
import gzip
import math
import random
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import gemmi
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

RADIUS_IFACE = 30.0      # keep probe antibody atoms within this of probe epitope
EPI_CUTOFF   = 5.0       # antigen residue is "epitope" if within this of antibody

TOKEN_COLS   = ["token", "assay_scaffolded_epitope_id"]
AF3_DIR_COLS = ["af3_dir"]
ID_COLS      = ["id"]
DP2_EPI_COLS = ["assay_scaffolded_epitope_chunk_resindices"]
DP2_EPI_MPNN = ["scaffolded_epitope_chunk_resindices"]
WS_COLS      = ["af3_window_start"]
TRUE_CLASH   = ["af3_n_clash_res"]
AB_CHAINS    = ("B", "C")


# --------------------------------------------------------------------------- IO
def parse_index_list(x) -> List[int]:
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return []
    if isinstance(x, (list, tuple, np.ndarray)):
        return [int(i) for i in x]
    s = str(x).strip().replace("[", "").replace("]", "").replace(",", " ")
    out = []
    for t in s.split():
        try:
            out.append(int(t))
        except ValueError:
            pass
    return out


def first_present(cols, cands):
    for c in cands:
        if c in cols:
            return c
    return None


def read_gemmi(p: Path) -> gemmi.Structure:
    s = str(p)
    if s.endswith(".gz"):
        with gzip.open(p, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
        return gemmi.make_structure_from_block(doc.sole_block())
    if p.suffix.lower() == ".cif":
        doc = gemmi.cif.read(str(p))
        return gemmi.make_structure_from_block(doc.sole_block())
    return gemmi.read_structure(str(p))


def get_chainA(st):
    for ch in st[0]:
        if ch.name == "A":
            return ch
    return st[0][0]


def find_af3_cif(af3_dir: Path) -> Optional[Path]:
    for pat in ("*_model.cif", "*_model.cif.gz", "model.cif", "model.cif.gz"):
        h = next(af3_dir.glob(pat), None) or next(af3_dir.rglob(pat), None)
        if h:
            return h
    return None


def load_af3_heavy(af3_dir: Path):
    """chain-A heavy-atom coords, per-atom residue index, and a CA mask."""
    cif = find_af3_cif(Path(af3_dir))
    if cif is None:
        return None
    try:
        st = read_gemmi(cif)
    except Exception:
        return None
    chA = get_chainA(st)
    coords, ridx, is_ca = [], [], []
    for i, res in enumerate(chA):
        for a in res:
            if a.element == gemmi.Element("H"):
                continue
            coords.append([a.pos.x, a.pos.y, a.pos.z])
            ridx.append(i)
            is_ca.append(a.name == "CA")
    if not coords:
        return None
    return (np.asarray(coords, float), np.asarray(ridx, int),
            np.asarray(is_ca, bool))


def load_probe(pdb: Path):
    """Returns (ab_centered Nx3, approach_dir 3) for a complex, or None.
    ab atoms are the antibody heavy atoms near the interface, centered on the
    probe's epitope centroid; approach_dir points epitope -> antibody COM."""
    try:
        st = read_gemmi(pdb)
    except Exception:
        return None
    model = st[0]
    ag, ab = [], []
    for ch in model:
        sink = ab if ch.name in AB_CHAINS else (ag if ch.name == "A" else None)
        if sink is None:
            continue
        for res in ch:
            for a in res:
                if a.element == gemmi.Element("H"):
                    continue
                sink.append([a.pos.x, a.pos.y, a.pos.z])
    if len(ag) < 5 or len(ab) < 20:
        return None
    ag = np.asarray(ag, float); ab = np.asarray(ab, float)

    # epitope = antigen atoms within EPI_CUTOFF of any antibody atom
    tree = cKDTree(ab)
    d, _ = tree.query(ag, k=1)
    epi = ag[d <= EPI_CUTOFF]
    if len(epi) < 3:
        return None
    epi_centroid = epi.mean(axis=0)
    approach = ab.mean(axis=0) - epi_centroid
    nrm = np.linalg.norm(approach)
    if nrm < 1e-6:
        return None
    approach /= nrm

    # keep antibody atoms near the interface (the bulk that sits over the epitope)
    near = np.linalg.norm(ab - epi_centroid, axis=1) <= RADIUS_IFACE
    ab_if = ab[near]
    if len(ab_if) < 20:
        return None
    return ab_if - epi_centroid, approach


# --------------------------------------------------------------------- geometry
def align_R(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Rotation R with R @ a == b for unit vectors a, b (Rodrigues)."""
    a = a / np.linalg.norm(a); b = b / np.linalg.norm(b)
    v = np.cross(a, b); c = float(np.dot(a, b))
    if c > 1 - 1e-8:
        return np.eye(3)
    if c < -1 + 1e-8:                       # antiparallel: 180 deg about any perp
        perp = np.array([1.0, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1.0, 0])
        axis = np.cross(a, perp); axis /= np.linalg.norm(axis)
        return rot_about(axis, math.pi)
    vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + vx + vx @ vx * (1.0 / (1.0 + c))


def rot_about(axis: np.ndarray, theta: float) -> np.ndarray:
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis; ct, st = math.cos(theta), math.sin(theta)
    return np.array([
        [ct + x*x*(1-ct),   x*y*(1-ct)-z*st, x*z*(1-ct)+y*st],
        [y*x*(1-ct)+z*st,   ct + y*y*(1-ct), y*z*(1-ct)-x*st],
        [z*x*(1-ct)-y*st,   z*y*(1-ct)+x*st, ct + z*z*(1-ct)],
    ])


def fab_clash(scaf_tree, scaf_res, scaf_n_res, ab_centered, approach,
              centroid, normal, spins, cutoff) -> int:
    """Min clashing-residue count for one probe over spins."""
    R = align_R(approach, normal)
    base = ab_centered @ R.T                     # approach now aligned to normal
    best = scaf_n_res + 1
    for phi in spins:
        placed = base @ rot_about(normal, phi).T + centroid
        hit = scaf_tree.query_ball_point(placed, cutoff)
        res = set()
        for idxs in hit:
            for j in idxs:
                res.add(int(scaf_res[j]))
        best = min(best, len(res))
        if best == 0:
            break
    return best


# --------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metrics_csv", required=True)
    ap.add_argument("--dp2_parquet", required=True)
    ap.add_argument("--probe_dir", required=True,
                    help="dir of antibody-antigen complex PDBs to use as probes")
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--max_probes", type=int, default=8)
    ap.add_argument("--n_spin", type=int, default=12)
    ap.add_argument("--clash_cutoff", type=float, default=4.0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    spins = [2 * math.pi * k / args.n_spin for k in range(args.n_spin)]

    # ---- build probe panel ----
    probe_files = sorted(Path(args.probe_dir).glob("*.pdb"))
    random.Random(args.seed).shuffle(probe_files)
    probes: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
    for pf in probe_files:
        if len(probes) >= args.max_probes:
            break
        pr = load_probe(pf)
        if pr is not None:
            probes[pf.stem] = pr           # stem ~ complex id (e.g. 7ox3_0P)
    if not probes:
        raise SystemExit("no usable probes found in --probe_dir")
    print(f"Probe panel ({len(probes)}): {list(probes)}")

    df = pd.read_csv(args.metrics_csv, low_memory=False)
    if args.limit > 0:
        df = df.head(args.limit).copy()
    print(f"  {len(df):,} rows")

    tok_col = first_present(df.columns, TOKEN_COLS)
    dir_col = first_present(df.columns, AF3_DIR_COLS)
    id_col  = first_present(df.columns, ID_COLS)
    ws_col  = first_present(df.columns, WS_COLS)

    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"] = (dp2["assay_scaffolded_epitope_id"]
                                          .astype(str).str.lower())
    epi_col  = first_present(dp2.columns, DP2_EPI_COLS)
    epi_mpnn = first_present(dp2.columns, DP2_EPI_MPNN)
    dp2_lookup = dp2.drop_duplicates("assay_scaffolded_epitope_id").set_index(
        "assay_scaffolded_epitope_id")

    fab_min = np.full(len(df), np.nan)
    fab_mean = np.full(len(df), np.nan)
    fab_np = np.full(len(df), np.nan)
    n_ok = n_fail = 0

    for pos, (_, row) in enumerate(df.iterrows()):
        tok = str(row[tok_col]).lower()
        if tok not in dp2_lookup.index or pd.isna(row.get(dir_col)):
            n_fail += 1; continue
        drow = dp2_lookup.loc[tok]
        if epi_col is not None:
            epi_ris = parse_index_list(drow[epi_col])
        else:
            ws = int(row[ws_col]) if (ws_col and pd.notna(row.get(ws_col))) else 0
            epi_ris = [ws + i for i in parse_index_list(drow[epi_mpnn])]
        loaded = load_af3_heavy(Path(str(row[dir_col])))
        if loaded is None or not epi_ris:
            n_fail += 1; continue
        coords, res_idx, ca_mask = loaded
        if max(epi_ris) >= int(res_idx.max()) + 1:
            n_fail += 1; continue
        epi_set = set(epi_ris)

        epi_ca = coords[ca_mask & np.isin(res_idx, list(epi_set))]
        if len(epi_ca) < 3:
            n_fail += 1; continue
        centroid = epi_ca.mean(axis=0)
        _, _, Vt = np.linalg.svd(epi_ca - centroid)
        normal = Vt[-1]
        all_ca = coords[ca_mask]
        if np.dot(normal, all_ca.mean(axis=0) - centroid) > 0:
            normal = -normal

        scaf_mask = ~np.isin(res_idx, list(epi_set))
        scaf_coords = coords[scaf_mask]
        scaf_res = res_idx[scaf_mask]
        n_res = int(res_idx.max()) + 1
        if len(scaf_coords) < 1:
            n_fail += 1; continue
        scaf_tree = cKDTree(scaf_coords)

        own = str(row[id_col]).lower() if id_col else None
        per_probe = []
        for pid, (ab_c, appr) in probes.items():
            if own is not None and pid.lower() == own:
                continue           # no leakage: skip the design's own target
            per_probe.append(fab_clash(scaf_tree, scaf_res, n_res,
                                       ab_c, appr, centroid, normal,
                                       spins, args.clash_cutoff))
        if not per_probe:
            n_fail += 1; continue
        fab_min[pos] = min(per_probe)
        fab_mean[pos] = float(np.mean(per_probe))
        fab_np[pos] = len(per_probe)
        n_ok += 1
        if (pos + 1) % 200 == 0:
            print(f"  {pos+1:,}  ok={n_ok:,} fail={n_fail:,}")

    df["fab_min_clash"] = fab_min
    df["fab_mean_clash"] = fab_mean
    df["fab_n_probes"] = fab_np
    out = Path(args.out_csv); out.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out, index=False)
    print(f"\nWrote {out}  (ok={n_ok:,} fail={n_fail:,})")

    tc = first_present(df.columns, TRUE_CLASH)
    if tc is not None:
        truth = pd.to_numeric(df[tc], errors="coerce")
        print(f"\nCorrelation of each proxy vs true clash ({tc}):")
        print("  proxy                 Pearson   Spearman")
        ref = ["fab_min_clash", "fab_mean_clash"]
        if "cone_p10_clash" in df.columns:
            ref.append("cone_p10_clash")
        if "cylinder_ca_clashes" in df.columns:
            ref.append("cylinder_ca_clashes")
        for c in ref:
            p = pd.to_numeric(df[c], errors="coerce")
            m = p.notna() & truth.notna()
            if m.sum() < 10:
                continue
            tag = "  <- cylinder" if c == "cylinder_ca_clashes" else ""
            print(f"  {c:20s}  {p[m].corr(truth[m]):+.3f}    "
                  f"{p[m].corr(truth[m], method='spearman'):+.3f}{tag}")


if __name__ == "__main__":
    main()
