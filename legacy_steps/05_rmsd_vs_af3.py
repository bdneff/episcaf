#!/usr/bin/env python3
"""
rfd3_vs_af3_rmsd_v3.py

Purpose
-------
Compute CA RMSD between paired structures:
  - RFD3 model_X mmCIF.gz
  - AF3 prediction mmCIF (model.cif)

Your AF3 layout (as observed):
  03_af3/outputs/<stem>/<stem>/seed-42_sample-0/model.cif

This script is designed to:
- Find the correct AF3 model.cif path for each RFD3 *_model_X.cif.gz
- Parse CIF/CIF.GZ reliably with gemmi (not MDAnalysis CIF readers)
- Convert to temporary PDBs
- Use MDAnalysis to:
    * select CA atoms robustly (does NOT rely on "protein" selection)
    * align mobile -> reference by least-squares
    * compute RMSD on shared residues
- Handle common mismatch cases:
    * CA count mismatch -> try resid intersection; then sequence alignment fallback
- Never silently fail: writes status + notes for every row
- Write:
    * rmsd_vs_af3_all.csv
    * rmsd_vs_af3_best_per_run.csv (lowest RMSD model_idx per run_id)

Dependencies
------------
pip install gemmi MDAnalysis numpy pandas

Usage
-----
python rfd3_vs_af3_rmsd_v3.py \
  --rfd3_outputs_root runs/.../02_rfd3/outputs \
  --af3_outputs_root  runs/.../03_af3/outputs \
  --verbose

Notes
-----
- Selection uses "name CA" (not "protein and name CA") because protein typing often fails after CIF->PDB.
- If you want strict paired-equality behavior, set --strict_equal_ca to error on mismatch instead of mapping.
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple, Optional, Dict

import numpy as np
import pandas as pd

import gemmi
import MDAnalysis as mda
from MDAnalysis.analysis import align


MODEL_RE = re.compile(r"^(?P<stem>.*_model_(?P<mid>\d+))\.cif(?:\.gz)?$")


# ----------------------------- CIF -> PDB -----------------------------

def cif_to_pdb(cif_path: Path, pdb_path: Path) -> None:
    """
    Read .cif or .cif.gz with gemmi and write PDB.
    Keeps only model 0 if multiple models are present.
    """
    st = gemmi.read_structure(str(cif_path))  # handles .cif and .cif.gz
    st.setup_entities()
    if len(st) > 1:
        st0 = gemmi.Structure()
        st0.name = st.name
        st0.cell = st.cell
        st0.spacegroup_hm = st.spacegroup_hm
        st0.add_model(st[0])
        st = st0
    st.write_pdb(str(pdb_path))


# -------------------------- AF3 path finder ---------------------------

def find_af3_model_cif(af3_outputs_root: Path, stem: str, prefer: str = "seed-42_sample-0") -> Path:
    """
    Robustly locate AF3 model.cif for a given stem.

    Searches:
      af3_outputs_root/stem/**/model.cif

    Prefers a match containing `prefer` in the path.
    """
    base = af3_outputs_root / stem
    if not base.exists():
        raise RuntimeError(f"AF3 base folder missing: {base}")

    matches = sorted(base.glob("**/model.cif"))
    if len(matches) == 0:
        raise RuntimeError(f"No AF3 model.cif found under: {base}")

    preferred = [p for p in matches if prefer in str(p)]
    return preferred[0] if preferred else matches[0]


def diffusion_run_id_from_rfd3_path(rfd3_path: Path) -> str:
    """
    Grouping key for 'best per run' summaries.
    Adjust if your concept of a 'run' differs.
    """
    return rfd3_path.parent.name


# -------------------------- Mapping utilities -------------------------

AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
    # common variants
    "HID": "H", "HIE": "H", "HIP": "H",
    "CYX": "C", "CYM": "C",
    "MSE": "M",
}

def aa3_to_1(resname: str) -> str:
    r = (resname or "").strip().upper()
    return AA3_TO_1.get(r, "X")


def needleman_wunsch(seqA: str, seqB: str, match=2, mismatch=-1, gap=-2) -> Tuple[str, str]:
    """
    Simple global alignment returning aligned sequences with '-' gaps.
    """
    n, m = len(seqA), len(seqB)
    score = np.zeros((n + 1, m + 1), dtype=int)
    trace = np.zeros((n + 1, m + 1), dtype=np.int8)  # 0 diag, 1 up, 2 left

    for i in range(1, n + 1):
        score[i, 0] = score[i - 1, 0] + gap
        trace[i, 0] = 1
    for j in range(1, m + 1):
        score[0, j] = score[0, j - 1] + gap
        trace[0, j] = 2

    for i in range(1, n + 1):
        a = seqA[i - 1]
        for j in range(1, m + 1):
            b = seqB[j - 1]
            diag = score[i - 1, j - 1] + (match if a == b else mismatch)
            up = score[i - 1, j] + gap
            left = score[i, j - 1] + gap
            best = max(diag, up, left)
            score[i, j] = best
            trace[i, j] = 0 if best == diag else (1 if best == up else 2)

    i, j = n, m
    outA, outB = [], []
    while i > 0 or j > 0:
        t = trace[i, j] if (i > 0 and j > 0) else (1 if i > 0 else 2)
        if t == 0:
            outA.append(seqA[i - 1]); outB.append(seqB[j - 1])
            i -= 1; j -= 1
        elif t == 1:
            outA.append(seqA[i - 1]); outB.append("-")
            i -= 1
        else:
            outA.append("-"); outB.append(seqB[j - 1])
            j -= 1

    return "".join(reversed(outA)), "".join(reversed(outB))


def build_ca_trace(u: mda.Universe) -> Tuple[np.ndarray, List[int], List[str]]:
    """
    Build a CA coordinate array and parallel residue identifiers for mapping.

    Returns:
      coords: (N,3)
      resids: list of resid (int) per CA
      resnames: list of resname (str) per CA
    """
    ca = u.select_atoms("name CA")
    if len(ca) == 0:
        raise RuntimeError("Selection 'name CA' is empty (no CA atoms)")

    coords = ca.positions.copy()
    # Each CA belongs to a residue; use atoms.resids/resnames aligned with atoms
    resids = ca.resids.tolist()
    resnames = [str(x) for x in ca.resnames]
    return coords, resids, resnames


def map_by_resid_intersection(ref_resids: List[int], mob_resids: List[int], min_shared: int) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """
    Map indices by intersecting residue numbers (resids). Works if numbering overlaps.
    Returns (ref_idx, mob_idx) arrays or None.
    """
    ref_map: Dict[int, int] = {}
    for i, r in enumerate(ref_resids):
        # keep first occurrence if any duplicates
        if r not in ref_map:
            ref_map[r] = i
    mob_map: Dict[int, int] = {}
    for i, r in enumerate(mob_resids):
        if r not in mob_map:
            mob_map[r] = i

    shared = sorted(set(ref_map.keys()).intersection(set(mob_map.keys())))
    if len(shared) < min_shared:
        return None

    ref_idx = np.array([ref_map[r] for r in shared], dtype=int)
    mob_idx = np.array([mob_map[r] for r in shared], dtype=int)
    return ref_idx, mob_idx


def map_by_sequence_alignment(ref_resnames: List[str], mob_resnames: List[str], min_shared: int) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    """
    Map indices by aligning sequences derived from CA residue names.
    Returns (ref_idx, mob_idx) arrays or None.
    """
    ref_seq = "".join(aa3_to_1(rn) for rn in ref_resnames)
    mob_seq = "".join(aa3_to_1(rn) for rn in mob_resnames)

    a_ref, a_mob = needleman_wunsch(ref_seq, mob_seq)

    ref_idx = []
    mob_idx = []
    i_ref = -1
    i_mob = -1
    for ar, am in zip(a_ref, a_mob):
        if ar != "-":
            i_ref += 1
        if am != "-":
            i_mob += 1
        if ar != "-" and am != "-":
            ref_idx.append(i_ref)
            mob_idx.append(i_mob)

    if len(ref_idx) < min_shared:
        return None

    return np.array(ref_idx, dtype=int), np.array(mob_idx, dtype=int)


# ----------------------------- RMSD core ------------------------------

def ca_rmsd_pair(
    ref_pdb: Path,
    mob_pdb: Path,
    min_ca: int = 20,
    strict_equal_ca: bool = False,
    verbose: bool = False,
) -> Tuple[float, int, str, str]:
    """
    Compute CA RMSD after least-squares alignment.

    Returns:
      rmsd_A, n_shared, mapping_strategy, diag_notes

    Strategy:
      - Build CA traces (coords + resids + resnames) for both
      - If CA counts equal and strict_equal_ca=True or they match: use index match
      - Else map:
          1) resid intersection (fast)
          2) sequence alignment fallback
    """
    ref_u = mda.Universe(str(ref_pdb))
    mob_u = mda.Universe(str(mob_pdb))

    ref_coords, ref_resids, ref_resnames = build_ca_trace(ref_u)
    mob_coords, mob_resids, mob_resnames = build_ca_trace(mob_u)

    n_ref = ref_coords.shape[0]
    n_mob = mob_coords.shape[0]

    if verbose:
        print(f"  Loaded ref: atoms={ref_u.atoms.n_atoms} residues={ref_u.residues.n_residues} CA={n_ref}")
        print(f"  Loaded mob: atoms={mob_u.atoms.n_atoms} residues={mob_u.residues.n_residues} CA={n_mob}")

    if n_ref < min_ca or n_mob < min_ca:
        raise RuntimeError(f"Too few CA atoms: ref_ca={n_ref} mob_ca={n_mob} (min_ca={min_ca})")

    if strict_equal_ca and n_ref != n_mob:
        raise RuntimeError(f"CA count mismatch (strict): ref_ca={n_ref} mob_ca={n_mob}")

    # Mapping
    if n_ref == n_mob:
        ref_idx = np.arange(n_ref, dtype=int)
        mob_idx = np.arange(n_mob, dtype=int)
        strategy = "index_match"
    else:
        # 1) residue id intersection
        mapped = map_by_resid_intersection(ref_resids, mob_resids, min_shared=min_ca)
        if mapped is not None:
            ref_idx, mob_idx = mapped
            strategy = "resid_intersection"
        else:
            # 2) sequence alignment fallback
            mapped = map_by_sequence_alignment(ref_resnames, mob_resnames, min_shared=min_ca)
            if mapped is None:
                raise RuntimeError(
                    f"Could not map CA traces (ref_ca={n_ref}, mob_ca={n_mob}); "
                    f"resid overlap too small and sequence alignment shared < {min_ca}"
                )
            ref_idx, mob_idx = mapped
            strategy = "sequence_alignment"

    # Build AtomGroups for alignment using the mapped indices.
    # We align using MDAnalysis (least squares) on the selected CA subset.
    ref_ca = ref_u.select_atoms("name CA")[ref_idx]
    mob_ca = mob_u.select_atoms("name CA")[mob_idx]

    if len(ref_ca) != len(mob_ca) or len(ref_ca) == 0:
        raise RuntimeError(f"Mapped CA mismatch: ref={len(ref_ca)} mob={len(mob_ca)}")

    n_shared = len(ref_ca)

    # Align mobile Universe to reference using those CA atoms
    # We can't pass arbitrary AtomGroups directly to alignto's 'select' string,
    # so we do the math ourselves with Kabsch using the mapped coordinates.
    # But we still keep MDAnalysis for IO + selections.
    ref_xyz = ref_ca.positions.copy()
    mob_xyz = mob_ca.positions.copy()

    # Kabsch alignment
    ref_cent = ref_xyz.mean(axis=0)
    mob_cent = mob_xyz.mean(axis=0)
    ref0 = ref_xyz - ref_cent
    mob0 = mob_xyz - mob_cent

    C = mob0.T @ ref0
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(V @ Wt))
    D = np.diag([1.0, 1.0, d])
    R = V @ D @ Wt

    mob_aligned = mob0 @ R + ref_cent
    diff = mob_aligned - ref_xyz
    rmsd = np.sqrt((diff * diff).sum() / n_shared)

    diag = ""
    if n_ref != n_mob:
        diag = f"CA mismatch raw: ref_ca={n_ref} mob_ca={n_mob}; mapped n={n_shared} via {strategy}"

    return float(rmsd), int(n_shared), strategy, diag


# -------------------------------- Main --------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rfd3_outputs_root", required=True, type=Path,
                    help="Path to .../02_rfd3/outputs")
    ap.add_argument("--af3_outputs_root", required=True, type=Path,
                    help="Path to .../03_af3/outputs")
    ap.add_argument("--out_all", default="rmsd_vs_af3_all.csv", type=Path)
    ap.add_argument("--out_best", default="rmsd_vs_af3_best_per_run.csv", type=Path)
    ap.add_argument("--tmp_pdb_dir", default="tmp_pdb_for_rmsd", type=Path)
    ap.add_argument("--min_ca", default=20, type=int)
    ap.add_argument("--prefer_seed", default="seed-42_sample-0")
    ap.add_argument("--strict_equal_ca", action="store_true",
                    help="If set, fail when ref and mob CA counts differ (no mapping fallback).")
    ap.add_argument("--verbose", action="store_true")
    args = ap.parse_args()

    args.tmp_pdb_dir.mkdir(parents=True, exist_ok=True)

    rfd3_files = sorted(args.rfd3_outputs_root.glob("**/*_model_*.cif.gz"))
    if not rfd3_files:
        print(f"No RFD3 *_model_*.cif.gz found under {args.rfd3_outputs_root}", file=sys.stderr)
        sys.exit(2)

    rows = []
    ok = fail = 0

    for rfd3_cif_gz in rfd3_files:
        run_id = diffusion_run_id_from_rfd3_path(rfd3_cif_gz)

        m = MODEL_RE.match(rfd3_cif_gz.name)
        if not m:
            rows.append({
                "status": "fail",
                "run_id": run_id,
                "model_idx": np.nan,
                "rfd3_cif_gz": str(rfd3_cif_gz),
                "af3_cif": "",
                "n_ca_shared": np.nan,
                "rmsd_A": np.nan,
                "mapping_strategy": "",
                "notes": f"Unrecognized filename pattern: {rfd3_cif_gz.name}",
            })
            fail += 1
            continue

        stem = m.group("stem")  # includes _model_X
        model_idx = int(m.group("mid"))

        try:
            af3_cif = find_af3_model_cif(args.af3_outputs_root, stem, prefer=args.prefer_seed)

            # Convert to PDBs
            ref_pdb = args.tmp_pdb_dir / f"AF3__{stem}.pdb"
            mob_pdb = args.tmp_pdb_dir / f"RFD3__{stem}.pdb"

            cif_to_pdb(af3_cif, ref_pdb)
            cif_to_pdb(rfd3_cif_gz, mob_pdb)

            if args.verbose:
                print(f"\n[{stem}]")
                print(f"  RFD3: {rfd3_cif_gz}")
                print(f"  AF3 : {af3_cif}")
                print(f"  PDB : ref={ref_pdb.name} mob={mob_pdb.name}")

            rmsd, n_shared, strategy, diag = ca_rmsd_pair(
                ref_pdb,
                mob_pdb,
                min_ca=args.min_ca,
                strict_equal_ca=args.strict_equal_ca,
                verbose=args.verbose,
            )

            rows.append({
                "status": "ok",
                "run_id": run_id,
                "model_idx": model_idx,
                "rfd3_cif_gz": str(rfd3_cif_gz),
                "af3_cif": str(af3_cif),
                "n_ca_shared": n_shared,
                "rmsd_A": rmsd,
                "mapping_strategy": strategy,
                "notes": diag,
            })
            ok += 1

            if args.verbose:
                print(f"  -> OK: rmsd={rmsd:.3f} Å, n_shared={n_shared}, strategy={strategy}")
                if diag:
                    print(f"     note: {diag}")

        except Exception as e:
            rows.append({
                "status": "fail",
                "run_id": run_id,
                "model_idx": model_idx,
                "rfd3_cif_gz": str(rfd3_cif_gz),
                "af3_cif": "",
                "n_ca_shared": np.nan,
                "rmsd_A": np.nan,
                "mapping_strategy": "",
                "notes": str(e),
            })
            fail += 1
            if args.verbose:
                print(f"  -> FAIL: {e}", file=sys.stderr)

    df = pd.DataFrame(rows)
    df.to_csv(args.out_all, index=False)

    ok_df = df[df["status"] == "ok"].copy()
    if len(ok_df) > 0:
        best = ok_df.sort_values(["run_id", "rmsd_A"]).groupby("run_id", as_index=False).first()
        best.to_csv(args.out_best, index=False)
    else:
        pd.DataFrame(columns=df.columns).to_csv(args.out_best, index=False)

    print(f"\nTotal rows: {len(df)}")
    print(f"OK rows: {ok}")
    print(f"FAIL rows: {fail}")
    print(f"Wrote: {args.out_all}")
    print(f"Wrote: {args.out_best}")

    if fail:
        print("\nTop failure reasons:")
        top = (df[df.status == "fail"]["notes"].value_counts().head(10))
        for reason, count in top.items():
            print(f"  {count:4d}  {reason}")


if __name__ == "__main__":
    main()
