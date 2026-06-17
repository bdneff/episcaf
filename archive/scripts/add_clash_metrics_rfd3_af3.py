#!/usr/bin/env python3
from __future__ import annotations

import argparse
import ast
import gzip
from pathlib import Path
from tempfile import NamedTemporaryFile

import numpy as np
import pandas as pd
import gemmi
import MDAnalysis as mda
from MDAnalysis.analysis.align import alignto

# -----------------------------
# Helpers: read CIF(.gz) -> PDB temp -> MDAnalysis Universe
# -----------------------------
def _read_gemmi_structure(p: Path) -> gemmi.Structure:
    s = str(p)
    if s.endswith(".gz"):
        with gzip.open(p, "rt") as f:
            doc = gemmi.cif.read_string(f.read())
        return gemmi.make_structure_from_block(doc.sole_block())
    if p.suffix.lower() == ".cif":
        doc = gemmi.cif.read(str(p))
        return gemmi.make_structure_from_block(doc.sole_block())
    # fallback (PDB etc.)
    return gemmi.read_structure(str(p))

def universe_from_structure_file(p: Path) -> mda.Universe:
    """
    MDAnalysis can't read CIF reliably in many installs.
    Convert via gemmi -> temporary PDB, then load into MDAnalysis.
    """
    p = Path(p)
    if p.suffix.lower() in (".cif",) or str(p).endswith(".cif.gz"):
        st = _read_gemmi_structure(p)
        with NamedTemporaryFile(suffix=".pdb", delete=True) as tmp:
            st.write_pdb(tmp.name)
            return mda.Universe(tmp.name)
    # PDB etc
    return mda.Universe(str(p))

def pick_antigen_atoms(u: mda.Universe) -> mda.core.groups.AtomGroup:
    """
    Antigen is chain/segid A in our conventions.
    Try segid A first, then chainID A.
    """
    ag = u.select_atoms("segid A")
    if len(ag) == 0:
        ag = u.select_atoms("chainid A")
    if len(ag) == 0:
        # last resort: just protein
        ag = u.select_atoms("protein")
    return ag

def pick_true_antibody_atoms(u_true: mda.Universe) -> mda.core.groups.AtomGroup:
    ab = u_true.select_atoms("segid B or segid C")
    if len(ab) == 0:
        ab = u_true.select_atoms("chainid B or chainid C")
    return ab

def ensure_list(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return None
    if isinstance(x, list):
        return x
    # dp2 parquet will usually already be list; CSV sometimes stringifies
    if isinstance(x, str):
        try:
            return list(ast.literal_eval(x))
        except Exception:
            pass
    return None

def compute_clash_resindices(
    pred_path: Path,
    true_pdb: Path,
    pred_epi_ris0: list[int],
    true_epi_ris0: list[int],
    cutoff: float = 4.0,
) -> list[int] | None:
    """
    Lawson-style:
      - align pred epitope chunk backbone to true epitope chunk backbone
      - merge pred antigen + true antibody
      - unintended residues = antigen residues NOT in epitope chunk
      - find antibody heavy atoms within cutoff Å of unintended
      - return antibody residue resindices (in the merged universe)
    """
    u_pred = universe_from_structure_file(pred_path)
    u_true = mda.Universe(str(true_pdb))

    predA = pick_antigen_atoms(u_pred)
    trueA = pick_antigen_atoms(u_true)
    trueAB = pick_true_antibody_atoms(u_true)
    if len(trueAB) == 0:
        raise ValueError("Could not find antibody atoms as segid/chain B or C in true pdb.")

    # Align using backbone atoms of epitope chunk
    alignto(
        predA.residues[pred_epi_ris0].atoms.select_atoms("backbone"),
        trueA.residues[true_epi_ris0].atoms.select_atoms("backbone"),
    )

    # Merge pred antigen (all atoms) + true antibody (all atoms)
    u_comb = mda.Merge(u_pred.atoms, trueAB.atoms)

    # Identify unintended residues on antigen (segid/chain A) in merged universe
    combA = pick_antigen_atoms(u_comb)

    mask = np.full(len(combA.residues), False)
    mask[pred_epi_ris0] = True
    unintended_res = combA.residues[~mask]

    # Antibody heavy atoms within cutoff of unintended residues
    sel = u_comb.select_atoms(
        f"((segid B or segid C) and (not name H*)) and around {cutoff} group unintended",
        unintended=unintended_res.atoms,
    )
    # return residue resindices (0-based)
    return list(sel.residues.resindices)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metrics_csv", required=True, help="Input CSV (your combined metrics; must have id,tok,pred,rfd3_path,af3_path)")
    ap.add_argument("--dp2_parquet", required=True, help="dp2 parquet or subset parquet (needs token + epitope_chunk_resindices + scaffolded_epitope_chunk_resindices)")
    ap.add_argument("--true_dir", required=True, help="Directory with true complex PDBs named <id>.pdb")
    ap.add_argument("--out_csv", required=True, help="Output CSV with clash cols appended")
    ap.add_argument("--cutoff", type=float, default=4.0, help="Å cutoff (Lawson used 4)")
    ap.add_argument("--limit", type=int, default=0, help="If >0, limit rows for debugging")
    args = ap.parse_args()

    metrics_csv = Path(args.metrics_csv)
    dp2_parquet = Path(args.dp2_parquet)
    true_dir = Path(args.true_dir)

    df = pd.read_csv(metrics_csv)
    if args.limit and args.limit > 0:
        df = df.head(args.limit).copy()

    required = ["id","tok","pred"]
    miss = [c for c in required if c not in df.columns]
    if miss:
        raise SystemExit(f"metrics_csv missing required cols: {miss}")

    # load dp2 + build mapping tok -> (true_epi_ris0, pred_epi_ris0)
    dp2 = pd.read_parquet(dp2_parquet)
    need_dp2 = ["assay_scaffolded_epitope_id","epitope_chunk_resindices","scaffolded_epitope_chunk_resindices"]
    miss2 = [c for c in need_dp2 if c not in dp2.columns]
    if miss2:
        raise SystemExit(f"dp2 parquet missing required cols: {miss2}")

    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()
    tok2true = {}
    for _, r in dp2.iterrows():
        tok = str(r["assay_scaffolded_epitope_id"]).lower()
        tok2true[tok] = (
            ensure_list(r["epitope_chunk_resindices"]),
            ensure_list(r["scaffolded_epitope_chunk_resindices"]),
        )

    # output columns
    df["rfd3_n_clash_res"] = 0
    df["rfd3_clash_resindices"] = None
    df["af3_n_clash_res"] = 0
    df["af3_clash_resindices"] = None
    df["clash_status"] = "ok"

    epi_fail = 0
    missing_tok = 0

    for i, row in df.iterrows():
        tok = str(row["tok"]).lower()
        pid = str(row["id"])
        true_pdb = true_dir / f"{pid}.pdb"

        if tok not in tok2true:
            df.at[i, "clash_status"] = "missing_tok_in_dp2"
            missing_tok += 1
            continue

        true_epi_ris0, pred_epi_ris0 = tok2true[tok]
        if true_epi_ris0 is None or pred_epi_ris0 is None:
            df.at[i, "clash_status"] = "missing_epitope_indices"
            epi_fail += 1
            continue

        try:
            # RFD3 vs true antibody clashes
            rfd3_root = Path(args.metrics_csv).parent.parent / "02_rfd3" / "outputs"
            rfd3_path = next(rfd3_root.glob(f"{tok}*0_model_{row["pred"]}.cif.gz"), None)
            rfd3_clash = compute_clash_resindices(
                pred_path=rfd3_path,
                true_pdb=true_pdb,
                pred_epi_ris0=pred_epi_ris0,
                true_epi_ris0=true_epi_ris0,
                cutoff=args.cutoff,
            )
            df.at[i, "rfd3_clash_resindices"] = rfd3_clash
            df.at[i, "rfd3_n_clash_res"] = int(len(rfd3_clash)) if rfd3_clash is not None else 0

            # AF3 vs true antibody clashes
            af3_root = Path(args.metrics_csv).parent.parent / "03_af3" / "outputs"
            af3_dir = next(af3_root.glob(f"{tok}*0_model_{row["pred"]}"), None)
            af3_path = (af3_dir / f"{af3_dir.name}_model.cif") if af3_dir else None
            # For AF3, the antigen in the AF3 complex includes extra stuff; but our pred_epi_ris0 is in scaffold antigen coordinates.
            # We already aligned on the epitope chunk backbone, so we just use the *mapped* chunk indices if you have them.
            # If your CSV already contains af3_window_start, map indices: af3_epi = [ws + i for i in pred_epi_ris0]
            af3_epi = pred_epi_ris0
            if "af3_window_start" in df.columns and pd.notna(row["af3_window_start"]):
                ws = int(row["af3_window_start"])
                af3_epi = [ws + int(j) for j in pred_epi_ris0]

            af3_clash = compute_clash_resindices(
                pred_path=af3_path,
                true_pdb=true_pdb,
                pred_epi_ris0=af3_epi,
                true_epi_ris0=true_epi_ris0,
                cutoff=args.cutoff,
            )
            df.at[i, "af3_clash_resindices"] = af3_clash
            df.at[i, "af3_n_clash_res"] = int(len(af3_clash)) if af3_clash is not None else 0

        except Exception:
            df.at[i, "clash_status"] = "clash_fail"
            epi_fail += 1
            continue

        if (i+1) % 200 == 0:
            print(f"[progress] {i+1}/{len(df)} | missing_tok={missing_tok} fail={epi_fail}")

    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.out_csv, index=False)

    print("DONE")
    print("rows:", len(df))
    print("missing_tok_in_dp2:", missing_tok)
    print("fail:", epi_fail)
    print("wrote:", args.out_csv)

if __name__ == "__main__":
    main()
