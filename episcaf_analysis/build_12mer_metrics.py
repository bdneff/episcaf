#!/usr/bin/env python3
"""
build_12mer_metrics.py

Build the composite-score inputs for the no-antibody 12-mer scaffold set, in one pass.
Per AF3 design it emits:

  epitope_chunk_rmsd  : AF3 epitope CA vs crystal epitope CA (Kabsch RMSD)         [gate + 0.35]
  overall_rmsd        : AF3 backbone CA vs the design's MPNN PDB CA (Kabsch RMSD)  [0.15]
  epitope_pae         : mean PAE over the epitope tile residues                    [0.25 term]
  scaffold_pae        : mean PAE over non-epitope residues
  mean_pae            : mean PAE over all residues
  ptm                 : AF3 summary_confidences ptm                               [0.10]
  cylinder_ca_clashes : scaffold CAs in the antibody-approach cylinder            [0.15, plain]
  cylinder_native_aware : same, minus residues sitting in native-antigen volume   [0.15, adopted ed=1]

No antibody / no is_pass: this is the ranking deployment of the validated composite.
Epitope geometry comes from contigs.csv (scaffold positions [left_len, left_len+L),
native epitope = tile_start_resid..tile_end_resid in the antigen chain).

    python build_12mer_metrics.py \
        --contigs_csv 01_contigs/contigs.csv \
        --af3_out_dir 05_af3/outputs \
        --mpnn_dir    /home/bneff/rfd3/run_12mer_scaffolding/04_mpnn \
        --crystal_dir data \
        --out_csv     06_score/metrics_12mer.csv \
        --exclude_dist 1.0 --limit 20
"""
import argparse
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
import gemmi

import native_cylinder_core as C   # cylinder_frame, inside_cylinder, native_aware_scaffold_count

NAME_RE = re.compile(
    r"^(?P<ag>[^_]+)_res(?P<res>\d+)_split(?P<split>\d+)_model(?P<model>\d+)_fixed_dldesign_(?P<dl>\d+)$",
    re.IGNORECASE)


# ----------------------------------------------------------------------------- #
# small structure helpers
# ----------------------------------------------------------------------------- #
def kabsch_fit(P, Q):
    """R,t mapping P->Q (so (R@P.T).T + t ~= Q) and the fit RMSD."""
    Pc, Qc = P - P.mean(0), Q - Q.mean(0)
    U, _, Vt = np.linalg.svd(Pc.T @ Qc)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    t = Q.mean(0) - R @ P.mean(0)
    aligned = (R @ P.T).T + t
    rmsd = float(np.sqrt(np.mean(np.sum((aligned - Q) ** 2, axis=1))))
    return R, t, rmsd


def chain_ca(chain):
    """CA coords for every residue in a gemmi chain, in order."""
    xs = []
    for res in chain:
        ca = res.find_atom("CA", "*")
        if ca is not None:
            xs.append([ca.pos.x, ca.pos.y, ca.pos.z])
    return np.asarray(xs, float)


def first_chain(model, prefer="A"):
    if prefer in [c.name for c in model]:
        return model[prefer]
    return model[0]


def load_af3(cif_path):
    """AF3 model: CA array (all residues) + heavy-atom xyz (unused now, kept simple)."""
    st = gemmi.read_structure(str(cif_path))
    ch = first_chain(st[0], "A")
    return chain_ca(ch)


def load_mpnn_ca(pdb_path):
    st = gemmi.read_structure(str(pdb_path))
    return chain_ca(first_chain(st[0], "A"))


def load_crystal(pdb_path, chain_id):
    """Return per-residue (seqid -> CA xyz) dict and all heavy-atom coords of the chain."""
    st = gemmi.read_structure(str(pdb_path))
    st.remove_alternative_conformations()
    model = st[0]
    ch = model[chain_id] if chain_id in [c.name for c in model] else model[0]
    ca_by_resid, heavy = {}, []
    for res in ch:
        ca = res.find_atom("CA", "*")
        if ca is not None:
            ca_by_resid[res.seqid.num] = [ca.pos.x, ca.pos.y, ca.pos.z]
        for atom in res:
            if atom.element.name != "H":
                heavy.append([atom.pos.x, atom.pos.y, atom.pos.z])
    return ca_by_resid, np.asarray(heavy, float)


def load_pae(conf_json):
    with open(conf_json) as fh:
        d = json.load(fh)
    for key in ("pae", "predicted_aligned_error"):
        if key in d:
            pae = np.asarray(d[key], float)
            if pae.ndim == 2:
                return pae
    return None


# ----------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--contigs_csv", required=True)
    ap.add_argument("--af3_out_dir", required=True)
    ap.add_argument("--mpnn_dir", required=True)
    ap.add_argument("--crystal_dir", required=True)
    ap.add_argument("--out_csv", required=True)
    ap.add_argument("--exclude_dist", type=float, default=1.0)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--progress_every", type=int, default=500)
    args = ap.parse_args()

    crystal_dir = Path(args.crystal_dir)
    out = Path(args.out_csv); out.parent.mkdir(parents=True, exist_ok=True)

    # ---- contigs: (antigen_lower, tile_start, split) -> row info ----
    print("loading contigs ...")
    contigs = pd.read_csv(args.contigs_csv)
    cmap = {}
    for _, r in contigs.iterrows():
        key = (str(r["antigen_id"]).lower(), int(r["tile_start_resid"]), int(r["split_id"]))
        cmap[key] = dict(antigen=str(r["antigen_id"]), chain=str(r["chain"]),
                         tstart=int(r["tile_start_resid"]), tend=int(r["tile_end_resid"]),
                         tile_seq=str(r["tile_seq"]), left_len=int(r["left_len"]))

    # ---- index AF3 outputs (real dirs = have a *_model.cif) ----
    print("indexing AF3 outputs ...")
    af3_idx = {}
    for d in Path(args.af3_out_dir).iterdir():
        if d.is_dir() and next(d.glob("*_model.cif"), None):
            af3_idx[d.name.lower()] = d
    print(f"  {len(af3_idx)} AF3 designs")

    # ---- index MPNN designed PDBs (for overall_rmsd), read in place from home ----
    print("indexing MPNN PDBs ...")
    mpnn_idx = {p.stem.lower(): p for p in Path(args.mpnn_dir).rglob("*.pdb")}
    print(f"  {len(mpnn_idx)} MPNN PDBs")

    # ---- crystal antigens loaded once per antigen ----
    crystals = {}   # antigen_lower -> (ca_by_resid, heavy_xyz, chain_id)

    def get_crystal(antigen, chain_id):
        k = antigen.lower()
        if k not in crystals:
            pdb = crystal_dir / f"{antigen.upper()}.pdb"
            if not pdb.exists():
                pdb = crystal_dir / f"{antigen}.pdb"
            crystals[k] = (*load_crystal(pdb, chain_id), chain_id)
        return crystals[k]

    rows = []
    names = sorted(af3_idx)
    if args.limit > 0:
        names = names[:args.limit]
    n_ok = n_fail = 0
    for i, name in enumerate(names):
        rec = dict(token=name, af3_dir=str(af3_idx[name]), status="ok")
        try:
            m = NAME_RE.match(name)
            if not m:
                rec["status"] = "name_parse_fail"; rows.append(rec); n_fail += 1; continue
            ag, res, split = m["ag"], int(m["res"]), int(m["split"])
            rec.update(antigen=ag, res=res, split=split,
                       model=int(m["model"]), dldesign=int(m["dl"]),
                       id=f"{ag.lower()}_res{res:04d}")          # epitope id (pools splits)
            cm = cmap.get((ag.lower(), res, split))
            if cm is None:
                rec["status"] = "no_contig"; rows.append(rec); n_fail += 1; continue
            L = len(cm["tile_seq"])
            epi_ris = list(range(cm["left_len"], cm["left_len"] + L))

            # AF3 structure
            cif = next(af3_idx[name].glob("*_model.cif"))
            af3_ca = load_af3(cif)
            if max(epi_ris) >= len(af3_ca):
                rec["status"] = f"epi_oob(n={len(af3_ca)})"; rows.append(rec); n_fail += 1; continue
            af3_epi = af3_ca[epi_ris]
            scaf_mask = np.ones(len(af3_ca), bool); scaf_mask[epi_ris] = False
            scaf_ca = af3_ca[scaf_mask]

            # crystal epitope (by residue number) + native heavy atoms
            ca_by_resid, heavy, _ = get_crystal(ag, cm["chain"])
            cry_resids = [rid for rid in range(cm["tstart"], cm["tend"] + 1) if rid in ca_by_resid]
            if len(cry_resids) != L:
                rec["status"] = f"crystal_epi_incomplete({len(cry_resids)}/{L})"
                rows.append(rec); n_fail += 1; continue
            cry_epi = np.array([ca_by_resid[r] for r in cry_resids], float)

            # epitope_chunk_rmsd + transform crystal antigen into AF3 frame
            R, t, epi_rmsd = kabsch_fit(cry_epi, af3_epi)      # crystal -> AF3
            rec["epitope_chunk_rmsd"] = epi_rmsd
            native_heavy_af3 = (R @ heavy.T).T + t

            # PAE decomposition
            conf = next(af3_idx[name].glob("*_confidences.json"), None)
            pae = load_pae(conf) if conf else None
            if pae is not None:
                n = pae.shape[0]
                epi = [r for r in epi_ris if r < n]
                scaf = [r for r in range(n) if r not in set(epi)]
                rec["mean_pae"] = float(np.nanmean(pae))
                rec["epitope_pae"] = float(np.nanmean(pae[np.ix_(epi, epi)])) if epi else None
                rec["scaffold_pae"] = float(np.nanmean(pae[np.ix_(scaf, scaf)])) if scaf else None

            # ptm
            sc = next(af3_idx[name].glob("*_summary_confidences.json"), None)
            if sc:
                with open(sc) as fh:
                    rec["ptm"] = json.load(fh).get("ptm")

            # overall_rmsd vs MPNN designed PDB
            mp = mpnn_idx.get(name)
            if mp is not None:
                mca = load_mpnn_ca(mp)
                if len(mca) == len(af3_ca):
                    _, _, rec["overall_rmsd"] = kabsch_fit(mca, af3_ca)
                else:
                    rec["overall_rmsd_status"] = f"len_mismatch({len(mca)}/{len(af3_ca)})"
            else:
                rec["overall_rmsd_status"] = "no_mpnn_pdb"

            # cylinder (plain + native-aware), epitope on the AF3 frame
            base, normal = C.cylinder_frame(af3_epi, af3_ca)
            aware, plain = C.native_aware_scaffold_count(
                scaf_ca, base, normal, native_heavy_af3, exclude_dist=args.exclude_dist)
            rec["cylinder_ca_clashes"] = plain
            rec["cylinder_native_aware"] = aware

            rows.append(rec); n_ok += 1
        except Exception as e:
            rec["status"] = f"error:{type(e).__name__}:{e}"
            rows.append(rec); n_fail += 1
        if args.progress_every and (i + 1) % args.progress_every == 0:
            print(f"  {i+1}/{len(names)}  ok={n_ok} fail={n_fail}")

    df = pd.DataFrame(rows)
    df.to_csv(out, index=False)
    print(f"\nwrote {out}  (ok={n_ok} fail={n_fail} of {len(names)})")
    if n_fail:
        print("fail reasons:")
        print(df.loc[df.status != "ok", "status"].value_counts().head(10).to_string())
    ok = df[df.status == "ok"]
    if len(ok):
        for c in ("epitope_chunk_rmsd", "overall_rmsd", "epitope_pae", "ptm",
                  "cylinder_ca_clashes", "cylinder_native_aware"):
            if c in ok:
                v = pd.to_numeric(ok[c], errors="coerce")
                print(f"  {c:22s} mean={v.mean():.3f}  median={v.median():.3f}  "
                      f"nonnull={v.notna().sum()}")


if __name__ == "__main__":
    main()
