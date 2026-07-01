#!/usr/bin/env python3
"""Why does the cylinder flag clashes that are not real steric clashes? (8pww probe)

For one assayed design, this locates the scaffold Calpha atoms the cylinder counts as
"clashing", then transforms the REAL antibody (from the AbDb complex, chains B/C) into the
same frame and measures how far each flagged atom is from the nearest antibody heavy atom.
The story behind af3_n_clash_res=0 with cylinder~10 is that the flagged atoms sit inside the
cylinder VOLUME but far from where the antibody actually is -- the cylinder (R=16, H=40 A) is a
coarse over-approximation of the antibody footprint. This prints that gap in numbers and writes
the aligned antibody + the flagged-atom residue list so you can see it in VMD.

Runs on Gemini (needs gemmi + the structures). Example for the top 8pww outlier (DP2_0804):
    conda activate ~/rfd3/env/rfd3_py312
    B=/tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/sourced_antibody_v1/no_antibody
    python3 scripts/cylinder_fp_probe.py \
        --worklist    results/assayed_cylinder_worklist.csv \
        --library_member DP2_0804 \
        --dp2_parquet $B/assay_scaffold_simple_metrics_403.parquet \
        --native_dir  /tgen_labs/altin/alphafold3/workspace/episcaf-experiments/data/abdb/complex_pdbfiles/cleaned \
        --af3_root    $B/af3_predictions \
        --out_dir     results/cylinder_fp/DP2_0804
"""
from __future__ import annotations
import argparse
from pathlib import Path
import sys
import numpy as np
import gemmi
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parent))
from dp3_native_cylinder import (  # identical primitives; no re-implementation  # noqa: E402
    read_gemmi, chain_by_name, load_af3, cylinder_frame, inside_cylinder, match_epitope,
    parse_index_list, first_present, DP2_EPI_COLS, RADIUS, HEIGHT, THREE2ONE)

CLASH_CUT = 4.0  # A; same definition as af3_n_clash_res (scaffold heavy vs antibody heavy)


CARVE_DIST = 1.0  # native-aware carve distance used in production (docs/CYLINDER_PARAMS.md)


def load_complex(pdb: Path):
    """Antigen chain A (CA + seq; heavy atoms + their CA-residue index) and the antibody
    (every non-A chain's heavy atoms)."""
    st = read_gemmi(pdb)
    ag = chain_by_name(st[0], "A")
    aca, aseq, ag_heavy, ag_heavy_resi = [], [], [], []
    for res in ag:
        a = res.find_atom("CA", altloc="*")
        if not a:
            continue
        ri = len(aca)
        aca.append([a.pos.x, a.pos.y, a.pos.z]); aseq.append(THREE2ONE.get(res.name.upper(), "X"))
        for at in res:
            if at.element != gemmi.Element("H"):
                ag_heavy.append([at.pos.x, at.pos.y, at.pos.z]); ag_heavy_resi.append(ri)
    ab_heavy = []
    for ch in st[0]:
        if ch.name == "A":
            continue
        for res in ch:
            for at in res:
                if at.element != gemmi.Element("H"):
                    ab_heavy.append([at.pos.x, at.pos.y, at.pos.z])
    return (np.asarray(aca, float), "".join(aseq), np.asarray(ag_heavy, float),
            np.asarray(ag_heavy_resi, int), np.asarray(ab_heavy, float))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--worklist", required=True)
    ap.add_argument("--library_member", required=True)
    ap.add_argument("--dp2_parquet", required=True)
    ap.add_argument("--native_dir", required=True)
    ap.add_argument("--af3_root", required=True)
    ap.add_argument("--out_dir", default="")
    args = ap.parse_args()

    import pandas as pd
    wl = pd.read_csv(args.worklist)
    row = wl[wl.library_member == args.library_member]
    if row.empty:
        sys.exit(f"{args.library_member} not in worklist")
    row = row.iloc[0]
    epi_hash = str(row.assay_scaffolded_epitope_id).lower()
    epi_id = str(row.id)

    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"] = dp2.assay_scaffolded_epitope_id.astype(str).str.lower()
    epic = first_present(dp2.columns, DP2_EPI_COLS)
    epi_ris = parse_index_list(
        dp2.set_index("assay_scaffolded_epitope_id").loc[epi_hash, epic])

    af3_dir = Path(args.af3_root) / epi_hash
    native_pdb = Path(args.native_dir) / f"{epi_id}.pdb"
    print(f"design={args.library_member}  epitope={epi_id}  af3={af3_dir}\n"
          f"native complex={native_pdb}")

    ca, res_idx, seq = load_af3(af3_dir)
    epi_pos = [p for p, ri in enumerate(res_idx) if ri in set(epi_ris)]
    epi_ca = ca[epi_pos]; epi_seq = "".join(seq[p] for p in epi_pos)
    centroid, normal, base = cylinder_frame(epi_ca, ca)

    scaf_pos = np.array([p for p in range(len(ca)) if res_idx[p] not in set(epi_ris)])
    scaf_ca = ca[scaf_pos]
    ins = inside_cylinder(scaf_ca, base, normal)
    flagged_xyz = scaf_ca[ins]
    flagged_res = [int(res_idx[scaf_pos[i]]) for i in np.where(ins)[0]]  # 0-based design residues
    print(f"\ncylinder: R={RADIUS} H={HEIGHT} A ; scaffold CAs inside = {int(ins.sum())} (the plain count)")

    aca, aseq, ag_heavy, ag_heavy_resi, ab_heavy = load_complex(native_pdb)
    m = match_epitope(epi_seq, epi_ca, aseq, aca)
    if m is None:
        sys.exit("could not align design epitope to native antigen")
    R, t, nepi_set, rmsd = m
    ab_al = (R @ ab_heavy.T).T + t          # antibody, in the design frame
    ag_al = (R @ ag_heavy.T).T + t          # native antigen heavy, in the design frame
    ag_nonepi = ag_al[~np.isin(ag_heavy_resi, list(nepi_set))]   # non-epitope antigen volume
    print(f"epitope alignment RMSD = {rmsd:.2f} A ; antibody heavy atoms = {len(ab_al)}")

    # native-aware carve: which flagged CAs sit in native-antigen volume (carved) vs not (survive)
    d_native, _ = cKDTree(ag_nonepi).query(flagged_xyz, k=1) if len(ag_nonepi) else (np.full(len(flagged_xyz), 1e9), None)
    carved = d_native <= CARVE_DIST
    survive = ~carved
    print(f"native-aware carve ({CARVE_DIST} A): {int(survive.sum())} survive (the native-aware "
          f"count), {int(carved.sum())} carved (sit on the native antigen)")

    # the key question: how close are the cylinder-flagged scaffold atoms to the REAL antibody?
    tree = cKDTree(ab_al)
    d_flag, _ = tree.query(flagged_xyz, k=1)
    n_real = int((d_flag <= CLASH_CUT).sum())
    ab_in_cyl = int(inside_cylinder(ab_al, base, normal).sum())
    print(f"\n--- why the flags are not steric clashes ---")
    print(f"flagged scaffold CAs within {CLASH_CUT} A of the real antibody (= real clashes): {n_real}/{len(d_flag)}")
    print(f"distance of flagged CAs to nearest antibody atom (A): "
          f"min={d_flag.min():.1f} median={np.median(d_flag):.1f} max={d_flag.max():.1f}")
    print(f"antibody heavy atoms that actually fall INSIDE the cylinder: {ab_in_cyl}/{len(ab_al)} "
          f"({100*ab_in_cyl/len(ab_al):.0f}%) -- how much of the cylinder the antibody really fills")
    # axial position (up the approach normal) of flagged atoms vs antibody, to show the mismatch
    ax_flag = (flagged_xyz - base) @ normal
    ax_ab = (ab_al - base) @ normal
    print(f"axial position along the approach normal (A):  flagged CAs {ax_flag.min():.0f}-{ax_flag.max():.0f}"
          f" ; antibody atoms {ax_ab.min():.0f}-{ax_ab.max():.0f}  (cylinder spans 0-{HEIGHT:.0f})")

    if args.out_dir:
        out = Path(args.out_dir); out.mkdir(parents=True, exist_ok=True)
        def write_points_pdb(path, coords, chain, resname):
            with open(path, "w") as fh:
                for i, (x, y, z) in enumerate(coords, 1):
                    fh.write(f"ATOM  {i % 100000:5d}  C   {resname:<3s} {chain}{(i % 9999) + 1:4d}    "
                             f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00           C\n")
                fh.write("END\n")
        # all in the design frame -> load together in VMD (see scripts/visualize_cylinder_fp.tcl)
        from dp3_native_cylinder import find_af3_cif
        cif = find_af3_cif(af3_dir)
        if cif is not None:
            read_gemmi(cif).write_pdb(str(out / "design.pdb"))            # full design model (ribbon)
        aca_al = (R @ aca.T).T + t
        ag_ca_nonepi = aca_al[[i for i in range(len(aca)) if i not in nepi_set]]
        write_points_pdb(out / "epitope_cas.pdb", epi_ca, "E", "EPI")      # epitope CAs (red)
        write_points_pdb(out / "antibody_aligned.pdb", ab_al, "Y", "AB")   # real antibody, point cloud
        write_points_pdb(out / "native_antigen.pdb", ag_ca_nonepi, "N", "NAG")  # native antigen (non-epi CAs)
        write_points_pdb(out / "native_antigen_heavy.pdb", ag_nonepi, "H", "NAH")  # the carve volume (heavy atoms)
        write_points_pdb(out / "flagged_cas.pdb", flagged_xyz, "Z", "FLG")       # all cylinder-flagged CAs
        write_points_pdb(out / "flagged_survive.pdb", flagged_xyz[survive], "S", "SUR")  # native-aware count
        write_points_pdb(out / "flagged_carved.pdb", flagged_xyz[carved], "C", "CRV")    # sit on native antigen
        (out / "flagged_scaffold_residues.txt").write_text(
            "0-based design residue indices of cylinder-flagged scaffold CAs:\n"
            + " ".join(map(str, flagged_res)) + "\n")
        bx, by, bz = base; nnx, nny, nnz = normal
        (out / "cylinder_frame.txt").write_text(
            f"base {bx} {by} {bz}\nnormal {nnx} {nny} {nnz}\nR {RADIUS}\nH {HEIGHT}\n")
        print(f"\nwrote {out}/  (design.pdb, epitope_cas.pdb, antibody_aligned.pdb, "
              f"native_antigen.pdb, flagged_survive.pdb, flagged_carved.pdb, flagged_cas.pdb, "
              f"cylinder_frame.txt)")
        print(f"VISUALIZE:  python scripts/plot_cylinder_fp_3d.py {out} --label '...'")


if __name__ == "__main__":
    main()
