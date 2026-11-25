"""Per-design viewer assets: align an AF2-folded design onto the native epitope frame, then emit a
small scaffold-only PDB per design plus ONE shared native PDB (antigen + antibody) for the whole run,
and a compact config JSON the results page renders with 3Dmol.

The native complex is fixed — every design is Kabsch-moved onto the native epitope, the native never
moves — so `viewer/native.pdb` is written once and shared by all designs; each `viewer/<id>.pdb` holds
just the aligned scaffold (chain S). This keeps designs.zip small and lets the viewer load the native
once and swap scaffolds per row.

This is a *display* companion to the scored metrics, not a re-score. Clashes are geometric (scaffold
heavy atoms near the antibody) and use a detected antibody even when the campaign scored with none (the
growPatch avenue passes `antibody_chains=[]` but the uploaded complex may still contain the Fab). The
alignment matches `campaign.evaluate_design`'s epitope RMSD (both Kabsch-superpose the design epitope
backbone onto the native), so `epitopeRMSD` reproduces the CSV's `epitope_rmsd` to rounding — verified
on the 2026-11-25 benchmark (3HFM, Y:14-21;Y:74-100): 0.585 / 0.600 / 0.844, antibody detected L+H.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from episcaf_service.rfd4 import topology

AA3 = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
       "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
       "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}
_BB = ("N", "CA", "C", "O")
SCAFFOLD_CHAIN = "S"


def parse_pdb(path) -> list[dict]:
    atoms = []
    for ln in Path(path).read_text().splitlines():
        if ln[:6] not in ("ATOM  ", "HETATM"):
            continue
        try:
            atoms.append(dict(name=ln[12:16].strip(), resn=ln[17:20].strip(), chain=ln[21],
                              resi=int(ln[22:26]), x=float(ln[30:38]), y=float(ln[38:46]),
                              z=float(ln[46:54]), elem=(ln[76:78].strip() or ln[12:16].strip()[:1])))
        except ValueError:
            continue
    return atoms


def _residues(atoms, chain) -> dict:
    res = {}
    for a in atoms:
        if a["chain"] != chain or a["resn"] not in AA3:
            continue
        r = res.setdefault(a["resi"], {"aa": AA3[a["resn"]], "bb": {}})
        if a["name"] in _BB:
            r["bb"][a["name"]] = np.array([a["x"], a["y"], a["z"]])
    return res


def _bb_stack(res, resis):
    pts, kept = [], []
    for i in resis:
        bb = res[i]["bb"]
        if all(k in bb for k in _BB):
            pts.extend(bb[k] for k in _BB)
            kept.append(i)
    return (np.array(pts) if pts else np.empty((0, 3))), kept


def kabsch(mob: np.ndarray, ref: np.ndarray):
    """Rotation+translation mapping mob onto ref (both Nx3)."""
    cm, cr = mob.mean(0), ref.mean(0)
    H = (mob - cm).T @ (ref - cr)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    R = Vt.T @ np.diag([1, 1, d]) @ U.T
    return R, cr - R @ cm


def _heavy_xyz(atoms, chains) -> np.ndarray:
    chains = set(chains)
    pts = [(a["x"], a["y"], a["z"]) for a in atoms
           if a["chain"] in chains and a["resn"] in AA3 and a["elem"].upper() != "H"]
    return np.array(pts, float) if pts else np.empty((0, 3))


def detect_antibody_chains(atoms, antigen_chains, epitope_heavy: np.ndarray,
                           cutoff: float = 4.5, min_contacts: int = 10) -> list[str]:
    """Contact-based antibody detection (John §1.4 fallback): a non-antigen protein chain is the
    antibody if it has at least `min_contacts` heavy-atom pairs within `cutoff` of the epitope. Used
    for display when the campaign scored without an antibody."""
    if epitope_heavy.size == 0:
        return []
    out = []
    for ch in [c for c in dict.fromkeys(a["chain"] for a in atoms) if c not in set(antigen_chains)]:
        hv = _heavy_xyz(atoms, [ch])
        if hv.size == 0:
            continue
        d2 = ((hv[:, None, :] - epitope_heavy[None, :, :]) ** 2).sum(-1)
        if int((d2 < cutoff * cutoff).sum()) >= min_contacts:
            out.append(ch)
    return out


def _fmt_atom(serial, name, resn, chain, resi, xyz, elem) -> str:
    nm = name if len(name) == 4 else " " + name
    return ("ATOM  %5d %-4s %3s %s%4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s\n"
            % (serial, nm.ljust(4), resn, chain, resi, xyz[0], xyz[1], xyz[2], elem[:2].rjust(2)))


def _write_pdb(atoms, path, chain_of=None):
    with open(path, "w") as fh:
        for s, a in enumerate(atoms, 1):
            ch = chain_of(a) if chain_of else a["chain"]
            fh.write(_fmt_atom(s, a["name"], a["resn"], ch, a["resi"],
                               (a["x"], a["y"], a["z"]), a["elem"]))
        fh.write("END\n")


def align_design(native, design, epitope, antibody_chains=None):
    """Move `design` atoms onto the native epitope frame. Returns (cfg, antibody_chains). `native`
    and `design` are atom lists from parse_pdb; `epitope` is an islands string or Island list.
    Raises ValueError if an island's native sequence is not found in the design."""
    islands = topology.parse_islands(epitope) if isinstance(epitope, str) else list(epitope)
    antigen_chains = list(dict.fromkeys(i.chain for i in islands))
    nat_res = {ch: _residues(native, ch) for ch in antigen_chains}
    model_chain = max({a["chain"] for a in design if a["resn"] in AA3},
                      key=lambda c: sum(1 for a in design if a["chain"] == c))
    des_res = _residues(design, model_chain)
    des_order = sorted(des_res)
    des_seq = "".join(des_res[i]["aa"] for i in des_order)

    mob, ref, epi_native, epi_design = [], [], [], []
    for isl in islands:
        res = nat_res[isl.chain]
        nres = [i for i in range(isl.start, isl.end + 1) if i in res]
        isl_seq = "".join(res[i]["aa"] for i in nres)
        pos = des_seq.find(isl_seq)
        if pos < 0 or not isl_seq:
            raise ValueError(f"island {isl.chain}:{isl.start}-{isl.end} ({isl_seq!r}) not in design")
        dres = des_order[pos:pos + len(nres)]
        epi_design += dres
        epi_native += [[isl.chain, i] for i in nres]
        dm, _ = _bb_stack(des_res, dres)
        nm, _ = _bb_stack(res, nres)
        k = min(len(dm), len(nm))
        mob.append(dm[:k]); ref.append(nm[:k])
    mob, ref = np.vstack(mob), np.vstack(ref)
    R, t = kabsch(mob, ref)
    rmsd = float(np.sqrt((((R @ mob.T).T + t - ref) ** 2).sum(1).mean()))
    for a in design:
        a["x"], a["y"], a["z"] = (R @ np.array([a["x"], a["y"], a["z"]]) + t).tolist()

    epi_native_set = {(c, i) for c, i in epi_native}
    epi_heavy = np.array([[a["x"], a["y"], a["z"]] for a in native
                          if a["chain"] in antigen_chains and (a["chain"], a["resi"]) in epi_native_set
                          and a["elem"].upper() != "H"], float)
    ab_chains = list(antibody_chains) if antibody_chains else detect_antibody_chains(
        native, antigen_chains, epi_heavy)
    ab_heavy = _heavy_xyz(native, ab_chains)

    epi_set = set(epi_design)
    clash = set()
    if ab_heavy.size:
        for a in design:
            if a["resn"] not in AA3 or a["elem"].upper() == "H" or a["resi"] in epi_set:
                continue
            p = np.array([a["x"], a["y"], a["z"]])
            if ((ab_heavy - p) ** 2).sum(1).min() < 16.0:
                clash.add(a["resi"])

    epi_ca = np.array([[a["x"], a["y"], a["z"]] for a in design
                       if a["resi"] in epi_set and a["name"] == "CA"])
    epi_c = epi_ca.mean(0) if len(epi_ca) else mob.mean(0)
    if ab_heavy.size:
        dirn = ab_heavy.mean(0) - epi_c
        dirn = dirn / (np.linalg.norm(dirn) or 1.0)
        aom_center = (epi_c + dirn * 12.0).tolist()
    else:
        aom_center = epi_c.tolist()

    cfg = dict(scaffoldChain=SCAFFOLD_CHAIN, antigenChains=antigen_chains, antibodyChains=ab_chains,
               epitopeScaffold=sorted(epi_design), epitopeNative=epi_native, clash=sorted(clash),
               aomCenter=[round(v, 2) for v in aom_center], aomRadius=13.0,
               epitopeRMSD=round(rmsd, 3), sequence=des_seq, length=len(des_seq),
               nativePdb="native.pdb")
    return cfg, ab_chains


def write_native(native, antigen_chains, antibody_chains, path):
    """Write the shared native complex (antigen + antibody chains, original ids) once per run."""
    keep = set(antigen_chains) | set(antibody_chains)
    _write_pdb([a for a in native if a["chain"] in keep and a["resn"] in AA3], path)


def write_viewer_files(out_dir, design_id, native_path, model_pdb, epitope, antibody_chains=None):
    """Emit viewer/<id>.pdb (aligned scaffold) + viewer/<id>.json, and viewer/native.pdb once.
    Returns the config, or None if the design could not be aligned (missing epitope sequence)."""
    vdir = Path(out_dir) / "viewer"
    vdir.mkdir(parents=True, exist_ok=True)
    native = parse_pdb(native_path)
    design = parse_pdb(model_pdb)
    try:
        cfg, ab = align_design(native, design, epitope, antibody_chains or None)
    except (ValueError, KeyError):
        return None
    _write_pdb([a for a in design if a["resn"] in AA3], vdir / f"{design_id}.pdb",
               chain_of=lambda a: SCAFFOLD_CHAIN)
    (vdir / f"{design_id}.json").write_text(json.dumps(cfg, indent=2))
    native_out = vdir / "native.pdb"
    if not native_out.is_file():
        write_native(native, cfg["antigenChains"], ab, native_out)
    return cfg


if __name__ == "__main__":
    import sys
    out = sys.argv[5] if len(sys.argv) > 5 else "."
    cfg = write_viewer_files(out, "design", sys.argv[1], sys.argv[2], sys.argv[3],
                             antibody_chains=(sys.argv[4].split(",") if len(sys.argv) > 4 and sys.argv[4] else None))
    print(f"epitope RMSD {cfg['epitopeRMSD']} A | antibody {cfg['antibodyChains']} | "
          f"{len(cfg['clash'])} clash residues | aom {cfg['aomCenter']}" if cfg else "could not align")
