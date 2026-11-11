"""Per-design AOMclash features in the epitope frame.

Faithful port of John Altin's `aomclash/aom_features.py` and the epitope frame from
`aomclash/aom_build.py` (AltinLab/epitope-universe @ feat/aomclash, 27ed71d), restructured so the
geometry/feature math takes coordinate arrays (I/O is a separate concern, so the same code runs on
our AF3 outputs rather than John's Tamarind campaign folders).

This file is built bottom-up and each piece is validated before the next is added:
  1. geometry core: `frame_from` (epitope frame) + `cylinder_counts` (the embedded cylinder), gated
     against episcaf_analysis.native_cylinder_core on a design with a recorded cylinder value.
  2. `design_feats` + `target_features`: the ~109 occupancy-map / shape / histogram features. The
     feature math is byte-for-byte John's; only the PDB reading and the contig->epitope-position
     resolution (specific to his Tamarind numbering) are lifted out, so the same code runs on our AF3
     outputs. Cross-checked exactly against John's `aom_features.design_feats` on a synthetic design
     (tests/test_aomclash_features.py), which also gates the cylinder sub-feature against
     native_cylinder_core.

Cylinder settings are Brandon's, from episcaf_analysis/native_cylinder_core.py + docs/CYLINDER_PARAMS.md.
"""
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree

# --- epitope-frame grid, matching aom_build.py (33^3 grid at 2 A, half-width 32) ---
HALF, STEP = 32.0, 2.0
AX = np.arange(-HALF, HALF + 1e-9, STEP)          # 33 samples per axis

# radial and height bins (A) for the scaffold-atom histograms in the epitope frame (aom_features.py)
RB = [0, 8, 16, 24, 32, 1e9]
ZB = [-1e9, -4, 0, 4, 8, 14, 22, 1e9]

# cylinder settings from bdneff/episcaf native_cylinder_core.py and docs/CYLINDER_PARAMS.md
CYL_R, CYL_H, CYL_OFF, CYL_CARVE = 16.0, 40.0, -4.0, 1.0


def frame_from(ca, ag_ca):
    """Epitope frame (port of aom_build.frame_from). Returns (origin, R, ratio); R rows are x, y, z.

    origin = epitope-CA centroid; z = outward normal (origin minus the centroid of antigen CA within
    25 A); x = epitope main axis in the plane perpendicular to z; y = z cross x. Returns None if the
    local antigen context has fewer than 20 CA atoms.
    """
    o = ca.mean(0)
    near = ag_ca[np.linalg.norm(ag_ca - o, axis=1) < 25.0]
    if len(near) < 20:
        return None
    n = o - near.mean(0)
    n /= np.linalg.norm(n)
    d = ca - o
    d = d - np.outer(d @ n, n)
    w, v = np.linalg.eigh(d.T @ d)
    x = v[:, 2]
    x = x - (x @ n) * n
    x /= np.linalg.norm(x)
    y = np.cross(n, x)
    return o, np.stack([x, y, n], 0), w[2] / max(w[1], 1e-6)


def kabsch(ref, mob):
    """Rotation and centres s.t. ref_point = R @ (mob_point - mob_centre) + ref_centre.

    Port of aom_features.kabsch.
    """
    rc, mc = ref.mean(0), mob.mean(0)
    U, S, Vt = np.linalg.svd((mob - mc).T @ (ref - rc))
    D = np.diag([1, 1, np.sign(np.linalg.det(Vt.T @ U.T))])
    return Vt.T @ D @ U.T, mc, rc


def hist2d(r, z):
    """2-D (radial x height) histogram of scaffold atoms in the epitope frame (aom_features.hist2d)."""
    h, _, _ = np.histogram2d(r, z, bins=[RB, ZB])
    return h.ravel()


def cylinder_counts(scaf_ca, epi_ca, ag_tree):
    """Native-aware cylinder count, port of aom_features.cylinder_counts.

    Counts scaffold CA atoms inside a cylinder standing on the epitope, minus those within CYL_CARVE
    of a native antigen atom. All coordinates share one rigid frame. `ag_tree` is a scipy cKDTree of
    the native antigen atoms in that same frame. Returns (native_aware_count, plain_count).

    This reproduces episcaf_analysis.native_cylinder_core.native_aware_scaffold_count; the two are
    cross-checked in the module smoke test.
    """
    cen = epi_ca.mean(0)
    _, _, Vt = np.linalg.svd(epi_ca - cen)
    nrm = Vt[-1]
    if nrm @ (np.vstack([epi_ca, scaf_ca]).mean(0) - cen) > 0:
        nrm = -nrm                                   # point away from the scaffold body
    base = cen + CYL_OFF * nrm
    v = scaf_ca - base
    proj = v @ nrm
    perp = np.linalg.norm(v - np.outer(proj, nrm), axis=1)
    ins = (proj >= 0) & (proj <= CYL_H) & (perp <= CYL_R)
    if not ins.any():
        return 0, 0
    d, _ = ag_tree.query(scaf_ca[ins])
    return int((d > CYL_CARVE).sum()), int(ins.sum())


def read_pdb_atoms(path):
    """Heavy atoms of a PDB as (name, chain, resid, x, y, z) tuples (port of aom_features.pdb_atoms).

    Keeps the altloc filter (`" A"`) and the hydrogen filter John uses, so a PDB read here matches the
    atom set his extractor sees. The AF3 CIF loader builds the same tuples via gemmi; this reader exists
    for PDB-based callers and the cross-check test.
    """
    out = []
    for l in open(path):
        if l.startswith("ATOM"):
            if l[16] not in " A":
                continue
            el = l[76:78].strip()
            name = l[12:16].strip()
            if (el == "H") if el else name.startswith("H"):
                continue
            out.append((name, l[21], int(l[22:26]), float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return out


def target_features(epi_ca, ag_ca, ag_atoms, n_isl):
    """Per-target geometry: the epitope frame, its KD-tree, and the `target_feats` block.

    Port of the geometry portion of aom_features.one_target (the campaign I/O is lifted out). Inputs are
    native-antigen coordinates (epitope chains only, never antibody):
      epi_ca   -- (E,3) epitope CA
      ag_ca    -- (M,3) antigen CA of the epitope chains (for the outward-normal frame)
      ag_atoms -- (K,3) antigen heavy atoms of the epitope chains (the clash reference)
      n_isl    -- number of epitope islands
    Returns (o, Rf, ag_tree, epi_L, target_feats) or None if the frame is undefined (<20 local CA).
    """
    fr = frame_from(epi_ca, ag_ca)
    if fr is None:
        return None
    o, Rf, ratio = fr
    ag_L = (ag_atoms - o) @ Rf.T
    ag_tree = cKDTree(ag_L)
    epi_L = (epi_ca - o) @ Rf.T
    out_ag = ag_L[(ag_L[:, 2] > 0) & (np.linalg.norm(ag_L, axis=1) < 25)]
    target_feats = dict(frame_ratio=ratio, n_epi=len(epi_ca), n_isl=n_isl,
                        ag_out25=len(out_ag), ag_out25_hi=int((out_ag[:, 2] > 5).sum()),
                        ag_near15=int((np.linalg.norm(ag_L, axis=1) < 15).sum()))
    return o, Rf, ag_tree, epi_L, target_feats


def design_feats(atoms, epi_pos, segs, ref, o, Rf, ag_tree, epi_L, P1, P2):
    """Features of one design in the epitope frame (faithful port of aom_features.design_feats).

    The feature math is byte-for-byte John's. Lifted out of his version: the PDB read (`atoms` passed
    in) and the contig->epitope-position derivation (`epi_pos` passed in). Inputs:
      atoms   -- design heavy atoms as (name, chain, resid, x, y, z); chain "A" is the design chain.
      epi_pos -- ordered design residue ids that are the epitope, in the SAME order as `ref`'s residues
                 (island order, ascending within island) so the Kabsch pairing is correct.
      segs    -- scaffold segment lengths (for n_seg / seg_max / seg_sum).
      ref     -- native epitope backbone atoms (N,CA,C,O per residue, same order as epi_pos).
      o, Rf   -- epitope frame from `target_features`.
      ag_tree -- KD-tree of native antigen atoms in the epitope frame.
      epi_L   -- epitope CA in the epitope frame (for the cylinder).
      P1, P2  -- occupancy map (33^3) smoothed at 1.0 and 2.5 voxels (antibody-geometry only).
    Returns a dict of ~105 features, or None if the design lacks the epitope backbone / scaffold atoms.
    """
    epi_pos = list(epi_pos)
    epi_set = set(epi_pos)
    bb = {(s, n): np.array(xyz) for n, c, s, *xyz in atoms if c == "A" and n in ("N", "CA", "C", "O")}
    try:
        mob = np.array([bb[(p, n)] for p in epi_pos for n in ("N", "CA", "C", "O")])
    except KeyError:
        return None
    Rk, mc, rc = kabsch(ref, mob)
    sc = [(s, n, np.array(xyz)) for n, c, s, *xyz in atoms if c == "A" and s not in epi_set]
    if not sc:
        return None
    res = np.array([s for s, _, _ in sc])
    isca = np.array([n == "CA" for _, n, _ in sc])
    X = np.array([x for _, _, x in sc])
    L = ((X - mc) @ Rk.T + rc - o) @ Rf.T                  # scaffold atoms in the epitope frame
    gi = ((L + HALF) / STEP).T
    p1 = map_coordinates(P1, gi, order=1, mode="constant", cval=0.0)   # P(antibody), smoothed 2 A
    p2 = map_coordinates(P2, gi, order=1, mode="constant", cval=0.0)   # smoothed 5 A
    dag, _ = ag_tree.query(L)
    dist = np.linalg.norm(L, axis=1)
    rr = np.hypot(L[:, 0], L[:, 1])
    z = L[:, 2]
    free = dag > 4.0                                        # not inside antigen space
    ures, inv = np.unique(res, return_inverse=True)
    lq = np.zeros(len(ures)); np.add.at(lq, inv, np.log1p(-np.minimum(p1, 0.999)))
    lqf = np.zeros(len(ures)); np.add.at(lqf, inv, np.log1p(-np.minimum(p1 * free, 0.999)))
    pmax = np.zeros(len(ures)); np.maximum.at(pmax, inv, p1)
    cax = L[isca]
    cyl_aware, cyl_plain = cylinder_counts(cax, epi_L, ag_tree)
    feat = dict(
        cyl_aware=cyl_aware, cyl_plain=cyl_plain,
        apr=len(X) / len(ures), n_scaf=len(ures), n_seg=len(segs), seg_max=max(segs), seg_sum=sum(segs),
        s_p1=p1.sum(), s_p2=p2.sum(), s_p1_free=(p1 * free).sum(), s_p1_near=(p1 * ~free).sum(),
        n_p05=(p1 > .05).sum(), n_p15=(p1 > .15).sum(), n_p30=(p1 > .30).sum(), n_p50=(p1 > .5).sum(),
        n_p15_free=((p1 > .15) & free).sum(), n_p30_free=((p1 > .30) & free).sum(),
        r_sum=(1 - np.exp(lq)).sum(), r_sum_free=(1 - np.exp(lqf)).sum(),
        r_n30=(pmax > .3).sum(), r_n50=(pmax > .5).sum(), r_n70=(pmax > .7).sum(),
        b_near20=(dist < 20).sum(), b_out25=((dist < 25) & (z > 0)).sum(), b_out15=((dist < 15) & (z > 0)).sum(),
        b_out35_hi=((dist < 35) & (z > 4)).sum(),
        ag3=(dag < 3).sum(), ag4=(dag < 4).sum(), ag6=(dag < 6).sum(),
        ag_free_out25=(((dist < 25) & (z > 0)) & free).sum(),
        z_mean=z.mean(), z_p90=np.percentile(z, 90), z_max=z.max(), r_mean=rr.mean(),
        ca_rg=np.sqrt(((cax - cax.mean(0)) ** 2).sum(1).mean()))
    for j, v in enumerate(hist2d(rr, z)):
        feat[f"h{j}"] = v
    for j, v in enumerate(hist2d(rr[free], z[free])):
        feat[f"hf{j}"] = v
    return feat
