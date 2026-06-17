#!/usr/bin/env python3
"""
native_cylinder_core.py

Native-antigen-aware cylinder clash metric.

Drop these functions into the existing clash code, right after the native complex
is aligned onto the AF3 epitope (i.e. once you have the native antigen heavy-atom
coordinates in the AF3 frame). They reuse the same epitope-plane / cylinder geometry
as add_cylinder_metric so the numbers stay consistent.

Two outputs per design:
  native_in_cylinder      -- TEST: how many native-antigen (non-epitope) CAs fall
                             inside the cylinder. Large => the cylinder is flagging
                             volume nature filled and the real antibody tolerated.
  cylinder_native_aware   -- SOLUTION: scaffold CAs inside the cylinder, EXCLUDING
                             those sitting within `exclude_dist` of a native-antigen
                             heavy atom (i.e. in nature-allowed space).

The native antigen is always available (the epitope is taken from it), so both
generalize to no-antibody targets.
"""
import numpy as np
from scipy.spatial import cKDTree

RADIUS, HEIGHT, OFFSET = 16.0, 40.0, -4.0   # locked cylinder params


def cylinder_frame(epi_ca, all_ca, radius=RADIUS, height=HEIGHT, offset=OFFSET):
    """Epitope-plane cylinder, identical to add_cylinder_metric."""
    centroid = epi_ca.mean(axis=0)
    _, _, Vt = np.linalg.svd(epi_ca - centroid)
    normal = Vt[-1]
    if np.dot(normal, all_ca.mean(axis=0) - centroid) > 0:
        normal = -normal
    base = centroid + offset * normal
    return base, normal


def inside_cylinder(points, base, normal, radius=RADIUS, height=HEIGHT):
    v = points - base
    proj = v @ normal
    perp = v - np.outer(proj, normal)
    dist = np.linalg.norm(perp, axis=1)
    return (proj >= 0.0) & (proj <= height) & (dist <= radius)


def count_native_in_cylinder(native_ca, native_is_epitope, base, normal,
                             radius=RADIUS, height=HEIGHT):
    """TEST: native non-epitope CAs inside the cylinder."""
    ins = inside_cylinder(native_ca, base, normal, radius, height)
    return int(np.sum(ins & ~native_is_epitope))


def native_aware_scaffold_count(scaf_ca, base, normal, native_heavy_xyz,
                                exclude_dist=4.0, radius=RADIUS, height=HEIGHT):
    """SOLUTION: scaffold CAs inside the cylinder but NOT in native-occupied space.
    Returns (native_aware_count, plain_count) so you can see what the carve-out removed.
    native_heavy_xyz: native-antigen heavy-atom coords in the AF3 frame (epitope
    atoms may be included or not; they sit at the epitope so they rarely matter)."""
    ins = inside_cylinder(scaf_ca, base, normal, radius, height)
    plain = int(ins.sum())
    if plain == 0:
        return 0, 0
    inside_pts = scaf_ca[ins]
    tree = cKDTree(native_heavy_xyz)
    d, _ = tree.query(inside_pts, k=1)
    keep = d > exclude_dist               # not sitting in native volume
    return int(keep.sum()), plain


# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    # ---- unit tests on synthetic geometry --------------------------------- #
    rng = np.random.default_rng(0)

    # epitope plane in z=0, scaffold body at -z so normal points +z
    epi = np.array([[0,0,0],[3,0,0],[0,3,0],[3,3,0]], float)
    body = np.array([[1,1,-30]], float)
    base, normal = cylinder_frame(epi, np.vstack([epi, body]))
    assert normal[2] > 0, "normal should point away from scaffold body (+z)"

    # inside/outside checks (base is at centroid + (-4)*normal => z=-4)
    pts = np.array([
        [1,1, 10],    # on-axis, within height -> inside
        [1,1, 50],    # beyond height -> outside
        [1,1,-10],    # behind base -> outside
        [1+20,1, 10], # outside radius -> outside
    ], float)
    ins = inside_cylinder(pts, base, normal)
    assert list(ins) == [True, False, False, False], list(ins)

    # TEST counter: 2 native non-epitope CAs inside, 1 epitope inside (excluded), 1 outside
    nat = np.array([[1,1,8],[0,0,12],[1,1,5],[40,0,8]], float)
    nat_epi = np.array([False, False, True, False])
    assert count_native_in_cylinder(nat, nat_epi, base, normal) == 2

    # SOLUTION carve-out: 3 scaffold CAs inside; one sits on a native atom -> removed
    scaf = np.array([[1,1,8],[2,1,12],[0,1,20]], float)
    native_heavy = np.array([[1,1,8.0]], float)   # within 4A of first scaffold CA only
    aware, plain = native_aware_scaffold_count(scaf, base, normal, native_heavy,
                                               exclude_dist=4.0)
    assert plain == 3 and aware == 2, (plain, aware)

    # Kabsch sanity (alignment you'd use to bring native epitope onto AF3 epitope)
    def kabsch(P, Q):
        Pc, Qc = P - P.mean(0), Q - Q.mean(0)
        U, _, Vt = np.linalg.svd(Pc.T @ Qc)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        D = np.diag([1, 1, d])
        R = Vt.T @ D @ U.T
        return R, Q.mean(0) - R @ P.mean(0)
    P = rng.normal(size=(8, 3))
    Rtrue = np.linalg.qr(rng.normal(size=(3, 3)))[0]
    if np.linalg.det(Rtrue) < 0: Rtrue[:, 0] *= -1
    ttrue = np.array([5.0, -2.0, 1.0])
    Q = (Rtrue @ P.T).T + ttrue
    R, t = kabsch(P, Q)
    assert np.allclose((R @ P.T).T + t, Q, atol=1e-9)

    print("all geometry unit tests passed")
