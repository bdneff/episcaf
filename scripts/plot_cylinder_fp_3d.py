#!/usr/bin/env python3
"""Render a cylinder_fp_probe.py output directory into two clean panels, in the cylinder frame.

Transforms every atom into the cylinder's own coordinates -- axial = distance up the approach
normal, (r1, r2) = the two in-plane directions -- so two views tell the whole story:
  SIDE  (r1 vs axial): the antibody sits high in the cylinder, the flagged scaffold at the base.
  TOP   (r1 vs r2, looking down the normal): whether the flagged scaffold spreads LATERALLY in
        the epitope plane (the 8pww false-positive) or stays tight on-axis.

Usage:
  python scripts/plot_cylinder_fp_3d.py results/cylinder_fp/DP2_0804 --label "8pww  DP2_0804 (cyl 10)"
"""
from __future__ import annotations
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams.update({"font.size": 14, "axes.titlesize": 16, "axes.labelsize": 15,
                     "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12})


def read_pdb(path: Path, ca_only=False):
    pts = []
    for ln in Path(path).read_text().splitlines():
        if not ln.startswith(("ATOM", "HETATM")):
            continue
        if ca_only and ln[12:16].strip() != "CA":
            continue
        pts.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return np.asarray(pts, float) if pts else np.empty((0, 3))


def frame(d: Path):
    base = normal = None; R = H = None
    for ln in (d / "cylinder_frame.txt").read_text().splitlines():
        t = ln.split()
        if t[0] == "base":   base = np.array(list(map(float, t[1:4])))
        if t[0] == "normal": normal = np.array(list(map(float, t[1:4])))
        if t[0] == "R":      R = float(t[1])
        if t[0] == "H":      H = float(t[1])
    normal = normal / np.linalg.norm(normal)
    ax = np.array([1.0, 0, 0]) if abs(normal[0]) < 0.9 else np.array([0, 1.0, 0])
    p1 = np.cross(normal, ax); p1 /= np.linalg.norm(p1)
    p2 = np.cross(normal, p1)
    return base, normal, p1, p2, R, H


def to_frame(pts, base, normal, p1, p2):
    if len(pts) == 0:
        return np.empty((0,)), np.empty((0,)), np.empty((0,))
    v = pts - base
    return v @ p1, v @ p2, v @ normal   # r1, r2, axial


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("indir")
    ap.add_argument("--label", default="")
    ap.add_argument("--out", default="")
    args = ap.parse_args()
    d = Path(args.indir)
    base, normal, p1, p2, R, H = frame(d)

    def rp(name, ca=False):
        p = d / name
        return read_pdb(p, ca_only=ca) if p.exists() else np.empty((0, 3))

    CARVE = 1.0  # native-aware carve distance: a scaffold CA within CARVE of a heavy atom is carved

    sets = {
        "scaffold":  (rp("design.pdb", ca=True), dict(s=8, c="0.82", zorder=1)),
        "antibody":  (rp("antibody_aligned.pdb"), dict(s=3, c="tab:blue", alpha=0.10, zorder=2)),
        "epitope":   (rp("epitope_cas.pdb"), dict(s=55, c="crimson", edgecolors="k", linewidths=0.4, zorder=5)),
    }
    if (d / "flagged_survive.pdb").exists():
        sets["flagged, carved"] = (rp("flagged_carved.pdb"),
                                   dict(s=60, c="0.35", marker="x", linewidths=1.8, zorder=6))
        sets["flagged, counts"] = (rp("flagged_survive.pdb"),
                                   dict(s=90, c="orange", edgecolors="k", linewidths=0.6, zorder=7))
    else:
        sets["flagged"] = (rp("flagged_cas.pdb"), dict(s=70, c="orange", edgecolors="k", linewidths=0.5, zorder=7))
    fr = {k: to_frame(v[0], base, normal, p1, p2) for k, v in sets.items()}

    # the carve VOLUME: native-antigen heavy atoms, each a CARVE-radius disk (union = the carved zone)
    heavy = rp("native_antigen_heavy.pdb")
    if len(heavy) == 0:
        heavy = rp("native_antigen.pdb")   # fallback: CAs (coarser)
    hr1, hr2, hax = to_frame(heavy, base, normal, p1, p2)
    # slab the heavy atoms to the axial band the flagged scaffold lives in (that's what carves them)
    fax = np.concatenate([fr[k][2] for k in sets if k.startswith("flagged")]) if any(
        k.startswith("flagged") for k in sets) else np.array([0.0, 6.0])
    lo, hi = (fax.min() - CARVE - 1, fax.max() + CARVE + 1) if len(fax) else (-2, 8)
    slab = (hax >= lo) & (hax <= hi)

    def carve_outline(ax, x, y, xr, yr, n=220):
        """Project the exclusion SPHERES (radius CARVE around each heavy atom) onto this view:
        a sphere projects to a CARVE-radius circle at any depth, so the union of those circles
        is the carve region in this plane -- circular cutouts in both side and top-down views."""
        from scipy.spatial import cKDTree
        if len(x) == 0:
            return
        X, Y = np.meshgrid(np.linspace(*xr, n), np.linspace(*yr, n))
        dmin, _ = cKDTree(np.c_[x, y]).query(np.c_[X.ravel(), Y.ravel()], k=1)
        Z = (dmin.reshape(X.shape) <= CARVE).astype(float)
        ax.contourf(X, Y, Z, levels=[0.5, 2], colors=["mediumseagreen"], alpha=0.18, zorder=1)
        ax.contour(X, Y, Z, levels=[0.5], colors=["seagreen"], linewidths=1.6, zorder=3)

    fig, (axS, axT) = plt.subplots(1, 2, figsize=(12, 6))
    # SIDE view: r1 vs axial -- carve region (antigen sphere-projection) reaching up to the flags
    carve_outline(axS, hr1, hax, (-R - 4, R + 4), (hax.min() - 2, 12))
    for k, (_, st) in sets.items():
        r1, r2, ax_ = fr[k]
        axS.scatter(r1, ax_, label=k, **st)
    axS.add_patch(plt.Rectangle((-R, 0), 2 * R, H, fill=False, ec="tab:cyan", lw=1.5))
    axS.axhline(4.0, ls="--", c="0.4", lw=1)   # epitope plane (offset -4)
    axS.set_xlabel("in-plane distance from axis (Å)")
    axS.set_ylabel("height up the approach normal (Å)")
    axS.set_title("side view"); axS.set_xlim(-R - 4, R + 4)
    # TOP-DOWN: same carve region, looking down the normal (slabbed to the flags' axial band)
    carve_outline(axT, hr1[slab], hr2[slab], (-R - 4, R + 4), (-R - 4, R + 4))
    for k, (_, st) in sets.items():
        r1, r2, ax_ = fr[k]
        axT.scatter(r1, r2, **st)
    th = np.linspace(0, 2 * np.pi, 100)
    axT.plot(R * np.cos(th), R * np.sin(th), c="tab:cyan", lw=1.5)
    axT.set_xlabel("in-plane x (Å)"); axT.set_ylabel("in-plane y (Å)")
    axT.set_title("top-down: green = carve volume (native antigen ±1Å)")
    axT.set_aspect("equal"); axT.set_xlim(-R - 4, R + 4); axT.set_ylim(-R - 4, R + 4)
    axS.legend(loc="upper right", framealpha=0.9, markerscale=1.2, fontsize=11)

    fig.suptitle(args.label or d.name, fontsize=17)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    out = Path(args.out) if args.out else d / "fp_3d.png"
    fig.savefig(out, dpi=150)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
