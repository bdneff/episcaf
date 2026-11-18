#!/usr/bin/env python3
"""Render a 2D slice of the REAL antibody-occupancy map for the AOMclash doc figure.

Loads the stored map (episcaf_analysis/aomclash/aom_map.pkl), smooths it as the scorer does
(P1), finds the outward axis empirically (the direction the occupancy is displaced toward),
and plots the plane that contains the outward axis as a heatmap of P(antibody) in the epitope
frame, in angstroms, with a color scale, the epitope origin, the outward direction, the 4 A
native-antigen free boundary, and two example scaffold atoms (one counted, one excluded).

    python scripts/render_aom_slice.py            # -> docs/img/fig/aom_slice.png

Not schematic: the field is the actual stored occupancy map. The native-antigen surface and the
two atom markers are illustrative overlays (the map does not store a specific antigen); they are
labeled as such in the figure.
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrow

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from episcaf_analysis.aomclash import aom_map as AM

STEP = 2.0  # A / voxel (aom_map.py)
OUT = os.path.join(ROOT, "docs/img/fig/aom_slice.png")
MAP = os.path.join(ROOT, "episcaf_analysis/aomclash/aom_map.pkl")


def main():
    art = AM.load_map(MAP)
    field, n_used = AM.mean_field(art)
    P1, _ = AM.smoothed_maps(field)
    n = P1.shape[0]
    c = n // 2
    coord = (np.arange(n) - c) * STEP  # A, origin at center

    # outward axis = axis whose occupancy-weighted centroid is most displaced from center
    disp = []
    for ax in range(3):
        w = P1.sum(axis=tuple(i for i in range(3) if i != ax))
        com = (np.arange(n) * w).sum() / w.sum()
        disp.append(abs(com - c))
    zaxis = int(np.argmax(disp))
    inax = [a for a in range(3) if a != zaxis]
    xax = inax[0]
    # orient outward = positive
    w = P1.sum(axis=tuple(i for i in range(3) if i != zaxis))
    com = (np.arange(n) * w).sum() / w.sum()
    flip = com < c

    # slice at the center of the third axis
    sl = [slice(None)] * 3
    sl[inax[1]] = c
    plane = P1[tuple(sl)]              # shape (n, n) over (xax, zaxis) in some order
    # arrange as [outward(rows), in-plane(cols)]
    plane = np.moveaxis(P1, (zaxis, xax), (0, 1))[:, :, c]
    if flip:
        plane = plane[::-1, :]
    peak = plane.max()
    peak_out = coord[np.unravel_index(plane.argmax(), plane.shape)[0]]

    fig, axp = plt.subplots(figsize=(5.2, 4.4), dpi=150)
    extent = [coord[0], coord[-1], coord[0], coord[-1]]  # in-plane x, outward y
    im = axp.imshow(plane, origin="lower", extent=extent, cmap="Blues",
                    vmin=0, vmax=peak, aspect="equal")
    cb = fig.colorbar(im, ax=axp, fraction=0.046, pad=0.04)
    cb.set_label("P(antibody occupies this point)", fontsize=9)
    cb.ax.tick_params(labelsize=8)

    # epitope origin + outward arrow
    axp.plot(0, 0, "o", color="#c0392b", ms=8, zorder=5)
    axp.annotate("epitope origin", (0, 0), (0, -13), fontsize=8, color="#c0392b",
                 ha="center", va="top", annotation_clip=False)
    axp.add_patch(FancyArrow(0, 2, 0, 18, width=0.6, head_width=3, head_length=4,
                             color="#c0392b", alpha=0.8, length_includes_head=True))
    axp.annotate("outward\n(antibody side)", (3, 24), fontsize=8, color="#c0392b",
                 va="center", ha="left", annotation_clip=False)

    # 4 A free boundary (native-antigen proximity rule) -- illustrative overlay
    axp.axhline(4.0, ls="--", color="#555555", lw=1)
    axp.annotate("4 Å free boundary", (coord[0] + 2, 5.6),
                 fontsize=7.5, color="#555555", annotation_clip=False)

    # two example scaffold atoms: one counted (free, in occupancy), one excluded (near antigen)
    axp.plot(-6, 14, "o", color="#1b4a86", ms=9, mec="white", zorder=6)
    axp.annotate("counted\n(free atom)", (-6, 14), (-26, 23),
                 fontsize=7.5, color="#1b4a86", ha="center", va="center",
                 annotation_clip=False,
                 arrowprops=dict(arrowstyle="->", color="#1b4a86", lw=0.8))
    axp.plot(7, 1.5, "o", color="#9aa3ad", ms=9, mec="white", zorder=6)
    axp.annotate("excluded\n(< 4 Å to antigen)", (7, 1.5), (19, -12),
                 fontsize=7.5, color="#5c636b", ha="center", va="center",
                 annotation_clip=False,
                 arrowprops=dict(arrowstyle="->", color="#9aa3ad", lw=0.8))

    axp.set_xlabel("in-plane position (Å)", fontsize=9)
    axp.set_ylabel("outward from epitope (Å)", fontsize=9)
    axp.tick_params(labelsize=8)
    axp.set_xlim(coord[0], coord[-1]); axp.set_ylim(coord[0], coord[-1])
    fig.tight_layout()
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fig.savefig(OUT, transparent=True, bbox_inches="tight")
    print("wrote %s  (map over %d complexes, peak P=%.3f at ~%.0f A out)" % (OUT, n_used, peak, peak_out))


if __name__ == "__main__":
    main()
