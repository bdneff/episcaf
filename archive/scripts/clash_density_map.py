#!/usr/bin/env python3
"""
clash_density_map.py

Build a 3D spatial density map of antibody clash coordinates for a single epitope.

For each design, the clash detection already aligned the AF3 structure onto the
true complex epitope, so clashing antibody residues are in the true complex
coordinate frame. We can therefore directly use the true PDB atom coordinates
for the clashing residue indices.

The density map is built in an epitope-centered, PCA-aligned coordinate frame
so that it is consistent across all designs and interpretable relative to the
epitope surface.

Usage:
    python scripts/clash_density_map.py \
        --metrics_csv  runs/run_rfd3_mpnn/04_filter/metrics_partial.csv \
        --true_pdb     /tgen_labs/altin/.../cleaned/4xwo_5P.pdb \
        --epitope_id   4xwo_5P \
        --dp2_parquet  datasets/dp2.parquet \
        --out_dir      runs/run_rfd3_mpnn/clash_density/4xwo_5P
"""

import argparse
import ast
from pathlib import Path

import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter
from scipy.stats import gaussian_kde
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_resindices(x) -> list:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return []
    if isinstance(x, list):
        return [int(i) for i in x]
    s = str(x).strip()
    if not s or s == "[]" or s == "nan":
        return []
    try:
        parsed = ast.literal_eval(s)
        return [int(i) for i in parsed]
    except Exception:
        return []


def pca_frame(coords: np.ndarray):
    """
    Return (center, axes) where axes is a 3x3 matrix whose rows are the
    PCA principal axes (sorted by descending variance).
    Used to define a consistent coordinate frame for the epitope.
    """
    center = coords.mean(axis=0)
    centered = coords - center
    cov = centered.T @ centered
    eigenvalues, eigenvectors = np.linalg.eigh(cov)
    # sort descending
    order = np.argsort(eigenvalues)[::-1]
    axes = eigenvectors[:, order].T  # rows are axes
    return center, axes


def to_local_frame(coords: np.ndarray, center: np.ndarray, axes: np.ndarray) -> np.ndarray:
    """Project coords into local PCA frame centered at epitope centroid."""
    return (coords - center) @ axes.T


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--metrics_csv",  required=True)
    parser.add_argument("--true_pdb",     required=True)
    parser.add_argument("--epitope_id",   required=True, help="e.g. 4xwo_5P")
    parser.add_argument("--dp2_parquet",  required=True)
    parser.add_argument("--out_dir",      required=True)
    parser.add_argument("--voxel_size",   type=float, default=1.0,
                        help="Voxel size in Angstroms (default: 1.0)")
    parser.add_argument("--kde_bw",       type=float, default=2.0,
                        help="KDE bandwidth in Angstroms (default: 2.0)")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # --- Load data ---
    print(f"Loading metrics CSV ...")
    df = pd.read_csv(args.metrics_csv, low_memory=False)
    sub = df[df["id"] == args.epitope_id].copy()
    print(f"  {len(sub)} designs for {args.epitope_id}")

    print(f"Loading true complex PDB ...")
    u = mda.Universe(args.true_pdb)

    # --- Get epitope residue indices from dp2 ---
    print(f"Loading dp2 for epitope indices ...")
    dp2 = pd.read_parquet(args.dp2_parquet)
    dp2["assay_scaffolded_epitope_id"] = dp2["assay_scaffolded_epitope_id"].astype(str).str.lower()
    dp2_row = dp2[dp2["id"] == args.epitope_id].iloc[0]
    true_epi_ris = [int(i) for i in dp2_row["epitope_chunk_resindices"]]
    print(f"  Epitope chunk residues: {len(true_epi_ris)}")

    # --- Build epitope coordinate frame ---
    # Use CA atoms of true epitope residues to define the local frame
    epi_ca = u.residues[true_epi_ris].atoms.select_atoms("name CA").positions
    epi_center, epi_axes = pca_frame(epi_ca)
    print(f"  Epitope centroid: {epi_center.round(2)}")

    # All epitope heavy atom positions (for reference surface in plots)
    epi_heavy = u.residues[true_epi_ris].atoms.select_atoms("not name H*").positions
    epi_local = to_local_frame(epi_heavy, epi_center, epi_axes)

    # Antibody atoms (chains B + C)
    ab_atoms = u.select_atoms("(segid B or segid C) and not name H*")
    if len(ab_atoms) == 0:
        ab_atoms = u.select_atoms("(chainid B or chainid C) and not name H*")
    print(f"  Antibody heavy atoms: {len(ab_atoms)}")

    # --- Collect clash atom coordinates ---
    print(f"Collecting clash coordinates ...")
    clash_coords_global = []
    n_with_clash = 0
    n_no_clash   = 0

    for _, row in sub.iterrows():
        ris = parse_resindices(row["af3_clash_resindices"])
        if not ris:
            n_no_clash += 1
            continue
        n_with_clash += 1
        # get heavy atom positions of clashing antibody residues from true PDB
        try:
            clash_atoms = u.residues[ris].atoms.select_atoms("not name H*")
            clash_coords_global.append(clash_atoms.positions)
        except Exception:
            continue

    if not clash_coords_global:
        print("No clash coordinates found — exiting")
        return

    clash_coords_global = np.vstack(clash_coords_global)
    print(f"  Designs with clash: {n_with_clash} / {len(sub)}")
    print(f"  Total clash atom instances: {len(clash_coords_global)}")

    # --- Project into epitope local frame ---
    clash_local = to_local_frame(clash_coords_global, epi_center, epi_axes)

    # Save raw point cloud
    np.save(out_dir / "clash_coords_local.npy", clash_local)
    np.save(out_dir / "epitope_coords_local.npy", epi_local)
    print(f"Saved point clouds to {out_dir}")

    # --- Voxel density grid ---
    print("Building voxel density grid ...")
    vs = args.voxel_size
    pad = 5.0  # Angstrom padding around data

    # define grid bounds from clash + epitope coords combined
    all_local = np.vstack([clash_local, epi_local])
    mins = all_local.min(axis=0) - pad
    maxs = all_local.max(axis=0) + pad

    nx = int(np.ceil((maxs[0] - mins[0]) / vs))
    ny = int(np.ceil((maxs[1] - mins[1]) / vs))
    nz = int(np.ceil((maxs[2] - mins[2]) / vs))
    print(f"  Grid: {nx} x {ny} x {nz} voxels at {vs}Å resolution")

    grid = np.zeros((nx, ny, nz), dtype=float)
    ix = np.floor((clash_local[:, 0] - mins[0]) / vs).astype(int).clip(0, nx - 1)
    iy = np.floor((clash_local[:, 1] - mins[1]) / vs).astype(int).clip(0, ny - 1)
    iz = np.floor((clash_local[:, 2] - mins[2]) / vs).astype(int).clip(0, nz - 1)
    np.add.at(grid, (ix, iy, iz), 1)

    # smooth with gaussian filter
    grid_smooth = gaussian_filter(grid, sigma=args.kde_bw / vs)
    np.save(out_dir / "density_grid.npy", grid_smooth)
    np.save(out_dir / "grid_mins.npy", mins)
    np.save(out_dir / "grid_voxel_size.npy", np.array([vs]))

    # --- 2D slice plots through density ---
    print("Generating 2D slice plots ...")
    fig, axes_2d = plt.subplots(1, 3, figsize=(15, 5))
    fig.suptitle(f"Clash density map — {args.epitope_id}\n"
                 f"({n_with_clash} designs with clashes, {len(clash_coords_global)} atom instances)",
                 fontsize=12)

    axis_labels = ["PC1 (Å)", "PC2 (Å)", "PC3 (Å)"]
    slice_planes = [
        (0, 1, 2, "PC1-PC2 plane (slice at PC3=0)"),
        (0, 2, 1, "PC1-PC3 plane (slice at PC2=0)"),
        (1, 2, 0, "PC2-PC3 plane (slice at PC1=0)"),
    ]

    for ax, (xi, yi, zi, title) in zip(axes_2d, slice_planes):
        # project by summing along the third axis
        density_2d = grid_smooth.sum(axis=zi)
        if zi == 0:
            density_2d = grid_smooth.sum(axis=0)
        elif zi == 1:
            density_2d = grid_smooth.sum(axis=1)
        else:
            density_2d = grid_smooth.sum(axis=2)

        x_edges = mins[xi] + np.arange(density_2d.shape[0] + 1) * vs if zi != 0 else \
                  mins[1] + np.arange(density_2d.shape[0] + 1) * vs
        y_edges = mins[yi] + np.arange(density_2d.shape[1] + 1) * vs if zi != 0 else \
                  mins[2] + np.arange(density_2d.shape[1] + 1) * vs

        im = ax.pcolormesh(
            mins[xi] + np.arange(grid_smooth.shape[xi] + (1 if xi < 2 else 0)) * vs
            if False else np.linspace(mins[xi], maxs[xi], density_2d.shape[0] + 1),
            np.linspace(mins[yi], maxs[yi], density_2d.shape[1] + 1),
            density_2d.T,
            cmap="hot_r", shading="flat"
        )

        # overlay epitope projection
        ax.scatter(epi_local[:, xi], epi_local[:, yi],
                   s=5, c="cyan", alpha=0.4, label="epitope atoms", zorder=3)

        ax.set_xlabel(axis_labels[xi])
        ax.set_ylabel(axis_labels[yi])
        ax.set_title(title, fontsize=9)
        plt.colorbar(im, ax=ax, label="clash density")
        ax.legend(fontsize=7, loc="upper right")

    plt.tight_layout()
    plot_path = out_dir / "clash_density_2d_slices.png"
    plt.savefig(plot_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved 2D slice plot: {plot_path}")

    # --- 3D scatter of top-density clash points ---
    print("Generating 3D scatter plot ...")
    fig3d = plt.figure(figsize=(10, 8))
    ax3d = fig3d.add_subplot(111, projection="3d")

    # subsample clash points for readability
    n_sample = min(5000, len(clash_local))
    idx = np.random.choice(len(clash_local), n_sample, replace=False)
    pts = clash_local[idx]

    sc = ax3d.scatter(pts[:, 0], pts[:, 1], pts[:, 2],
                      c=np.linalg.norm(pts, axis=1),
                      cmap="hot_r", s=2, alpha=0.3, label="clash atoms")
    ax3d.scatter(epi_local[:, 0], epi_local[:, 1], epi_local[:, 2],
                 c="cyan", s=8, alpha=0.6, label="epitope atoms")

    ax3d.set_xlabel("PC1 (Å)")
    ax3d.set_ylabel("PC2 (Å)")
    ax3d.set_zlabel("PC3 (Å)")
    ax3d.set_title(f"Clash atom cloud — {args.epitope_id}\n({n_sample} sampled points)")
    ax3d.legend()
    plt.colorbar(sc, ax=ax3d, label="distance from epitope center (Å)", shrink=0.5)

    plot3d_path = out_dir / "clash_density_3d_scatter.png"
    plt.savefig(plot3d_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved 3D scatter plot: {plot3d_path}")

    # --- Summary stats ---
    print(f"\n=== SUMMARY ===")
    print(f"Epitope: {args.epitope_id}")
    print(f"Total designs:        {len(sub)}")
    print(f"Designs with clash:   {n_with_clash} ({100*n_with_clash/len(sub):.1f}%)")
    print(f"Designs no clash:     {n_no_clash} ({100*n_no_clash/len(sub):.1f}%)")
    print(f"Total clash instances:{len(clash_coords_global)}")
    print(f"Peak density voxel:   {grid_smooth.max():.1f} counts")
    print(f"Output dir: {out_dir}")


if __name__ == "__main__":
    main()
