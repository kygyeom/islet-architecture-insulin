"""Draw one islet structure (PDB file from generate_structures.f90) as a 3D image in
the style of Fig 1(a): alpha cells red, beta cells blue, delta cells green, black
edges, no axes, and a cut-away quarter (cells with x > x_mid and y < y_mid are hidden)
so that the interior is visible.

Clean script version of islets_3D_image.ipynb (kept in this folder).  Same drawing
settings as the notebook: marker size 150, black edges, axis limits x 75-165,
y 55-145, z 55-145 (lattice units), figure 8 x 8 in.  The notebook saved at dpi 1000.

Usage
  python structure/draw_islet_3d.py structure/pdb/core/config_ini_3D4,6.pdb -o results/core_P40_D203.png
  python structure/draw_islet_3d.py FILE.pdb -o out.png --dpi 300 --no-cut
"""
import argparse
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

COLOR = {"C": "red", "S": "blue", "Cd": "green"}      # alpha, beta, delta; other names black


def read_pdb(path):
    """Return atom names and coordinates (fixed PDB columns, as in the notebook)."""
    names, xyz = [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith("ATOM"):
            xyz.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            names.append(line[12:16].strip())
    return np.array(names), np.array(xyz)


def draw(pdb, out, dpi=300, cut=True):
    names, xyz = read_pdb(pdb)
    x_mid, y_mid = xyz[:, 0].mean(), xyz[:, 1].mean()
    keep = ~((xyz[:, 0] > x_mid) & (xyz[:, 1] < y_mid)) if cut else np.ones(len(xyz), bool)
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")
    # one collection for all cells, so matplotlib sorts every marker by depth
    # (the notebook drew one scatter per cell, which gives the same result, slowly)
    colors = [COLOR.get(n, "black") for n in names[keep]]
    ax.scatter(*xyz[keep].T, c=colors, s=150, edgecolors="black", depthshade=False)
    ax.set_axis_off()
    ax.set_xlim(75, 165); ax.set_ylim(55, 145); ax.set_zlim(55, 145)
    fig.savefig(out, dpi=dpi)
    plt.close(fig)
    counts = {k: int((names == k).sum()) for k in np.unique(names)}
    print(f"{pdb}: {len(names)} cells {counts} -> {out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("pdb")
    ap.add_argument("-o", "--out", required=True)
    ap.add_argument("--dpi", type=int, default=300)
    ap.add_argument("--no-cut", action="store_true", help="draw all cells (no cut-away quarter)")
    a = ap.parse_args()
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    draw(a.pdb, a.out, a.dpi, not a.no_cut)


if __name__ == "__main__":
    main()
