"""Figures for README.md.

  results/readme_timeseries.png    glucose and the three hormone outputs of one Core islet
  results/readme_interactions.png  cell-to-cell interaction network of the model

Example run (data/example_trace/): Core islet, P_beta = 60 %, D = 114 um, K = 0.4,
glucose protocol 7.0 -> 15.4 -> 4.2 mM (inputs at 0-60, 60-120 and 140-200 min).
Files are <stem>_islet1_6_2_1_0_ext2g0.dat with columns (time [min], value):
gt = G/G0, ha / hb / hd = mean output per alpha / beta / delta cell.
"""
import pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch

ROOT = pathlib.Path(__file__).resolve().parent.parent
DATA = ROOT / "data" / "example_trace"
OUT = ROOT / "results"
G0 = 7.0                                            # baseline glucose (mM)
C = {"alpha": "#d62728", "beta": "#1f77b4", "delta": "#2ca02c", "glu": "#333333"}
INK, MUTED, GRID = "#222222", "#666666", "#dddddd"


def load(stem):
    t, v = np.loadtxt(DATA / f"{stem}_islet1_6_2_1_0_ext2g0.dat", unpack=True)
    return t, v


def timeseries():
    rows = [("gt", "Glucose (mM)", C["glu"], G0),
            ("ha", "Glucagon\n(α cells)", C["alpha"], 1.0),
            ("hb", "Insulin\n(β cells)", C["beta"], 1.0),
            ("hd", "Somatostatin\n(δ cells)", C["delta"], 1.0)]
    fig, axes = plt.subplots(4, 1, figsize=(7.5, 6.2), sharex=True)
    for ax, (stem, label, col, scale) in zip(axes, rows):
        t, v = load(stem)
        ax.plot(t, v * scale, color=col, lw=1.2)
        ax.set_ylabel(label, fontsize=10, color=INK, rotation=0, ha="right", va="center", labelpad=8)
        for x in (60, 120, 140):
            ax.axvline(x, color=GRID, lw=0.9, zorder=0)
        ax.spines[["top", "right"]].set_visible(False)
        ax.tick_params(labelsize=9, colors=MUTED)
        ax.yaxis.set_major_locator(plt.MaxNLocator(3))
    for x0, x1, txt in [(0, 60, "7.0 mM"), (60, 120, "15.4 mM"), (140, 200, "4.2 mM")]:
        axes[0].text((x0 + x1) / 2, 1.04, txt, transform=axes[0].get_xaxis_transform(),
                     ha="center", va="bottom", fontsize=9, color=MUTED)
    axes[-1].set_xlabel("Time (min)", fontsize=10, color=INK)
    axes[-1].set_xlim(0, 200)
    fig.align_ylabels(axes)
    fig.tight_layout()
    fig.savefig(OUT / "readme_timeseries.png", dpi=200)
    plt.close(fig)


def interactions():
    # Model topology (A_ss' in the Methods): alpha->beta, alpha->delta, beta->delta excitatory (+1),
    # the three reciprocal links inhibitory (-1); no same-type coupling (A_ss = 0).
    pos = {"alpha": (0.0, 0.0), "beta": (2.2, 0.0), "delta": (1.1, 1.75)}
    name = {"alpha": "α", "beta": "β", "delta": "δ"}
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    r = 0.36
    for k, (x, y) in pos.items():
        ax.add_patch(Circle((x, y), r, facecolor=C[k], edgecolor="white", lw=2, zorder=3))
        ax.text(x, y, name[k], ha="center", va="center", fontsize=20, color="white", zorder=4)

    def link(a, b, sign, off):
        (x0, y0), (x1, y1) = pos[a], pos[b]
        d = np.array([x1 - x0, y1 - y0]); L = np.linalg.norm(d); u = d / L; n = np.array([-u[1], u[0]])
        p0 = np.array([x0, y0]) + u * (r + 0.06) + n * off
        p1 = np.array([x1, y1]) - u * (r + 0.06) + n * off
        col = INK if sign > 0 else MUTED
        if sign > 0:                                   # excitatory: arrow head
            ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=16, lw=1.6, color=col, zorder=2))
        else:                                          # inhibitory: dashed line ending in a flat bar
            ax.plot(*zip(p0, p1), color=col, lw=1.6, ls=(0, (4, 2)), zorder=2)
            ax.plot(*zip(p1 - n * 0.11, p1 + n * 0.11), color=col, lw=2.2, zorder=2)
        mid = (p0 + p1) / 2 + n * 0.16 * np.sign(off)
        ax.text(*mid, "+" if sign > 0 else "−", ha="center", va="center", fontsize=13, color=col)

    for a, b in [("alpha", "beta"), ("alpha", "delta"), ("beta", "delta")]:
        link(a, b, +1, 0.12)
        link(b, a, -1, 0.12)
    ax.text(1.1, -0.75, "solid: excitatory (+)    dashed: inhibitory (−)\nno coupling between cells of the same type",
            ha="center", va="top", fontsize=9, color=MUTED)
    ax.set_xlim(-0.6, 2.8); ax.set_ylim(-1.2, 2.3); ax.set_aspect("equal"); ax.axis("off")
    fig.tight_layout()
    fig.savefig(OUT / "readme_interactions.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    timeseries()
    interactions()
    print("saved", OUT / "readme_timeseries.png", OUT / "readme_interactions.png")
