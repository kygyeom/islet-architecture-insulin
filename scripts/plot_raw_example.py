"""Plot the raw simulation output in data/raw_example/ (one Core and one Mix islet).

    python scripts/plot_raw_example.py            -> results/raw_example.png

Each folder holds the four raw files of one 200-min run:
    gt_*.dat  glucose G/G0
    ha_*.dat  glucagon  (mean output per alpha cell)
    hb_*.dat  insulin   (mean output per beta cell)
    hd_*.dat  somatostatin (mean output per delta cell)
Two columns: time (min), value.  One line every 0.1 min (2000 lines).
"""
import pathlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw_example"
G0 = 7.0
RUNS = [("core", "Core islet  ($P_\\beta$ 60 %, D 114 µm)"), ("mix", "Mix islet  ($P_\\beta$ 40 %, D 203 µm)")]
ROWS = [("gt", "Glucose (mM)", "#333333", G0), ("ha", "Glucagon", "#d62728", 1.0),
        ("hb", "Insulin", "#1f77b4", 1.0), ("hd", "Somatostatin", "#2ca02c", 1.0)]


def load(folder, stem):
    f = next((RAW / folder).glob(f"{stem}_*.dat"))
    return np.loadtxt(f, unpack=True)


fig, axes = plt.subplots(len(ROWS), 2, figsize=(10, 6.5), sharex=True, sharey="row")
for j, (folder, title) in enumerate(RUNS):
    axes[0, j].set_title(title, fontsize=10)
    for i, (stem, label, col, scale) in enumerate(ROWS):
        t, v = load(folder, stem)
        ax = axes[i, j]
        ax.plot(t, v * scale, color=col, lw=1.0)
        for x in (60, 120, 140):
            ax.axvline(x, color="#dddddd", lw=0.8, zorder=0)
        ax.spines[["top", "right"]].set_visible(False)
        if j == 0:
            ax.set_ylabel(label, fontsize=9)
for ax in axes[-1]:
    ax.set_xlabel("Time (min)")
    ax.set_xlim(0, 200)
fig.tight_layout()
out = ROOT / "results" / "raw_example.png"
fig.savefig(out, dpi=150)
print("saved", out)
