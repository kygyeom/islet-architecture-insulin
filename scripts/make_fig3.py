"""Step 6 - Fig 3: random-forest feature importance and pair interactions.

  (a) single-feature permutation importance for IPA (all data)
  (b) the same for IS
  (c) IPA importance within Core and within Mix
  (d) IS importance within Core and within Mix
  (e) pair interaction (Friedman-H additive excess, target units) for IPA
  (f) the same for IS
  Error bars: SD over the 5 cross-validation folds.  Dashed line in (e)(f): 2 x noise
  floor (interaction criterion; floor from results/rf/summary.csv).

Input   results/rf/*.csv  (run rf_analysis.py first)
Output  results/fig3.png
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import RESULTS  # noqa: E402

R = RESULTS / "rf"
FS = 20
LAB = {"G": "G", "K": "K", "Pbeta": r"$P_{\beta}$", "D": "D", "S": "S"}
EK = dict(lw=1.2, capsize=3, ecolor="black")


def style(ax, xlabel, ylabel):
    ax.set_xlabel(xlabel, fontsize=FS, labelpad=6); ax.set_ylabel(ylabel, fontsize=FS, labelpad=6)
    ax.xaxis.set_major_locator(MaxNLocator(nbins=5))
    ax.tick_params(axis="x", labelsize=FS - 2, pad=2); ax.tick_params(axis="y", labelsize=FS, pad=2)
    ax.set_xlim(left=0)


def xerr0(v, e):
    """Clip the lower error bar at zero so the axis does not extend to negative values."""
    v = np.clip(np.asarray(v, float), 0, None); e = np.asarray(e, float)
    return [np.minimum(e, v), e]


def series(ax, df, ylabel):
    df = df.sort_values("perm_importance")                      # largest on top
    y = np.arange(len(df))
    ax.barh(y, df["perm_importance"].clip(lower=0), xerr=xerr0(df["perm_importance"], df["sd"]), error_kw=EK)
    ax.set_yticks(y); ax.set_yticklabels([LAB[f] for f in df["feature"]], fontsize=FS)
    ax.set_ylim(-0.55, len(df) + 0.5); style(ax, "Importance", ylabel)


def grouped(ax, d, order, ylabel, legend=False):
    y = np.arange(len(order)); h = 0.32
    for lab, off in (("Core", -h / 2), ("Mix", h / 2)):
        s = d[d["structure"] == lab].set_index("feature").reindex(order)
        ax.barh(y + off, s["perm_importance"].clip(lower=0), height=h,
                xerr=xerr0(s["perm_importance"], s["sd"]), error_kw=EK, label=lab)
    ax.set_yticks(y); ax.set_yticklabels([LAB[f] for f in order], fontsize=FS)
    if legend:
        ax.legend(fontsize=FS - 2, frameon=False, loc="lower right")
    ax.set_ylim(-0.85, len(order) - 0.45); ax.invert_yaxis(); style(ax, "Importance", ylabel)


def pairs(ax, tag, floor, ylabel):
    p = pd.read_csv(R / f"{tag}_pairs.csv").sort_values("interaction_abs")
    y = np.arange(len(p))
    ax.barh(y, p["interaction_abs"], height=0.6,
            xerr=xerr0(p["interaction_abs"], p["interaction_abs_sd"]), error_kw=EK)
    ax.axvline(2 * floor, ls="--", color="gray", lw=1.5)          # interaction criterion
    ax.set_yticks(y); ax.set_yticklabels([f"{{{LAB[a]}, {LAB[b]}}}" for a, b in zip(p.A, p.B)], fontsize=FS - 2)
    ax.set_ylim(-0.55, len(p) + 0.6); style(ax, "Importance", ylabel)
    ax.tick_params(axis="y", labelsize=FS - 2)


def main():
    floor = pd.read_csv(R / "summary.csv").set_index("tag").noise_floor
    fig = plt.figure(figsize=(3532 / 300, 3726 / 300)); gs = gridspec.GridSpec(3, 2, figure=fig)
    bys = pd.concat([pd.read_csv(R / f"{tg}_{st}_single.csv").assign(structure=st.capitalize(), target=tg)
                     for tg in ("IPA", "IS") for st in ("core", "mix")])
    # rows: all data / by structure / pairs; columns: IPA (left), IS (right)
    ax_a = fig.add_subplot(gs[0, 0]); series(ax_a, pd.read_csv(R / "IPA_single.csv"), "Feature for IPA")
    ax_b = fig.add_subplot(gs[0, 1]); series(ax_b, pd.read_csv(R / "IS_single.csv"), "Feature for IS")
    ax_c = fig.add_subplot(gs[1, 0]); grouped(ax_c, bys[bys.target == "IPA"], ["G", "K", "Pbeta", "D"], "Feature for IPA", legend=True)
    ax_d = fig.add_subplot(gs[1, 1]); grouped(ax_d, bys[bys.target == "IS"], ["K", "Pbeta", "D"], "Feature for IS")
    ax_e = fig.add_subplot(gs[2, 0]); pairs(ax_e, "IPA", floor["IPA"], "Feature Pair for IPA")
    ax_f = fig.add_subplot(gs[2, 1]); pairs(ax_f, "IS", floor["IS"], "Feature Pair for IS")
    for ax, lab in [(ax_a, "(a)"), (ax_b, "(b)"), (ax_c, "(c)"), (ax_d, "(d)"), (ax_e, "(e)"), (ax_f, "(f)")]:
        ax.text(0.02, 0.98, lab, transform=ax.transAxes, fontsize=FS + 2, va="top", ha="left")
    gs.update(left=0.155, right=0.985, bottom=0.065, top=0.995, wspace=0.55, hspace=0.36)
    # left-align the pair labels of (e)(f): pad = width of the longest label
    fig.canvas.draw(); rend = fig.canvas.get_renderer()
    for ax in (ax_e, ax_f):
        labs = ax.get_yticklabels()
        wmax = max(l.get_window_extent(rend).width for l in labs) * 72 / fig.dpi
        ax.tick_params(axis="y", pad=wmax + 6)
        for l in labs:
            l.set_ha("left")
    fig.savefig(RESULTS / "fig3.png", dpi=300)
    print("saved", RESULTS / "fig3.png")


if __name__ == "__main__":
    main()
