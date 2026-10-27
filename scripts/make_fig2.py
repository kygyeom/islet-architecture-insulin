"""Step 4 - Fig 2: IPA and IS of Core and Mix islets.

  (a) IPA vs G/G0 at K = 0.4
  (b) IPA vs coupling K (all six glucose levels, equal weight)
  (c) IS vs P_beta (K = 0.4, P_beta 20-60 %)
  (d) IS vs D      (K = 0.4, D 114-203 um)

Bars / lines: mean over the 36 islet configurations (6 D x 6 P_beta) of one structure;
error bars / bands: standard deviation over those configurations (ddof = 1).
Each configuration value is the realization mean (data/conditions.csv).  The normal
(1.0) and low (0.6) levels are simulated in all four batches; they are averaged over
the batches per configuration, so every glucose level contributes one value per
configuration.  IS is taken from data/is_values.csv (one batch, X = 0.6).

Output  results/fig2.png, results/fig2_values.csv (the plotted means and SDs)
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, RESULTS  # noqa: E402

G_LEVELS = [0.6, 1.0, 1.4, 1.8, 2.2, 2.6]
KS = [0.0, 0.4, 0.8]
STRUCS = ["core", "mix"]
COLOR = {"core": "#1f77b4", "mix": "#ff7f0e"}
LABEL = {"core": "Core", "mix": "Mix"}
FS = 20


def level_vals(cond, struc, g):
    """One IPA per configuration (D, P_beta) for one structure and glucose level."""
    d = cond[(cond.struc == struc) & np.isclose(cond.G, g)]
    return d.groupby(["D", "Pbeta"]).IPA.mean().values          # 36 values


def main():
    cond = pd.read_csv(DATA / "conditions.csv", dtype={"batch": str})
    is_ = pd.read_csv(DATA / "is_values.csv")
    is_ = is_[np.isclose(is_.K, 0.4)]
    RESULTS.mkdir(exist_ok=True)
    rows = []

    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    plt.subplots_adjust(hspace=0.45, wspace=0.42)
    w = 0.35

    # (a) K = 0.4, IPA vs G/G0
    ax = axes[0, 0]; k4 = cond[np.isclose(cond.K, 0.4)]; x = np.arange(len(G_LEVELS))
    for i, s in enumerate(STRUCS):
        v = [level_vals(k4, s, g) for g in G_LEVELS]
        m = [a.mean() for a in v]; sd = [a.std(ddof=1) for a in v]
        ax.bar(x + (i - 0.5) * w, m, w, yerr=sd, capsize=4, color=COLOR[s], label=LABEL[s])
        rows += [dict(panel="a", struc=s, x=g, mean=mm, sd=ss, n=len(a)) for g, mm, ss, a in zip(G_LEVELS, m, sd, v)]
    ax.set_xticks(x); ax.set_xticklabels([f"{g}" for g in G_LEVELS])
    ax.set_xlabel(r"$G/G_0$", fontsize=FS); ax.set_ylabel("IPA", fontsize=FS)
    ax.set_ylim(0, 1.35); ax.tick_params(labelsize=FS)
    ax.legend(loc="upper right", ncol=2, fontsize=15, frameon=False, handlelength=1.2, columnspacing=0.8)

    # (b) IPA vs K, six glucose levels pooled with equal weight (216 configurations x levels)
    ax = axes[0, 1]; x = np.arange(len(KS))
    for i, s in enumerate(STRUCS):
        m, sd = [], []
        for K in KS:
            v = np.concatenate([level_vals(cond[np.isclose(cond.K, K)], s, g) for g in G_LEVELS])
            m.append(v.mean()); sd.append(v.std(ddof=1))
            rows.append(dict(panel="b", struc=s, x=K, mean=m[-1], sd=sd[-1], n=len(v)))
        ax.bar(x + (i - 0.5) * w, m, w, yerr=sd, capsize=4, color=COLOR[s])
    ax.set_xticks(x); ax.set_xticklabels([f"{k}" for k in KS])
    ax.set_xlabel(r"$K$", fontsize=FS); ax.set_ylabel("IPA", fontsize=FS)
    ax.set_ylim(0, 1.35); ax.tick_params(labelsize=FS)

    # (c)(d) IS vs P_beta and vs D, mean +- SD over the other parameter
    def line(ax, key, lo, hi, xlabel, panel):
        for s in STRUCS:
            d = (is_[(is_.struc == s) & (is_[key] >= lo) & (is_[key] <= hi)]
                 .groupby(key).IS.agg(["mean", "std", "size"]).reset_index())
            ax.plot(d[key], d["mean"], marker="o", color=COLOR[s], label=LABEL[s])
            ax.fill_between(d[key], d["mean"] - d["std"], d["mean"] + d["std"], alpha=0.2, color=COLOR[s], linewidth=0)
            rows.extend(dict(panel=panel, struc=s, x=r[key], mean=r["mean"], sd=r["std"], n=r["size"]) for _, r in d.iterrows())
        ax.set_xlabel(xlabel, fontsize=FS); ax.set_ylabel("IS", fontsize=FS)
        ax.set_ylim(-0.03, 0.43); ax.tick_params(labelsize=FS)

    line(axes[1, 0], "Pbeta", 20, 60, r"$P_{\beta}$ (%)", "c")
    axes[1, 0].legend(loc="upper right", ncol=2, fontsize=15, frameon=False, columnspacing=0.8,
                      handlelength=1.5, borderpad=0.2, labelspacing=0.3)
    line(axes[1, 1], "D", 114, 203, "D (μm)", "d")

    for ax, lab in zip(axes.flatten(), ["(a)", "(b)", "(c)", "(d)"]):
        ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        ax.text(0.02, 0.97, lab, transform=ax.transAxes, fontsize=FS, va="top")
    fig.savefig(RESULTS / "fig2.png", dpi=300, bbox_inches="tight")

    vals = pd.DataFrame(rows)
    vals.to_csv(RESULTS / "fig2_values.csv", index=False)
    a = vals[vals.panel == "a"].pivot(index="x", columns="struc", values="mean").round(2)
    print("Fig 2(a) IPA at K = 0.4 (mean over 36 configurations):"); print(a.to_string())
    print("saved", RESULTS / "fig2.png")


if __name__ == "__main__":
    main()
