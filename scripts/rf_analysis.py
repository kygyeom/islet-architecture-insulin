"""Step 5 - random-forest analysis behind Fig 3.

Targets
  IPA : one value per (K, S, P_beta, D, G/G0) - 1296 rows.  The normal (1.0) and low
        (0.6) levels exist in every batch and are averaged over the batches.
        Features G, K, Pbeta, D, S  (S = 0 Core, 1 Mix).
  IS  : one value per (K, S, P_beta, D) from data/is_values.csv - 216 rows.
        Features K, Pbeta, D, S.
  The same analyses are repeated within each structure (Fig 3 (c)(d)).

Method (scripts/interaction_importance.py, used without modification)
  * RandomForestRegressor(n_estimators=200, random_state=0), 5-fold GroupKFold with one
    group per row (config_id), held-out R^2.
  * single-feature importance: held-out permutation importance (n_repeats = 30).
  * pair interaction: Friedman-H additive excess of the 2-way partial dependence,
    absolute size in target units (interaction_abs), n_bg = 400 background rows.
  * noise floor: the same pipeline on a target shuffled once across rows (seed 0);
    floor = largest interaction_abs of the shuffled run.  A pair counts as an
    interaction when interaction_abs > 2 x floor.

Output  results/rf/<tag>_single.csv, <tag>_pairs.csv, results/rf/summary.csv
Run time: about 15 min on a laptop (RF uses all cores).
"""
import sys, time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA, RESULTS  # noqa: E402
import interaction_importance as ii  # noqa: E402

KW = dict(n_estimators=200, n_repeats=30, n_bg=400)
SEED = 0
THRESHOLD = 2.0                      # interaction if interaction_abs > 2 x noise floor
OUT = RESULTS / "rf"


def ipa_table(cond):
    cfg = ["K", "struc", "D", "Pbeta"]
    ipa = cond.groupby(cfg + ["G"], as_index=False).IPA.mean()    # average normal/low over batches
    ipa["S"] = (ipa.struc == "mix").astype(int)
    ipa["Pbeta"] = ipa.Pbeta.astype(float); ipa["D"] = ipa.D.astype(float)
    ipa["config_id"] = np.arange(len(ipa))                          # one group per row
    return ipa


def is_table(is_):
    d = is_.copy()
    d["S"] = (d.struc == "mix").astype(int)
    d["Pbeta"] = d.Pbeta.astype(float); d["D"] = d.D.astype(float)
    d["config_id"] = np.arange(len(d))
    return d


def analyse(tag, d, target, feats, summary):
    t0 = time.time()
    single, pairs, r2 = ii.run(d, target, feats, "config_id", seed=SEED, **KW)
    _, null, r2_null = ii.null_floor(d, target, feats, "config_id", seed=SEED, **KW)
    floor = float(null["interaction_abs"].max())
    pairs["above_null"] = pairs["interaction_abs"] / floor
    pairs["interaction"] = pairs["above_null"] > THRESHOLD
    single.to_csv(OUT / f"{tag}_single.csv", index=False)
    pairs.to_csv(OUT / f"{tag}_pairs.csv", index=False)
    summary.append(dict(tag=tag, rows=len(d), heldout_r2=r2, null_r2=r2_null, noise_floor=floor,
                        interactions=";".join(f"{a}x{b}" for a, b in pairs.loc[pairs.interaction, ["A", "B"]].values)))
    print(f"=== {tag}: rows={len(d)}  held-out R2={r2:.3f}  null R2={r2_null:.3f}  "
          f"floor={floor:.4f}  ({time.time() - t0:.0f} s)")
    print(single.round(4).to_string(index=False))
    print(pairs[["A", "B", "interaction_abs", "above_null", "interaction"]].round(3).to_string(index=False), flush=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    cond = pd.read_csv(DATA / "conditions.csv", dtype={"batch": str})
    ipa = ipa_table(cond)
    is_ = is_table(pd.read_csv(DATA / "is_values.csv"))
    summary = []
    analyse("IPA", ipa, "IPA", ["G", "K", "Pbeta", "D", "S"], summary)
    analyse("IS", is_, "IS", ["K", "Pbeta", "D", "S"], summary)
    for st in ("core", "mix"):
        analyse(f"IPA_{st}", ipa[ipa.struc == st], "IPA", ["G", "K", "Pbeta", "D"], summary)
        analyse(f"IS_{st}", is_[is_.struc == st], "IS", ["K", "Pbeta", "D"], summary)
    pd.DataFrame(summary).to_csv(OUT / "summary.csv", index=False)
    print("saved", OUT / "summary.csv")


if __name__ == "__main__":
    main()
