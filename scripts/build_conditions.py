"""Step 3 - average the realizations into one IPA per condition and compute the IS.

Input   data/ipa_realizations.csv  (one row per run x glucose phase; columns described
        in README.md).  Further extract_ipa.py outputs can be appended with --extra.
Output  data/conditions.csv   batch, X, struc, K, Pbeta, D, ncell, phase, G,
                              IPA (mean over realizations), IPA_sd, n_real
        data/is_values.csv    struc, K, Pbeta, D, IS
                              IS = [IPA(G/G0 = 2.2) - IPA(G/G0 = 0.6)] / (2.2 - 0.6),
                              both taken from the SAME batch (default 1204, X = 0.6),
                              so that the high and low phases come from one simulation.

A "condition" is (batch, structure, K, P_beta, D, glucose phase).  The normal (1.0) and
low (0.6) phases are simulated in every batch; downstream scripts average them over the
batches (see rf_analysis.py / make_fig2.py).
"""
import argparse, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import DATA  # noqa: E402

KEY = ["batch", "X", "struc", "K", "Pbeta", "D", "ncell", "phase", "G"]


def read_realizations(paths):
    df = pd.concat([pd.read_csv(p, dtype={"batch": str}, float_precision="round_trip") for p in paths],
                   ignore_index=True)
    df["batch"] = df.batch.str.zfill(4)
    df["G"] = df.G.round(1); df["K"] = df.K.round(1)
    return df


def conditions(real):
    c = (real.groupby(KEY, as_index=False)
             .agg(IPA=("IPA", "mean"), IPA_sd=("IPA", "std"), n_real=("IPA", "size")))
    return c.sort_values(["batch", "struc", "K", "D", "Pbeta", "G"]).reset_index(drop=True)


def insulin_sensitivity(cond, batch="1204", g_hi=2.2, g_lo=0.6):
    b = cond[cond.batch == batch]
    cfg = ["struc", "K", "D", "Pbeta"]
    hi = b[np.isclose(b.G, g_hi)].set_index(cfg).IPA
    lo = b[np.isclose(b.G, g_lo)].set_index(cfg).IPA
    if hi.empty or lo.empty:
        raise SystemExit(f"batch {batch} has no G/G0 = {g_hi} and/or {g_lo} rows")
    is_ = ((hi - lo) / (g_hi - g_lo)).rename("IS").reset_index()
    return is_.sort_values(cfg).reset_index(drop=True)[["struc", "K", "Pbeta", "D", "IS"]]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--realizations", default=str(DATA / "ipa_realizations.csv"))
    ap.add_argument("--extra", nargs="*", default=[], help="additional extract_ipa.py outputs")
    ap.add_argument("--is-batch", default="1204")
    ap.add_argument("--outdir", default=str(DATA))
    a = ap.parse_args()
    real = read_realizations([a.realizations] + a.extra)
    dup = real.duplicated(["batch", "struc", "K", "D", "Pbeta", "phase", "nt"])
    if dup.any():
        raise SystemExit(f"{dup.sum()} duplicated (batch, condition, nt) rows")
    cond = conditions(real)
    is_ = insulin_sensitivity(cond, a.is_batch)
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    cond.to_csv(out / "conditions.csv", index=False)
    is_.to_csv(out / "is_values.csv", index=False)
    print(f"{len(real)} realization rows -> {len(cond)} conditions "
          f"(n_real: {cond.n_real.value_counts().to_dict()}) -> {out / 'conditions.csv'}")
    print(f"IS from batch {a.is_batch}: {len(is_)} configurations -> {out / 'is_values.csv'}")


if __name__ == "__main__":
    main()
