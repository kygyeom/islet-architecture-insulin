"""Step 2 - extract the IPA (insulin-production activity) of every simulation run.

IPA of one run in one glucose phase = time average of hb(t) over that phase, where
hb(t) is the mean output of the beta cells, (1/N_beta) * sum_j A_j (1 + cos theta_j),
written by the Fortran code to  hb_islet<nt>_<pb>_<sz>_<ep>_<kd>_ext2g0.dat.

    normal  rows [0, 600)     t =   0-60  min   G/G0 label 1.0
    high    rows [600, 1200)  t =  60-120 min   G/G0 label 1 + 2X
    low     rows [1400, 2000) t = 140-200 min   G/G0 label 0.6
    (rows 1200-1399, t = 120-140 min, are the recovery interval and are not used)

This is the same window average as the original analysis notebook (function
`counsumption`, beta / 'consum').

Usage
    python scripts/extract_ipa.py RUN_DIR --batch 1204 --x 0.6 -o results/ipa_1204.csv

RUN_DIR is searched recursively; files under a folder named DATA_c are Core runs,
files under DATA are Mix runs (the Fortran programs write there).
"""
import argparse, re, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import CELLS, D_UM, PHASES, g_label, k_from_ep  # noqa: E402

NAME = re.compile(r"hb_islet(\d+)_(\d+)_(\d+)_(\d+)_(\d+)_ext2g0\.dat$")
NROWS = 2000                      # itmax / 10 = 20000 / 10
COLUMNS = ["batch", "X", "struc", "K", "Pbeta", "size_idx", "ncell", "D",
           "nt", "phase", "G", "IPA", "source"]


def ipa_of_file(path):
    """Return {phase: window mean of column 2} for one hb file."""
    hb = np.loadtxt(path, usecols=1)
    if hb.shape[0] != NROWS:
        raise ValueError(f"{path}: {hb.shape[0]} rows, expected {NROWS} (run incomplete?)")
    return {ph: float(hb[a:b].mean()) for ph, (a, b, _) in PHASES.items()}


def extract(run_dir, batch, x, source="this_run"):
    rows = []
    for f in sorted(Path(run_dir).rglob("hb_islet*_ext2g0.dat")):
        struc = {"DATA_c": "core", "DATA": "mix"}.get(f.parent.name)
        m = NAME.search(f.name)
        if struc is None or m is None:
            continue
        nt, pb, sz, ep, kd = map(int, m.groups())
        if kd != 0:
            continue                                  # self coupling is 0 in the paper
        for phase, ipa in ipa_of_file(f).items():
            rows.append(dict(batch=batch, X=x, struc=struc, K=k_from_ep(ep), Pbeta=10 * pb,
                             size_idx=sz, ncell=CELLS[sz], D=D_UM[sz], nt=nt, phase=phase,
                             G=g_label(phase, x), IPA=ipa, source=source))
    if not rows:
        raise SystemExit(f"no hb_islet*_ext2g0.dat files found under {run_dir}")
    df = pd.DataFrame(rows, columns=COLUMNS)
    dup = df.duplicated(["struc", "K", "Pbeta", "D", "nt", "phase"])
    if dup.any():
        raise SystemExit(f"{dup.sum()} duplicated runs (same condition and nt) under {run_dir}")
    return df.sort_values(["struc", "K", "D", "Pbeta", "G", "nt"]).reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--batch", required=True, help="batch id, e.g. 1204")
    ap.add_argument("--x", type=float, required=True, help="high-phase coefficient X of this batch")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    df = extract(a.run_dir, a.batch, a.x)
    Path(a.out).parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(a.out, index=False)
    n_runs = len(df) // len(PHASES)
    print(f"{n_runs} runs -> {len(df)} rows ({df.struc.value_counts().to_dict()} rows by structure) -> {a.out}")


if __name__ == "__main__":
    main()
