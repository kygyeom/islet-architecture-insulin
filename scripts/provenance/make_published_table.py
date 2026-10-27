"""Provenance only - how data/ipa_realizations.csv was assembled.

You do NOT need to run this script to reproduce the figures; it documents how the
published table was built from the authors' intermediate files, which are not part
of this repository.

Sources
  * Batches 1009 (X=0.8), 1011 (X=0.4), 1204 (X=0.6): legacy analysis tables
    result_<batch>_core_data.csv (Core) and result_<batch>_data.csv (Mix), written by
    the original analysis notebook from the hb_* time series.
      - Only rows type == 'beta' and re_type == 'consum' are the IPA.
      - Their G column uses 1.0 = normal, 2 = high, 0.5 = low phase (phase codes,
        not G/G0); they are relabelled 1.0, 1 + 2X and 0.6.
      - IMPORTANT: the row with nt = k holds the CUMULATIVE MEAN of realizations
        1..k (the notebook appended each realization to a list and stored its mean).
        Individual realizations are recovered as
            r1 = v1,   r2 = 2*v2 - v1,   r3 = 3*v3 - 2*v2.
      - Rows are filtered by the 'struc' column (the legacy Mix file of batch 0123
        contained Core rows; batch 0123 is not taken from legacy tables at all).
  * Batch 0123 (X=0.2): re-run on 2026-10-27 with the paper protocol
    (low phase -0.2, 3 realizations as separate jobs with distinct seeds) and
    extracted from the hb_* files with scripts/extract_ipa.py.

Usage
  python scripts/provenance/make_published_table.py LEGACY_DIR IPA_0123_CSV -o data/ipa_realizations.csv
"""
import argparse, sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from common import BATCHES, CELLS, D_UM, g_label  # noqa: E402
from extract_ipa import COLUMNS  # noqa: E402

LEGACY_D = {92: 1, 114: 2, 137: 3, 159: 4, 181: 5, 202: 6}     # legacy diameter label -> sml_size
PHASE_CODE = {1.0: "normal", 2.0: "high", 0.5: "low"}


def legacy_batch(legacy_dir, batch):
    x = BATCHES[batch]
    out = []
    for struc, fn in (("core", f"result_{batch}_core_data.csv"), ("mix", f"result_{batch}_data.csv")):
        d = pd.read_csv(Path(legacy_dir) / fn)
        d = d[(d.struc == struc) & (d.type == "beta") & (d.re_type == "consum") & (d.self == 0)].copy()
        d["phase"] = d.G.map(PHASE_CODE)
        assert d.phase.notna().all()
        key = ["eps", "diameter", "per", "phase"]
        v = d.pivot_table(index=key, columns="nt", values="result1", aggfunc="first")
        assert v.notna().all().all() and list(v.columns) == [1, 2, 3], "need nt = 1, 2, 3 for every condition"
        r = pd.DataFrame({1: v[1], 2: 2 * v[2] - v[1], 3: 3 * v[3] - 2 * v[2]})   # undo cumulative mean
        r = r.stack().rename("IPA").reset_index().rename(columns={"level_4": "nt"})
        r["nt"] = r["nt"].astype(int)
        sz = r.diameter.map(LEGACY_D)
        out.append(pd.DataFrame(dict(
            batch=batch, X=x, struc=struc, K=r.eps.round(1), Pbeta=(10 * r.per).astype(int),
            size_idx=sz, ncell=sz.map(CELLS), D=sz.map(D_UM), nt=r.nt, phase=r.phase,
            G=[g_label(p, x) for p in r.phase], IPA=r.IPA, source="legacy_table_cumulative_mean_undone")))
    return pd.concat(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("legacy_dir")
    ap.add_argument("ipa_0123_csv", help="output of extract_ipa.py for the 2026-10-27 re-run of batch 0123")
    ap.add_argument("-o", "--out", required=True)
    a = ap.parse_args()
    parts = [legacy_batch(a.legacy_dir, b) for b in ("1009", "1011", "1204")]
    b0 = pd.read_csv(a.ipa_0123_csv, dtype={"batch": str}, float_precision="round_trip")
    b0["batch"] = "0123"; b0["source"] = "rerun_2026-10-27"
    parts.append(b0[COLUMNS])
    df = pd.concat(parts, ignore_index=True)[COLUMNS]
    df["G"] = df.G.round(1); df["K"] = df.K.round(1)
    df = df.sort_values(["batch", "struc", "K", "D", "Pbeta", "G", "nt"]).reset_index(drop=True)
    # 4 batches x 2 structures x 3 K x 6 D x 6 Pbeta x 3 phases x 3 realizations
    assert len(df) == 4 * 2 * 3 * 6 * 6 * 3 * 3, len(df)
    assert not df.duplicated(["batch", "struc", "K", "D", "Pbeta", "phase", "nt"]).any()
    assert np.isfinite(df.IPA).all()
    df.to_csv(a.out, index=False)
    print(f"{len(df)} rows -> {a.out}")
    print(df.groupby(["batch", "source"]).size().to_string())


if __name__ == "__main__":
    main()
