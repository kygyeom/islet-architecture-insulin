
"""interaction_importance.py
Held-out permutation importance (replaces MDI) + Friedman H interaction strength
(additive excess of the 2-way partial dependence over the two 1-way ones).
Friedman & Popescu (2008), Ann. Appl. Stat. 2, 916.
"""
import numpy as np, pandas as pd, itertools
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold
from sklearn.metrics import r2_score


# ---------- 1. single-feature permutation importance (held-out) ----------
def _perm_imp(model, X, y, col, n_repeats, rng):
    base = r2_score(y, model.predict(X))
    out = []
    for _ in range(n_repeats):
        Xp = X.copy()
        Xp[col] = X[col].to_numpy()[rng.permutation(len(X))]
        out.append(base - r2_score(y, model.predict(Xp)))
    return float(np.mean(out))


# ---------- 2. Friedman H: additive excess of 2-way partial dependence ----------
def _pd(model, Xbg, cols, grids):
    """Brute-force partial dependence of `model` on `cols` over the cartesian `grids`."""
    pts = list(itertools.product(*grids))
    out = np.empty(len(pts))
    for i, pt in enumerate(pts):
        Xt = Xbg.copy()
        for c, v in zip(cols, pt):
            Xt[c] = v
        out[i] = model.predict(Xt).mean()
    return out.reshape([len(g) for g in grids])


def _h_pair(model, Xbg, a, b):
    ga, gb = np.sort(Xbg[a].unique()), np.sort(Xbg[b].unique())
    Fab = _pd(model, Xbg, [a, b], [ga, gb])
    Fa  = _pd(model, Xbg, [a], [ga])
    Fb  = _pd(model, Xbg, [b], [gb])
    Fab = Fab - Fab.mean(); Fa = Fa - Fa.mean(); Fb = Fb - Fb.mean()
    excess = Fab - Fa[:, None] - Fb[None, :]          # <-- the additive excess
    denom = np.sum(Fab ** 2)
    return (float(np.sum(excess ** 2) / denom) if denom > 0 else np.nan,   # H^2 (0-1)
            float(np.sqrt(np.mean(excess ** 2))))                          # absolute, target units


# ---------- 3. driver ----------
def run(df, target, features, group_col,
        n_estimators=200, n_repeats=30, n_bg=400, seed=0, n_splits=5):
    rng = np.random.default_rng(seed)
    X, y, g = df[features], df[target].to_numpy(), df[group_col].to_numpy()
    singles = {f: [] for f in features}
    H2  = {p: [] for p in itertools.combinations(features, 2)}
    Hab = {p: [] for p in itertools.combinations(features, 2)}
    r2 = []
    for tr, te in GroupKFold(n_splits=n_splits).split(X, y, groups=g):
        m = RandomForestRegressor(n_estimators=n_estimators, random_state=seed, n_jobs=-1)
        m.fit(X.iloc[tr], y[tr])
        Xte, yte = X.iloc[te], y[te]
        r2.append(r2_score(yte, m.predict(Xte)))
        for f in features:
            singles[f].append(_perm_imp(m, Xte, yte, f, n_repeats, rng))
        bg = Xte.sample(min(n_bg, len(Xte)), random_state=seed)
        for a, b in H2:
            h2, hab = _h_pair(m, bg, a, b)
            H2[(a, b)].append(h2); Hab[(a, b)].append(hab)
    s = (pd.DataFrame({"feature": list(singles),
                       "perm_importance": [np.mean(v) for v in singles.values()],
                       "sd": [np.std(v, ddof=1) for v in singles.values()]})
         .sort_values("perm_importance", ascending=False))
    rows = []
    for k in H2:
        h2, hab = np.array(H2[k]), np.array(Hab[k])
        rows.append({"A": k[0], "B": k[1],
                     "H2": h2.mean(), "H2_sd": h2.std(ddof=1),
                     "interaction_abs": hab.mean(), "interaction_abs_sd": hab.std(ddof=1),
                     "significant": hab.mean() - 2 * hab.std(ddof=1) / np.sqrt(len(hab)) > 0})
    p = pd.DataFrame(rows).sort_values("H2", ascending=False)
    return s, p, float(np.mean(r2))


# ---------- 4. null floor: same pipeline on a shuffled target ----------
def null_floor(df, target, features, group_col, seed=0, **kw):
    """Re-run `run()` with the target shuffled across configurations.
    The resulting `interaction_abs` column is the noise floor an interaction
    must clear to be reported."""
    d = df.copy()
    d[target] = d[target].sample(frac=1.0, random_state=seed).to_numpy()
    return run(d, target, features, group_col, seed=seed, **kw)
