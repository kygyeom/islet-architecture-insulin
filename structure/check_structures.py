"""Checks on the islet structure files (PDB) written by generate_structures.f90.

For every file  <root>/{core,mix}/config_ini_3D<b_ratio>,<sml_size>.pdb  it reports:
  1. cell number n, against an independent re-computation of the HCP lattice
     (same constants as the Fortran code, single precision) and the expected
     743, 1405, 2403, 3791, 5563, 7881;
  2. duplicate coordinates (expected 0);
  3. largest distance from the islet centre in cell diameters, against the radius
     is_size = 1.2*sml_size + 3.8 (all cells must lie inside).  Sites that lie exactly
     on the sphere (sizes 1 and 6: r = 5 and 11) are inside in the Fortran code because
     of single-precision rounding; with the 3-decimal PDB coordinates they read as
     5.0001 / 11.0001, so a tolerance of 1e-3 is used;
  4. sphere fill: n / [sqrt(2) * (4 pi / 3) * is_size^3]  (HCP packing has sqrt(2)
     cells per cubic cell diameter);
  5. alpha / beta / delta counts against the prescribed fractions
     (n_alpha = int(f_alpha*n), n_delta = int(f_delta*n), n_beta = the rest);
  6. Core only: every beta cell is at least as close to the centre as every non-beta
     cell (beta cells fill the innermost sites);
  7. neighbour composition: fraction of beta cells that have no non-beta neighbour
     (they are not coupled to any alpha or delta cell, since the model has only
     heterotypic coupling), and the mean number of non-beta neighbours per beta cell.
     Neighbours: the 12 HCP nearest neighbours (distance 4 lattice units = one cell
     diameter), the same pairs as neighbor_find in the Fortran code.

Usage
  python structure/check_structures.py structure/pdb [--delta 0.2]
--delta: the delta fraction the files were generated with (0.2 for the archived files,
0.1 for files made with DELTA = 0.1).
Atom names: C = alpha, S = beta, Cd = delta, H = cell number 2000 written with a marker
name instead of its type (legacy behaviour of the generator; counted separately).
"""
import argparse, re, sys
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree      # installed with scikit-learn

f32 = np.float32
A = 4.0                                 # lattice constant = one cell diameter
NATX, NATY, NATZ = 60, 60, 30
S3 = float(f32(np.sqrt(f32(3.0))))
S83 = float(f32(np.sqrt(f32(8.0) / f32(3.0))))
S23 = float(f32(np.sqrt(f32(2.0) / f32(3.0))))
CENTRE = np.array([A * NATX / 2, A * S3 / 2 * NATY / 2, A * S83 * NATZ / 2])
EXPECTED_N = {1: 743, 2: 1405, 3: 2403, 4: 3791, 5: 5563, 6: 7881}
D_UM = {s: int(round(2 * 5.1 * n ** (1 / 3))) for s, n in EXPECTED_N.items()}
NAME = re.compile(r"config_ini_3D(\d+),(\d+)\.pdb$")


def is_size(s):
    """is_size = sml_size*1.2 + 3.8, evaluated in single precision like the Fortran code."""
    return float(f32(f32(s) * f32(1.2)) + f32(3.8))


def lattice_count(s):
    """Independent count of the HCP sites inside radius is_size (same loops as init_structure)."""
    i, j, k = np.meshgrid(np.arange(1, NATX + 1), np.arange(1, NATY + 1), np.arange(1, NATZ + 1), indexing="ij")
    x0, y0, z0 = (i - 1) * A, (j - 1) * A * S3 / 2, (k - 1) * A * S83
    odd = (j % 2) == 1
    p1 = np.stack([np.where(odd, x0, x0 + A / 2), y0, z0], -1).reshape(-1, 3)
    p2 = np.stack([np.where(odd, x0 + A / 2, x0 + A), y0 + A * S3 / 6, z0 + A * S23], -1).reshape(-1, 3)
    r = np.linalg.norm(np.vstack([p1, p2]) - CENTRE, axis=1) / A
    return int((r < is_size(s)).sum())


def fractions(b, delta):
    """Prescribed (alpha, delta) fractions in single precision, as in the generator."""
    if abs(delta - 0.2) < 1e-9:
        fa = f32(f32(0.7) - f32(0.1) * f32(b - 1)); fd = f32(0.2)
    else:
        fa = f32(f32(0.8) - f32(0.1) * f32(b - 1)); fd = f32(0.1)
    return fa, fd


def read_pdb(path):
    names, xyz = [], []
    for line in path.read_text().splitlines():
        if line.startswith("ATOM"):
            t = line.split()
            names.append(t[1]); xyz.append([float(t[2]), float(t[3]), float(t[4])])
    return np.array(names), np.array(xyz)


def check(root, delta):
    rows, ok = [], True
    for struc in ("core", "mix"):
        for f in sorted((Path(root) / struc).glob("config_ini_3D*.pdb")):
            b, s = map(int, NAME.search(f.name).groups())
            names, xyz = read_pdb(f)
            n = len(names)
            dup = n - len(np.unique(np.round(xyz, 3), axis=0))
            r = np.linalg.norm(xyz - CENTRE, axis=1) / A
            fa, fd = fractions(b, delta)
            na_e, nd_e = int(fa * n), int(fd * n)
            nb_e = n - na_e - nd_e
            cnt = {k: int((names == k).sum()) for k in ("C", "S", "Cd", "H")}
            # the single 'H' cell (if any) belongs to one of the three types
            counts_ok = (abs(cnt["C"] - na_e) + abs(cnt["S"] - nb_e) + abs(cnt["Cd"] - nd_e)) == cnt["H"]
            fill = n / (np.sqrt(2) * 4 * np.pi / 3 * is_size(s) ** 3)
            beta = names == "S"
            core_ok = ""
            if struc == "core":
                other = (~beta) & (names != "H")    # the H cell has no type; leave it out
                core_ok = "yes" if r[beta].max() <= r[other].min() + 1e-3 else "NO"   # PDB coordinates have 3 decimals
            tree = cKDTree(xyz)
            nb = tree.query_ball_point(xyz[beta], r=A * 1.01)
            nonb = np.array([sum(1 for j in lst if not beta[j]) for lst in nb])
            row = dict(struc=struc, Pbeta=10 * b, size=s, D=D_UM[s], n=n, n_lattice=lattice_count(s),
                       dup=dup, rmax=r.max(), is_size=is_size(s), fill=fill,
                       alpha=cnt["C"], beta=cnt["S"], delta=cnt["Cd"], H=cnt["H"],
                       beta_exp=nb_e, alpha_exp=na_e, delta_exp=nd_e, counts_ok=counts_ok,
                       core_inner=core_ok, beta_isolated=float(np.mean(nonb == 0)), nonbeta_nb=float(nonb.mean()))
            ok &= (n == EXPECTED_N[s] == row["n_lattice"]) and dup == 0 and r.max() < is_size(s) + 1e-3 and counts_ok \
                and core_ok in ("", "yes")
            rows.append(row)
    return rows, ok


def report(rows, ok, delta, root):
    out = [f"Structure checks: {root}   (generated with delta = {delta})", ""]
    out.append("Per file")
    out.append("struc Pb% s  D(um)     n  n_lat dup  rmax  is_size  fill  alpha  beta delta H   beta_exp ok  core_inner  "
               "beta_no_nonbeta_nb  mean_nonbeta_nb")
    for r in rows:
        out.append(f"{r['struc']:5s} {r['Pbeta']:3d} {r['size']:2d} {r['D']:5d} {r['n']:6d} {r['n_lattice']:6d} {r['dup']:3d} "
                   f"{r['rmax']:5.2f} {r['is_size']:7.2f} {r['fill']:5.3f} {r['alpha']:6d} {r['beta']:5d} {r['delta']:5d} "
                   f"{r['H']:1d} {r['beta_exp']:9d} {'ok' if r['counts_ok'] else 'NO':3s} {r['core_inner']:>9s}  "
                   f"{r['beta_isolated']:18.3f} {r['nonbeta_nb']:16.2f}")
    out.append("")
    out.append("Summary by size (identical for every P_beta and both architectures)")
    out.append(" size  D(um)     n  max radius (cell diameters)  is_size  sphere fill")
    seen = set()
    for r in rows:
        if r["size"] in seen:
            continue
        seen.add(r["size"])
        same = {(x["n"], round(x["rmax"], 4)) for x in rows if x["size"] == r["size"]}
        out.append(f"{r['size']:5d} {r['D']:6d} {r['n']:6d} {r['rmax']:27.2f} {r['is_size']:8.2f} {r['fill']:12.3f}"
                   + ("" if len(same) == 1 else "   (DIFFERS between files)"))
    out.append("")
    out.append("Fraction of beta cells with no non-beta neighbour, mean over P_beta (by D)")
    out.append(" D(um)   core    mix")
    for s in sorted(seen):
        v = {st: np.mean([x["beta_isolated"] for x in rows if x["size"] == s and x["struc"] == st]) for st in ("core", "mix")}
        out.append(f"{D_UM[s]:6d} {v['core']:6.3f} {v['mix']:6.3f}")
    out.append("")
    out.append(f"files: {len(rows)}   all checks passed: {ok}")
    return "\n".join(out)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root")
    ap.add_argument("--delta", type=float, default=0.2)
    ap.add_argument("-o", "--out")
    a = ap.parse_args()
    rows, ok = check(a.root, a.delta)
    if not rows:
        sys.exit(f"no PDB files under {a.root}/core or {a.root}/mix")
    text = report(rows, ok, a.delta, a.root)
    print(text)
    if a.out:
        Path(a.out).write_text(text + "\n")


if __name__ == "__main__":
    main()
