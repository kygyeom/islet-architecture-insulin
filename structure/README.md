# 3D islet structures

This folder generates, draws and checks the islet structures shown in Fig 1(a).

| File | What it is |
|---|---|
| `generate_structures.f90` | Commented copy of the archived generator (`1009_structure_make/multi_islets_bio_struc_make.f90`). It writes one PDB file per (P_β, size) |
| `make_structures.sh` | Compiles the generator and writes all 72 files: core and mix × 6 P_β × 6 sizes |
| `system_parameter.in`, `structure_parameter.in` | Inputs; identical to `../fortran/` |
| `pdb/core/`, `pdb/mix/` | The 72 archived PDB files used for the figures (17.8 MB) |
| `draw_islet_3d.py` | Draws one PDB file in the style of Fig 1(a): α red, β blue, δ green, one quarter cut away |
| `islets_3D_image.ipynb` | The original drawing notebook (kept as is) |
| `check_structures.py` | Structure checks (see below) |
| `check_output.txt` | Output of the checks for the archived files and for regenerated files with δ = 10 % |

File names are `config_ini_3D<b_ratio>,<sml_size>.pdb`. The first index is the β index
(P_β = 10·b_ratio %). The second is the size index (743 … 7881 cells, D = 92 … 203 µm).
Atom names: `C` = α, `S` = β, `Cd` = δ.

## ⚠ The δ fraction differs from the simulations

The archived generator, and therefore every file in `pdb/`, uses these fractions:

* δ = 20 %
* α = 70 % − 10 %·(b_ratio − 1)

All simulations in the paper used δ = 10 % and α = 80 % − 10 %·(b_ratio − 1).
In both cases β = 10 %·b_ratio. Positions and β placement do not depend on this choice.
In Core files, β cells occupy the same sites. Only the α : δ split of the remaining
cells differs.

To make structures with the simulation fractions, run
`make_structures.sh 1 0.1`.

## Regenerate

```bash
structure/make_structures.sh 1 0.2      # archived fractions -> structure/pdb_regenerated/{core,mix}
structure/make_structures.sh 1 0.1      # simulation fractions (delta 10 %)
python structure/check_structures.py structure/pdb_regenerated --delta 0.1
python structure/draw_islet_3d.py structure/pdb/core/config_ini_3D4,4.pdb -o results/core_P40_D159.png
```

Changes from the archived generator (calculation unchanged):

1. English comments.
2. TYPE (core/mix), SEED and DELTA are read from the command line. The archived file
   had `structure_type = 1` hard-coded and was edited by hand for the Core run.
3. The output folder is created once.
   The archived code ran `rm -r config_ini_3D` inside the per-islet routine. On Linux,
   that left only the last file; we checked this with gfortran 9.4. The archived files
   were made on Windows. There, the `rm` call presumably failed without effect, because all 36 files of each run exist.

The code that places the cells (`init_structure` for Mix, `init_structure_core` for
Core) is identical to the placement code of the simulations.

We regenerated the files with DELTA 0.2 and compared them with the archived ones:

| Architecture | Result |
|---|---|
| Core | 36/36 files have the same coordinates and the same β-cell set |
| Mix | 36/36 files have the same coordinates. The cell types differ, because Mix placement is random and the archived run used a clock seed |

## Checks (`check_structures.py`)

The script checks these things for every file:

* the cell count, against an independent count of the lattice;
* duplicate coordinates;
* the largest radius, against is_size;
* sphere fill;
* the α/β/δ counts, against the prescribed fractions;
* for Core, that β fills the innermost sites;
* the fraction of β cells with no non-β neighbour.

Result for both the archived and the regenerated files: **all 72 pass**.

| Check | Result |
|---|---|
| Cell number | 743, 1405, 2403, 3791, 5563, 7881 (same as the lattice count) |
| Duplicate coordinates | 0 |
| Max radius (cell diameters) | 5.00, 6.19, 7.39, 8.58, 9.80, 11.00 (is_size = 5.0 … 11.0) |
| Sphere fill | 0.995 – 1.006 of √2·(4π/3)·is_size³ |
| β count | = n − int(f_α n) − int(f_δ n), that is P_β·n up to rounding |

**Legacy detail.** In islets with ≥ 2000 cells, cell number 2000 is written with atom
name `H` instead of its type. This is the "i-consider" marker of the original
Monte-Carlo code. The check counts it separately.

The share of β cells with no non-β neighbour is shown below. With heterotypic coupling
only, such cells are not coupled to any other cell. Values are means over P_β, for the
archived files.

| D (µm) | Core | Mix |
|---|---|---|
| 92 | 0.38 | 0.003 |
| 114 | 0.47 | 0.002 |
| 137 | 0.52 | 0.002 |
| 159 | 0.57 | 0.001 |
| 181 | 0.63 | 0.001 |
| 203 | 0.66 | 0.001 |
