# Reproduction guide and technical details

Code and data for the PLOS ONE manuscript on Core and Mix islets of Langerhans.
The islets are modelled as networks of coupled oscillators.

**Authors:** Yun Gyeom Kim, Junghyo Jo, Taegeun Song
<!-- TODO: final author list, affiliations, title, DOI and citation once accepted. -->

## Project overview

### Question

How does islet architecture shape insulin output and glucose responsiveness? The study
varies four factors:

* the arrangement of the cell types: **Core**, with β cells in the centre, or **Mix**,
  with all types placed at random;
* the β-cell fraction **P_β**;
* the islet size **D**;
* the intercellular coupling **K**.

### Model

The cells are α, β and δ cells on a 3D hexagonal close-packed lattice. Each cell is a
glucose-modulated Stuart–Landau-like oscillator with amplitude A and phase θ:

    dA_i/dt = (λ_σ(g) − A_i²) A_i + Σ_j K_σjσi A_j cos(θ_j − θ_i)
    dθ_i/dt = ω_i + Γ_σ (g − g0) cos θ_i + Σ_j K_σjσi (A_j / A_i) sin(θ_j − θ_i)
    dg/dt   = g0 · H_α − g · H_β + I(t)

The terms are:

* The sums run over the 12 nearest lattice neighbours.
* λ_α = ½[1 + tanh((7 − g)/5)], λ_β = ½[1 + tanh((g − 7)/5)],
  λ_δ = 0.5 · ½[1 + tanh((g − 5)/5)].
* Γ_α = +0.1 and Γ_β = Γ_δ = −0.1.
* ω_i ≈ 2π/5 + 0.1·N(0,1).

The model has **heterotypic coupling only** (A_σσ = 0). Cells of the same type are not
coupled. Between types, the sign follows the code `info = 121212`:

* α → β, α → δ and β → δ coupling is +K;
* β → α, δ → α and δ → β coupling is −K.

The definitions follow the code. **H_α and H_β are per-cell means**:
H_β = (1/N_β) Σ_{j∈β} A_j (1 + cos θ_j), and H_α is the same over the α cells.
The glucose input I(t) has three 60-min phases: normal 0–60 min, high 60–120 min and
low 140–200 min (table below).

### Metrics

* **IPA** (insulin-production activity) is the time average of H_β over one 60-min
  glucose phase.
* **IS** = [IPA(G/G0 = 2.2) − IPA(G/G0 = 0.6)] / 1.6,
  taken within one simulation batch.

### Analysis

A random forest is fitted to IPA (features G, K, P_β, D, S) and to IS (features K, P_β,
D, S), with 5-fold held-out validation. The analysis then measures:

* single-feature importance, as held-out permutation importance;
* pair interactions, as the Friedman-H additive excess of the 2-way partial dependence.

A pair counts as an interaction only if it exceeds **2 × the noise floor**. The noise
floor is the same statistic computed on a target shuffled across configurations.

### Main findings

1. Core islets have higher IPA and higher IS than Mix islets.
2. Core IPA rises steadily with glucose. Mix IPA peaks at G/G0 = 1.8 and then declines.
3. At K = 0 the Core–Mix difference vanishes.
4. Of all feature pairs, only **K × S** (coupling × structure) exceeds 2 × the noise
   floor (IPA, 2.06×). No pair does so for IS.

### Repository map

```
fortran/      simulation code used for the paper (Core, Mix) + inputs     -> Steps 1-2
scripts/      batch runner, IPA extraction, conditions/IS, Fig 2, RF, Fig 3
data/         processed results: IPA per run, per condition, IS
results/      figures and RF tables produced by the scripts
structure/    3D structure generator, PDB files (Fig 1a), drawing script, structure checks
original/     the original code as received (verbatim) + table of every change made
```

The repository reproduces the full chain:

```
Fortran simulation (one batch)  ->  IPA per run  ->  IPA per condition, IS  ->  Fig 2
                                                                         ->  random forest  ->  Fig 3
```

Every figure can be rebuilt from the processed tables in `data/` in a few minutes,
without running the simulation (Quick start).
The simulation code is also included. With it, every run can be repeated from scratch.

**License:** to be decided (see the end of this file).  **DOI:** TODO.

---

## Contents

| Path | What it is |
|---|---|
| `fortran/islet_core.f90` | Simulation, **Core** islet: beta cells fill the sites nearest the centre; alpha and delta cells are placed at random in the outer shell |
| `fortran/islet_mix.f90` | Simulation, **Mix** islet: all three cell types are placed at random |
| `fortran/system_parameter.in`, `structure_parameter.in` | Input parameters, unchanged from the paper runs |
| `scripts/run_batch.sh` | Step 1. Compiles the code and runs one batch (648 runs) in parallel |
| `scripts/extract_ipa.py` | Step 2. Reads the `hb_*` time series and writes the IPA of each run and glucose phase |
| `scripts/build_conditions.py` | Step 3. Averages the realizations and computes the IS |
| `scripts/make_fig2.py` | Step 4. Fig 2 |
| `scripts/rf_analysis.py` | Step 5. Random-forest importance and interaction analysis |
| `scripts/make_fig3.py` | Step 6. Fig 3 |
| `scripts/interaction_importance.py` | RF / Friedman-H module by Taegeun Song, included unmodified (also in `original/`) |
| `scripts/common.py` | How the loop indices map to the physical parameters, the batch table, and the phase windows |
| `scripts/provenance/make_published_table.py` | Documents how `data/ipa_realizations.csv` was built from the authors' intermediate files. You do not need it to reproduce the paper |
| `data/ipa_realizations.csv` | IPA of every run: 4 batches × 648 runs × 3 glucose phases = 7,776 rows |
| `data/conditions.csv` | Realization mean and SD for each condition (2,592 rows) |
| `data/is_values.csv` | IS for each configuration (216 rows) |
| `results/` | Figures and RF tables written by the scripts |
| `structure/` | 3D structure generator, 72 PDB files, Fig 1(a)-style drawing script, structure checks (see `structure/README.md`) |
| `original/` | The original code as received, verbatim, and a table of all changes to the paper-run code (see `original/README.md`) |

## Requirements

* Simulation: a Fortran 90 compiler. Tested with GNU Fortran 9.4 (Ubuntu 20.04) and bash.
* Analysis: Python ≥ 3.8 with numpy, pandas, scikit-learn and matplotlib.
  Tested with Python 3.12.3, numpy 2.5.1, pandas 3.0.5, scikit-learn 1.9.0 and
  matplotlib 3.11.1. The extraction step also ran with Python 3.8, numpy 1.19 and pandas 1.2.

## Quick start (figures from the processed data)

```bash
python scripts/build_conditions.py   # data/ipa_realizations.csv -> data/conditions.csv, data/is_values.csv
python scripts/make_fig2.py          # results/fig2.png (+ fig2_values.csv)
python scripts/rf_analysis.py        # results/rf/*.csv          (~15 min, uses all cores)
python scripts/make_fig3.py          # results/fig3.png
```

Expected key numbers:

| Quantity | Value |
|---|---|
| Fig 2(a), K = 0.4, Core IPA at G/G0 = 0.6 … 2.6 | 0.51, 0.64, 0.75, 0.85, 0.93, 1.01 |
| Fig 2(a), K = 0.4, Mix IPA at G/G0 = 0.6 … 2.6 | 0.47, 0.55, 0.64, 0.67, 0.65, 0.59 (peak at 1.8) |
| RF held-out R², IPA (1,296 rows) / IS (216 rows) | 0.999 / 0.990 |
| Pairs above 2 × noise floor | IPA: K × S only (2.06 × floor). IS: none |

The RF numbers can differ in the third decimal from the values in the manuscript working files. The reason is the last-bit differences in the IPA that come from rebuilding the realizations. Example: K × S is 2.061 × floor here and 2.064 × floor there. The conclusions do not change.

## Full reproduction (from the simulation)

### Step 1 — run a batch

```bash
scripts/run_batch.sh 0.6 1204          # X = 0.6, batch id 1204, 5 jobs in parallel
```

A batch is 2 structures × 6 sizes × 6 P_β × 3 K × 3 realizations = **648 runs**.
The script splits them into 36 jobs, one for each structure × size × realization.
Each job runs in its own folder, because the program writes scratch files into its
working directory.

**Timing:** on an 8-core Linux server with 5 jobs in parallel, one batch (648 runs)
took about **81 min**. Raw output is about 470 MB per batch.

To check the setup first, run a short test: `SIZES=1 NTS=1 scripts/run_batch.sh 0.6 test 2`
(36 runs, about 2 min).

Run the program directly with:

```
islet_core X NT SEED [SIZE] [BRATIO] [EPP]      (islet_mix: same arguments)
```

`SIZE`, `BRATIO` and `EPP` limit the run to one value of that loop. If you omit them
or give 0, the program runs all values. The header of each `.f90` file explains every
argument.

### Step 2 — extract the IPA

```bash
python scripts/extract_ipa.py runs/1204 --batch 1204 --x 0.6 -o results/ipa_1204.csv
```

### Step 3 onwards

To use your own runs, give their tables to `build_conditions.py`.
Each new table needs its own batch id:

```bash
python scripts/build_conditions.py --realizations results/ipa_0123.csv \
       --extra results/ipa_1011.csv results/ipa_1204.csv results/ipa_1009.csv --outdir results/mydata
```

`make_fig2.py` and `rf_analysis.py` read from `data/`. To use your own tables, copy
them into `data/` or edit `DATA` in `scripts/common.py`.

---

## Model in brief

Each cell *j* is a Stuart–Landau oscillator with amplitude *A_j* and phase *θ_j*.
The cells sit on a 3D hexagonal close-packed lattice, and nearest neighbours are
coupled with strength *K* (`eps`). Each cell type has a glucose-dependent growth rate
λ(*g*). Glucose follows

    dg/dt = g0 · h_α − g · h_β + I(t)

Here, h_β = (1/N_β) Σ_{j∈β} A_j (1 + cos θ_j) is the **mean output per β cell**
(`hb` in the code). h_α is the same quantity for the α cells.

**Glucose input I(t).** Total time is 200 min, with dt = 0.01 min (`itmax` = 20000).
The model uses `extamp` = 2 and g0 = 7 mM.

| Time (min) | I(t) | Phase | G/G0 label | Output rows averaged for the IPA |
|---|---|---|---|---|
| 0–60 | 0 | normal | 1.0 | 0–599 |
| 60–120 | +X · extamp · g0 | high | 1 + 2X | 600–1199 |
| 120–140 | 0 | recovery | – | not used |
| 140–200 | −0.2 · extamp · g0 | low | 0.6 (4.2 mM) | 1400–1999 |

The G/G0 labels are nominal. A label is the fixed point of the glucose equation when
the hormone output equals 1. The g(t) that the model produces is written to `gt_*.dat`.

**IPA** (insulin-production activity) is the time average of h_β over one phase.

**IS** = [IPA(2.2) − IPA(0.6)] / 1.6. Both IPA values are taken from batch 1204
(X = 0.6). In that batch, the high phase is G/G0 = 2.2 and the low phase is 0.6.

### Parameters

| Loop variable (code) | Values | Paper parameter |
|---|---|---|
| `sml_size` | 1 … 6 | islet radius `is_size` = 1.2·sml_size + 3.8 (lattice units). Number of cells N = 743, 1405, 2403, 3791, 5563, 7881. Diameter **D** = 2·5.1 µm·N^(1/3) = 92, 114, 137, 159, 181, 203 µm |
| `b_ratio` | 1 … 6 | β-cell percentage **P_β** = 10 … 60 %. δ is fixed at 10 %, so α = 90 % − P_β |
| `epp` | 1 … 3 | coupling **K** = 0.4·(epp − 1) = 0, 0.4, 0.8 |
| `kdd` | 1 | self (diagonal) coupling = 0 |
| `nt` (argument NT) | 1 … 3 | realization index |
| structure | Core / Mix | **S** (0 / 1 in the RF) |
| `system_parameter.in` | `itmax` 20000, `dt` 0.01, `extamp` 2, `Gam` 0.1, `info` 121212 (coupling-sign code), `g0` 7, `lamdmax` 0.5 | |

### Batches

The four batches are identical except for the high-phase coefficient X
(code line `st = X*extamp*g0` for 60 < t ≤ 120).

| Batch id | X | High phase G/G0 | Nominal G_high (mM) | Source of the published values |
|---|---|---|---|---|
| 0123 | 0.2 | 1.4 | 9.8 | Re-run on 2026-10-27 (separate jobs, distinct seeds) |
| 1011 | 0.4 | 1.8 | 12.6 | Original paper runs |
| 1204 | 0.6 | 2.2 | 15.4 | Original paper runs (also the source of the IS) |
| 1009 | 0.8 | 2.6 | 18.2 | Original paper runs |

Every batch also contains the normal (1.0) and low (0.6) phases. In Fig 2 and in the
RF, the values for these two levels are first averaged over the four batches for each
configuration. Each glucose level therefore gives one value per configuration.

### Output files of the simulation

The Core program writes to `DATA_c/`. The Mix program writes to `DATA/`.
File names have the form `<var>_islet<nt>_<b_ratio>_<sml_size>_<epp−1>_0_ext2g0.dat`.
Each file has 2 columns (time in min, value) and 2000 lines, one line every 10 steps.

| `<var>` | Content |
|---|---|
| `hb` | mean output per β cell (gives the IPA) |
| `ha`, `hd` | the same for α and δ cells |
| `gt` | g(t)/g0 |
| `order_islet…_nonbc_{alpha,beta,delta}` | phase order parameter of each cell type |

Each run writes 7 files (about 0.7 MB in total). The raw time series of the four
batches (2,592 runs) are **not** in this repository because of their size.
<!-- TODO: decide whether to deposit the raw time series (~1.9 GB) on Zenodo; give the DOI here. -->
The scripts in Steps 1–2 regenerate them.

### Processed data (`data/`)

`ipa_realizations.csv`: one row for each run and glucose phase.

| Column | Meaning |
|---|---|
| `batch`, `X` | batch id and its high-phase coefficient |
| `struc` | `core` or `mix` |
| `K`, `Pbeta`, `D` | coupling, β-cell percentage (%), diameter label (µm) |
| `size_idx`, `ncell` | `sml_size` and the number of cells |
| `nt` | realization (1–3) |
| `phase`, `G` | `normal` / `high` / `low` and the nominal G/G0 label |
| `IPA` | mean of `hb` over the phase window |
| `source` | `rerun_2026-10-27` or `legacy_table_cumulative_mean_undone` (see Provenance) |

`conditions.csv` has the same keys without `nt`. Its columns `IPA` (mean of the
3 realizations), `IPA_sd` and `n_real` replace the per-run IPA.
`is_values.csv` gives `struc, K, Pbeta, D, IS`.

### Provenance and data corrections

* **Batches 1009 / 1011 / 1204.** These values were taken from the authors' legacy
  analysis tables. In those tables, the row for realization *k* holds the **cumulative
  mean** of realizations 1…*k*, not the value of realization *k*. The individual
  realizations were recovered as r1 = v1, r2 = 2v2 − v1, r3 = 3v3 − 2v2.
  The realization means in `conditions.csv` equal the legacy *k* = 3 values to
  within 1e-15.
* **Batch 0123.** The legacy Mix table of this batch also contained Core rows, and its
  low phase was not run with the paper protocol. The whole batch was therefore re-run
  on 2026-10-27 with the archived code. The protocol had a low phase of −0.2 and 3
  realizations, each run as a separate job with its own seed. The values were then
  extracted with `extract_ipa.py`.
* The largest diameter label was written as 202 µm in the legacy tables.
  2·5.1·7881^(1/3) = 203.0, so this repository uses 203. The RF uses D only for its
  order, so this label change gives identical results.
* `scripts/provenance/make_published_table.py` contains the exact procedure.

## Equivalence test (public code vs. paper-run code)

The public Fortran files differ from the archived paper-run files in four ways:

1. English comments were added.
2. X, NT and SEED are read from the command line.
3. The loops can be limited to one value.
4. The seed can be fixed.

Every equation and constant is the same. A `diff` with all comment lines removed shows
only these changes.

The archived code seeds its random generators from the clock (1/100 s). For the test,
a copy of the archived code was made with two changes:

* The loops were limited to one condition: `sml_size` = 1, `b_ratio` = 1, `nt` = 1, `epp` = 2 (K = 0.4).
* The clock seed line was replaced by one fixed-seed line: `seed = 0; seed(1) = 12345`.

The test copy and the public code (`islet_core 0.6 1 12345 1 1 2`) were each compiled
with gfortran 9.4 and run. Results:

| Test | -O0 | -O2 |
|---|---|---|
| Core, X = 0.6: all 7 output files and `fort.46` | byte-identical | byte-identical |
| Mix, X = 0.6 | byte-identical | byte-identical |
| Mix, X = 0.2 (single-precision literal vs. argument) | byte-identical | byte-identical |
| Public code run twice, same seed | byte-identical | byte-identical |

Notes:

* **X is read as a default (single-precision) REAL.** The archived code used a
  single-precision literal such as `0.6`. That literal becomes 0.6000000238… when it
  is promoted to double precision, and reading X as single precision gives the same
  double value. If X is read as a double, the result is not byte-identical.
* **Compiler flags change the trajectories.** The oscillator network is sensitive to
  floating-point rounding. The same seed therefore gives different trajectories under
  -O0 and -O2. In the test run, the IPA of a single run changed by up to 2.7 %
  (Mix, low phase). To get identical output, use the same compiler and the same flags.
* The original paper runs used clock seeds, so those individual runs cannot be
  repeated exactly. The published values are realization means. With fixed seeds,
  new runs are reproducible.

## License

<!-- TODO: choose a license (e.g. MIT or BSD-3-Clause for the code, CC BY 4.0 for the data)
     and confirm the reuse terms of the structure-generation code by D.-T. Hoang (APCTP). -->
To be decided.

## Citation

<!-- TODO: manuscript citation and Zenodo DOI. -->
