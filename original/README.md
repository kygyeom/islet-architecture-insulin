# Original code (verbatim)

This folder keeps the starting code exactly as the authors received it. Nothing in it
was edited. The code used for the paper is in `../fortran/`. The table below lists every
change from this folder to that code.

| File | What it is | Author |
|---|---|---|
| `multi_islets_bio.f90` | The original islet model. One HCP islet of α, β and δ Stuart–Landau oscillators, coupled to a glucose variable. The cell types are arranged by a Monte-Carlo cell-sorting step | Taegeun Song. The HCP lattice and Monte-Carlo routines are by Danh-Tai Hoang (APCTP); see the comment block in the file |
| `system_parameter.in`, `structure_parameter.in` | Input files that came with it | same |
| `interaction_importance.py` | Random-forest permutation importance and Friedman-H interaction analysis. The same file is in `../scripts/` and is used without change | Taegeun Song |

## How the original file was identified

The candidates were the following files in the authors' working folder:

* the two top-level files `multi_islets_bio.f90` and `multi_islets_bio_core.f90`;
* the earliest run folders, for example `0825_G_confiem/1/multi_islets_bio.f90`;
* the later variants `multi_core_mixed_str.f90` and `multi_core_ad_mixed.f90`.

The top-level `multi_islets_bio.f90` was chosen for five reasons:

1. **It is the simplest file.** It simulates one islet with fixed fractions
   (α 0.3, β 0.6, δ 0.1) and a fixed radius (`is_size = 3.1`). It has none of the
   loops over size, β fraction, coupling or realization. Every other candidate has
   some of these loops.
2. **It has no lab changes to the cell placement.** It has no Core routine: no
   `init_structure_core` and no distance ranking. It also has no random shuffle
   (Fisher–Yates). In the derived files the shuffle has a Korean comment
   (`배열 요소 랜덤하게 섞기`), which marks it as a lab addition. The earliest run
   folder (`0825_G_confiem/1`, 2023-08-27) already has the shuffle.
3. **The Monte-Carlo cell sorting is still active.** `CALL load_config_ini`,
   `equi_lattice1` and `write_config_3D` run at lines 795–804. In every derived file
   these calls are commented out.
4. **Only author comments by others.** The only author comments are
   "edited Taegeun SONG" and the D.-T. Hoang block.
5. **It is the "base code" in git.** The authors' local git history stores it, with
   the top-level `multi_islets_bio_core.f90`, as the single commit "Add base code".
   The top-level `multi_islets_bio_core.f90` is not the original: it already has the
   size and β loops, the 3-phase input and K loops.

**Uncertainty.** The file date (2023-10-09) is later than the earliest derived run
folders (2023-08-27). The file was probably saved again or copied after the derived
versions were made. The identification therefore rests on the content, not on the date.
We have no copy of the file as first received, for example as an e-mail attachment.

## What we changed and why (original → paper-run code in `../fortran/`)

Line numbers refer to two files:

* **O** = `original/multi_islets_bio.f90`
* **P** = the archived paper-run code (`1011_result/생성/{core_0.6,mix_0.6}`), which
  `../fortran/islet_{core,mix}.f90` reproduce byte for byte.

PC = Core paper-run file, PM = Mix paper-run file.

| # | What changed | Original (O) | Paper-run code (P) | Why |
|---|---|---|---|---|
| 1 | Islet size | `is_size(is)=3.1` (O 59). Array limit `mxdim=2500` (O 6) | Loop `do sml_size = 1, 6` and `is_size(is)=sml_size*1.2 + 3.8` (PM 42, 58). `mxdim=100000` | Six islet sizes: 743 to 7881 cells, D = 92 to 203 µm |
| 2 | β fraction | Fixed `s_ratio = 0.3 / 0.6 / 0.1` (O 53–55) | Loop `do b_ratio = 1, 6`. `s_ratio(is,1)=0.8 - 0.1*(b_ratio-1)`, `s_ratio(is,2)=0.1 + 0.1*(b_ratio-1)`, `s_ratio(is,3)=0.1` (PM 43, 51–53) | P_β from 10 to 60 % with δ fixed at 10 % |
| 3 | Mix arrangement | Types are set in lattice order: all α first, then β, then δ (O 556–560). Then the Monte-Carlo cell sorting (`equi_lattice1`, `n_equi1` = 500000, O 795–804) rearranges them | Random Fisher–Yates shuffle of the type array (PM 672: `do i = icount, 2, -1 … random_number`). The Monte-Carlo calls are commented out (PM 916–925) | A well-mixed reference structure whose arrangement does not depend on sorting energies |
| 4 | Core arrangement | No Core option | The n_β lattice sites nearest the centre become β. The other sites take the shuffled α/δ types (PC 671–677 shuffle, PC 737–790 distance ranking). Earlier variants had the same algorithm as a separate routine `init_structure_core` (`multi_core_mixed_str.f90`, `multi_core_ad_mixed.f90`, and `../structure/generate_structures.f90`). In the paper-run Core file it is inside `init_structure` | The Core architecture: β cells form the inner core, α and δ cells the mantle |
| 5 | Glucose input I(t) | One pulse, `if(time.gt.300..and.time.le.360.) st=extamp*g0` (O 111). `itmax = 7000`, dt = 0.01, so the run lasts 70 min | Three phases: `st= X*extamp*g0` for 60–120 min and `st=-0.2*extamp*g0` for 140–200 min (PM 187–188). `itmax = 20000` (200 min). X = 0.2, 0.4, 0.6 or 0.8 per batch | Normal, high and low glucose phases. Their averages give the IPA, and the high and low phases give the IS. ⚠ In O the pulse starts at 300 min but the run ends at 70 min, so O as delivered never applies the pulse. This is an observation only; we do not know how the original author used the file |
| 6 | Coupling K | `eps` and `kdia` are read from `system_parameter.in` (O 201–204; the file gives eps = 0) | Loop `do epp = 1,3`, `eps = 0.4*(epp-1)`. `kdia = 0` from `do kdd = 1,1`. The values in the input file are now skipped (`read(900,*) chtmp!, eps`). `coupmat` is called in the main loop (PM 45, 68–87) | Coupling K = 0, 0.4, 0.8 |
| 7 | Realizations and random seed | One run. `ran_seed` is called after the structure is built (O 68) | Loop `do nt = 1,3`. `ran_seed` is called before `get_structure` in every iteration, so each structure has its own clock seed (PM 44, 62). Paper runs: one program call per realization; the 2026-10-27 re-run gave each job a distinct seed. `../fortran/` takes SEED as an argument | Independent realizations of placement and initial state |
| 8 | ran1 initialisation | `idum=max(idum,1)` (O 456) | `idum=max(-idum,1)` (PC 553, PM 558) | ⚠ **Unsure why this was changed.** With a positive seed, the paper-run code starts ran1's shuffle table from 1. So in the first loop iteration of each program call, the oscillator initial values do not depend on the seed. Realizations still differ through the cell placement, and later iterations continue from that iteration's seed. Kept unchanged so that the paper results are reproduced |
| 9 | δ-cell output | Only `ha` and `hb` are computed (O 151–156) | `hd` is added (PC 228) and written to `hd_*.dat`. The glucose equation is unchanged: `dg=(ha*g0-hb*g)+st` (O 159, PM 237) | Diagnostic output only |
| 10 | Output files | `DATA/{gt,ha,hb,order…}_islet1_ext2g0.dat` (O 40–49). `gt` holds g, written every step (O 163) | `DATA` (Mix) or `DATA_c` (Core) `/<var>_islet<nt>_<b_ratio>_<sml_size>_<epp-1>_<kdd-1>_ext2g0.dat`, every 10 steps. `gt` holds g/g0 (PM 151–164, 240). Structure writing is switched off (`info_structure_write=1`) | One file set per run. File names carry the parameters |
| 11 | Model equations and constants | — | Unchanged: oscillator equations, λ(g), pinning Γ, coupling sign code `info = 121212`, ω, initial-condition ranges, hormone definitions h = (1/N) Σ A(1+cos θ) | — |

## Things we are unsure of

* Whether `multi_islets_bio.f90` is byte-identical to the file as first received. It is
  the oldest content we have; see the uncertainty note above.
* Why ran1's initialisation was changed (row 8).
* Whether the original author meant the 300–360 min pulse for a longer `itmax` (row 5).
