#!/usr/bin/env bash
# =====================================================================================
# Step 1 - run one simulation batch (Core + Mix, all 648 runs) in parallel.
#
#   scripts/run_batch.sh X BATCH_ID [NPAR] [SEED_BASE] [WORKDIR]
#
#   X          high-phase input coefficient (0.2, 0.4, 0.6 or 0.8 in the paper)
#   BATCH_ID   label used for the output folder and the extracted table (e.g. 1204)
#   NPAR       number of jobs run at the same time            (default 5)
#   SEED_BASE  > 0: deterministic seeds, job j (0..35) uses SEED_BASE + 100*j
#                   (each job advances the seed by 1 per run, 18 runs per job)
#              0  : clock seeds as in the original runs (NOT reproducible)
#                                                               (default 1000)
#   WORKDIR    where the jobs are written                     (default runs/<BATCH_ID>)
#   Environment: FC (gfortran), FFLAGS (-O2); SIZES / NTS restrict the jobs for a quick
#   test, e.g.  SIZES=1 NTS=1 scripts/run_batch.sh 0.6 test   (2 jobs, 36 runs, ~2 min).
#
# One batch = 2 structures x 6 sizes x 6 P_beta x 3 K x 3 realizations = 648 runs.
# The work is split into 36 jobs (structure x size x realization); each job loops over
# P_beta (6) and K (3) inside one program call and writes into its own folder, because
# the Fortran code writes scratch files (tmp_fname, tmp_nt, tmp_ep, fort.46) into the
# working directory.  Largest islets are started first.
#
# Reference timing: 8-core Linux server, NPAR=5, gfortran 9.4: 648 runs in ~81 min.
# Raw output: ~13 MB per job folder (~470 MB per batch), 7 time-series files per run.
#
# Next step:  python scripts/extract_ipa.py runs/<BATCH_ID> --batch <BATCH_ID> --x X -o results/ipa_<BATCH_ID>.csv
# =====================================================================================
set -euo pipefail

X=${1:?usage: run_batch.sh X BATCH_ID [NPAR] [SEED_BASE] [WORKDIR]}
BATCH=${2:?usage: run_batch.sh X BATCH_ID [NPAR] [SEED_BASE] [WORKDIR]}
NPAR=${3:-5}
SEED_BASE=${4:-1000}
REPO=$(cd "$(dirname "$0")/.." && pwd)
WORK=${5:-$REPO/runs/$BATCH}
FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2}
SIZES=${SIZES:-6 5 4 3 2 1}      # override only for a quick test, e.g. SIZES=1 NTS=1
NTS=${NTS:-1 2 3}

mkdir -p "$WORK/bin"
echo "compiling with: $FC $FFLAGS"
$FC $FFLAGS -o "$WORK/bin/islet_core" "$REPO/fortran/islet_core.f90"
$FC $FFLAGS -o "$WORK/bin/islet_mix"  "$REPO/fortran/islet_mix.f90"

JOBS="$WORK/jobs.txt"; : > "$JOBS"
j=0
for size in $SIZES; do                  # largest islets first (longest jobs)
  for struc in core mix; do
    for nt in $NTS; do
      d="$WORK/${struc}_s${size}_nt${nt}"
      mkdir -p "$d"
      cp "$REPO/fortran/system_parameter.in" "$REPO/fortran/structure_parameter.in" "$d/"
      if [ "$SEED_BASE" -gt 0 ]; then seed=$((SEED_BASE + 100 * j)); else seed=0; fi
      # arguments: X NT SEED SIZE BRATIO(0 = all) EPP(0 = all)
      echo "cd '$d' && '$WORK/bin/islet_$struc' $X $nt $seed $size 0 0 > run.log 2>&1 && echo \"\$(date +%T) done $struc size=$size nt=$nt\"" >> "$JOBS"
      j=$((j + 1))
    done
  done
done

echo "batch $BATCH: X=$X, $j jobs (18 runs each), $NPAR in parallel, seed base $SEED_BASE -> $WORK"
date > "$WORK/started.txt"
tr '\n' '\0' < "$JOBS" | xargs -0 -P "$NPAR" -I{} bash -c '{}'
date > "$WORK/finished.txt"
echo "finished. Next: python scripts/extract_ipa.py '$WORK' --batch $BATCH --x $X -o results/ipa_$BATCH.csv"
