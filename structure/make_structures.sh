#!/usr/bin/env bash
# =====================================================================================
# Generate the 72 islet structures (2 architectures x 6 P_beta x 6 sizes) as PDB files.
#
#   structure/make_structures.sh [SEED] [DELTA] [OUTDIR]
#
#   SEED   > 0 deterministic (default 1); 0 = clock seed as in the archived run
#   DELTA  0.2 (default) = fractions of the archived generator (published PDB files)
#          0.1           = fractions used in the simulations (delta 10 %)
#   OUTDIR default structure/pdb_regenerated  -> OUTDIR/{core,mix}/config_ini_3D<b>,<s>.pdb
#
# Takes a few seconds per architecture.  Requires gfortran.
# =====================================================================================
set -euo pipefail
SEED=${1:-1}
DELTA=${2:-0.2}
HERE=$(cd "$(dirname "$0")" && pwd)
OUT=${3:-$HERE/pdb_regenerated}
FC=${FC:-gfortran}
FFLAGS=${FFLAGS:--O2}
WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

$FC $FFLAGS -o "$WORK/generate_structures" "$HERE/generate_structures.f90"
for type in core mix; do
  mkdir -p "$WORK/$type" "$OUT/$type"
  cp "$HERE/system_parameter.in" "$HERE/structure_parameter.in" "$WORK/$type/"
  (cd "$WORK/$type" && "$WORK/generate_structures" "$type" "$SEED" "$DELTA" > run.log 2>&1)
  cp "$WORK/$type"/config_ini_3D/*.pdb "$OUT/$type/"
  echo "$type: $(ls "$OUT/$type" | wc -l) PDB files -> $OUT/$type"
done
echo "check with: python structure/check_structures.py $OUT --delta $DELTA"
