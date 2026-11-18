#!/bin/bash
# ACA WA2 - validate an implementation against the reference.
#
#   scripts/validate.sh <variant> <reference_binary> <candidate command...>
#
# examples:
#   scripts/validate.sh poisson7 bin/x86_64-gnu/poisson7 ./my_poisson7
#   OMP_NUM_THREADS=8 scripts/validate.sh box27 bin/aarch64-gnu/box27 ./box27_omp
#   scripts/validate.sh varcoef bin/x86_64-gnu/varcoef srun -n 4 ./varcoef_mpi
#
# Generates small datasets with tools/gen_input.py (including odd and
# non-cubic grids), runs the reference and the candidate for 100 iterations
# on each, and compares the results with tools/compare. The candidate is
# called as:  <candidate command...> -i 100 INPUT OUTPUT
# Exit status 0 if all PASS. Needs python3 with numpy.
set -u
VARIANT=$1; REF=$2; shift 2
CAND=("$@")
HERE=$(cd "$(dirname "$0")/.." && pwd)
CMP=$(dirname "$REF")/compare
GEN="python3 $HERE/tools/gen_input.py"
TMP=${TMPDIR:-/tmp}/wa2-validate-$$
mkdir -p "$TMP"
STATUS=0
for DIMS in "32 32 32" "101 101 101" "128 128 128" "40 33 57"; do
    TAG=$(echo $DIMS | tr ' ' 'x')
    $GEN "$VARIANT" --dims $DIMS -o "$TMP/in_$TAG.dat" > /dev/null || exit 2
    "$REF" -i 100 "$TMP/in_$TAG.dat" "$TMP/ref_$TAG.bin" > /dev/null || exit 2
    if ! "${CAND[@]}" -i 100 "$TMP/in_$TAG.dat" "$TMP/cand_$TAG.bin" > "$TMP/run_$TAG.txt" 2>&1; then
        echo "$TAG: candidate failed"; tail -5 "$TMP/run_$TAG.txt"; STATUS=1; continue
    fi
    if "$CMP" "$TMP/ref_$TAG.bin" "$TMP/cand_$TAG.bin" > "$TMP/cmp_$TAG.txt"; then
        echo "$TAG: PASS  $(grep 'relative L2' "$TMP/cmp_$TAG.txt")"
    else
        echo "$TAG: FAIL"; cat "$TMP/cmp_$TAG.txt"; STATUS=1
    fi
done
rm -rf "$TMP"
exit $STATUS
