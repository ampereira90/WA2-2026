#!/bin/bash
# ACA WA2 - template job for the A64FX (ARM) partition.
#SBATCH --job-name=wa2-arm
#SBATCH --account=CHANGE_ME_a            # ARM accounts end in 'a' (see `billing`)
#SBATCH --partition=normal-arm           # dev-arm only for compilation / short tests
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --output=wa2-arm-%j.out

# Load the same toolchain used to build (check names with `ml spider GCC` on an ARM node)
ml GCC/13.3.0

cd "${SLURM_SUBMIT_DIR}"
BIN=bin/aarch64-gnu                      # or bin/aarch64-fujitsu, ...
VARIANT=poisson7                         # poisson7 | box27 | varcoef
DATA=${VARIANT}_256.dat                   # created beforehand, e.g.
                                         #   tools/gen_input.py ${VARIANT} --dims 256 256 256 -o ${VARIANT}_256.dat
mkdir -p results

# record the platform (Task 1.4: reproducibility)
lscpu > platform-${SLURM_JOB_ID}.txt
numactl -H >> platform-${SLURM_JOB_ID}.txt 2>/dev/null

export OMP_PLACES=cores
export OMP_PROC_BIND=close
for T in 1 12 24 48; do
    export OMP_NUM_THREADS=$T
    for rep in 1 2 3 4 5; do
        srun --cpu-bind=none ${BIN}/${VARIANT} -i 20 ${DATA} results/${VARIANT}_T${T}.bin | grep '^CSV' | sed "s/^/T=${T},/"
    done
done
