#!/bin/bash
# ACA WA2 - template job for the x86 (2x AMD EPYC 7742) partition.
#SBATCH --job-name=wa2-x86
#SBATCH --account=CHANGE_ME_x            # x86 accounts end in 'x' (see `billing`)
#SBATCH --partition=normal-x86           # dev-x86 only for compilation / short tests
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --exclusive
#SBATCH --time=00:30:00
#SBATCH --output=wa2-x86-%j.out

ml GCC/13.3.0                            # check with `ml spider GCC`

cd "${SLURM_SUBMIT_DIR}"
BIN=bin/x86_64-gnu
VARIANT=poisson7                         # poisson7 | box27 | varcoef
DATA=${VARIANT}_256.dat                   # created beforehand, e.g.
                                         #   tools/gen_input.py ${VARIANT} --dims 256 256 256 -o ${VARIANT}_256.dat
mkdir -p results

lscpu > platform-${SLURM_JOB_ID}.txt
numactl -H >> platform-${SLURM_JOB_ID}.txt 2>/dev/null

export OMP_PLACES=cores
export OMP_PROC_BIND=close
for T in 1 8 16 32 64 128; do
    export OMP_NUM_THREADS=$T
    for rep in 1 2 3 4 5; do
        srun --cpu-bind=none ${BIN}/${VARIANT} -i 20 ${DATA} results/${VARIANT}_T${T}.bin | grep '^CSV' | sed "s/^/T=${T},/"
    done
done
