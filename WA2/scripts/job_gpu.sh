#!/bin/bash
# ACA WA2 - template job for the GPU partition (4x NVIDIA A100 40GB + 2x AMD EPYC 7742 per node).
#SBATCH --job-name=wa2-gpu
#SBATCH --account=CHANGE_ME_g            # GPU accounts end in 'g' (see `billing`)
#SBATCH --partition=normal-a100-40       # dev-a100-40 for short tests
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --gpus=1                         # up to 4 for the SYCL scalability study (2.6)
#SBATCH --cpus-per-task=32
#SBATCH --time=00:30:00
#SBATCH --output=wa2-gpu-%j.out
##SBATCH --exclusive                     # uncomment for scalability / energy runs

ml CUDA/11.8.0                           # check with `ml spider CUDA` / `ml spider NVHPC`

cd "${SLURM_SUBMIT_DIR}"
nvidia-smi -q > platform-${SLURM_JOB_ID}.txt
lscpu >> platform-${SLURM_JOB_ID}.txt

DATA=varcoef_ram.dat    # e.g. tools/gen_input.py varcoef --size-mib 1024 --measure memory --round up -o varcoef_ram.dat

# your CUDA / SYCL binary; validate it against the CPU reference with tools/compare
srun ./my_jacobi_cuda -i 100 ${DATA} result_gpu.bin

# profiling example (Task 2.7):
# srun nsys profile -o wa2-%q{SLURM_JOB_ID} ./my_jacobi_cuda -i 20 ${DATA} result_gpu.bin
