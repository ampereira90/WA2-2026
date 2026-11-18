# ACA — Work Assignment 2 — 3D Jacobi code base

This directory contains the reference implementations for Task 2 of Work Assignment 2. There are three variants of the 3D Jacobi iteration; your team works on the variant it was assigned:

| Variant | Directory | Problem |
|---|---|---|
| `poisson7` | `poisson7/` | Poisson equation, constant coefficients, 7-point stencil, analytic sine source |
| `box27` | `box27/` | Poisson equation, compact 27-point stencil, Gaussian point sources |
| `varcoef` | `varcoef/` | Steady heat conduction in a heterogeneous composite (4 materials), 7-point finite-volume stencil with harmonic-mean face conductivities |

All three solve a linear system on a regular grid of `nx × ny × nz` interior points, plus a one-point halo that holds the Dirichlet boundary values. Each Jacobi iteration computes a new iterate from the previous one only: `u_new(p) = F(u(neighbours of p))`. The header file of each variant (`poisson7.h`, `box27.h`, `varcoef.h`) states the exact update formula and the contents of its input dataset.

Each program reads the dataset to process from an **input file**, created with `tools/gen_input.py` for a required size in MiB, and writes the final solution to an **output file**. The reference codes are sequential and store values in single precision (`float`). They are written in portable C++11 and compile with GCC, the Fujitsu compilers (trad and clang mode), Intel icpx, NVIDIA nvc++ and clang, on both x86-64 and A64FX.

## Layout

```
common/      command line, timer, input file reader, result output, report (shared)
poisson7/    main.cpp  grid.cpp  jacobi.cpp  poisson7.h
box27/       main.cpp  grid.cpp  jacobi.cpp  stencil.cpp  box27.h
varcoef/     main.cpp  grid.cpp  jacobi.cpp  varcoef.h
tools/       gen_input.py   (creates input datasets)
             compare.cpp    (validates results)
             plot_result.py (renders a result as an image)
scripts/     Slurm templates for the ARM, x86 and A100 partitions, validate.sh
Makefile
```

## Building

The Makefile selects a toolchain with `COMPILER=` and puts the binaries in `bin/<arch>-<toolchain>/`. The same directory is visible from the x86 login nodes and from the A64FX compute nodes, so binaries for different architectures never overwrite each other.

```bash
make                         # GCC, native architecture, reference flags (-O2 -g)
make COMPILER=fujitsu        # FCC trad mode  (A64FX compute node)
make COMPILER=fujitsu-clang  # FCC -Nclang    (A64FX compute node)
make COMPILER=fujitsu-cross  # FCCpx on an x86 login node, binaries for A64FX
make COMPILER=intel          # icpx (x86)
make COMPILER=nvhpc          # nvc++ (x86 / GPU nodes)
make OPT="-O3 -fopenmp"      # extra flags appended to the reference flags
make info                    # print the compiler and the flags in use
```

On Deucalion the x86 compilers work on the login nodes. Native A64FX compilation must be done on an ARM compute node (`salloc -N1 -p dev-arm -A <account>a -t 1:00:00`) or with the Fujitsu cross compilers on a login node. Load the toolchain module first; `ml spider GCC`, `ml spider CUDA`, `ml spider NVHPC` and `ml spider intel` list the available versions. The x86 and ARM module trees are different, so check on each partition.

The **reference build** is `make` with the default GCC flags. Use it for every "reference" measurement, and record the module versions and the `make info` output in your report.

## Creating a dataset

`tools/gen_input.py` needs Python 3 with numpy (on Deucalion: `ml SciPy-bundle`, or `python3 -m pip install --user numpy`). It picks a cubic grid from a required size in MiB:

```bash
tools/gen_input.py varcoef --size-mib 1024 -o varcoef_1G.dat         # input file of ~1024 MiB
tools/gen_input.py varcoef --size-mib 1024 --measure memory --round up \
                   -o varcoef_ram.dat                                 # reference uses > 1 GiB of memory
tools/gen_input.py varcoef --size-mib 0.03125 --measure memory --round down \
                   -o varcoef_l1.dat                                  # reference data fit in 32 KiB
tools/gen_input.py varcoef --dims 40 33 57 -o small.dat               # explicit grid
tools/gen_input.py varcoef --size-mib 8 --dry-run                     # only print the chosen size
```

`--measure file` (the default) makes the size refer to the input file. `--measure memory` makes it refer to the memory used by the data structures of the *reference* code, which is what matters when you choose sizes that fit L1, the LLC or only RAM in Tasks 1.4 and 2.3. `--round nearest|down|up` picks the grid closest to the target, the largest one not above it, or the smallest one above it. The script prints the grid, the file size and the reference memory footprint.

The input file and the reference memory footprint differ because the file stores only the data that define the problem: the initial guess and boundary values, plus the source term (`box27`) or the material map and table (`varcoef`). The approximate reference footprints are 8 B/point (`poisson7`, plus row-pointer overhead), 16 B/point (`box27`) and 12 B/point (`varcoef`). If you change the data structures, compute the footprint of *your* version when you justify the sizes. The A64FX node has 32 GB of HBM2 shared by 4 CMGs (8 GB each). Generating a 1 GiB dataset takes a few seconds; generate the large ones once and reuse them.

The input format is a 56-byte header (variant, grid size, grid spacing `h`), a table of named fields with their type and byte offset, an optional material table, and then the fields. Each field covers the whole grid including the halo, `(nx+2)·(ny+2)·(nz+2)` values with `k` fastest, then `j`, then `i`. `common/common.h` documents it, and `InputFile::read_planes()` reads any range of `i`-planes of a field. A contiguous array with the same ordering can therefore be filled with a single call, and an MPI rank can read just its own slab.

## Running

```
bin/<arch>-<toolchain>/<variant> [-i ITERS] [-t TOL] [-p EVERY] INPUT OUTPUT
```

`INPUT` is a dataset created for the same variant; the program refuses a dataset made for another variant. `OUTPUT` receives the final solution. `-i` sets the number of Jacobi iterations (default 100). `-t` stops earlier when `||u_new − u||₂ < TOL` (default 0, i.e. never). `-p` prints the convergence history.

The program reports the load time (allocation and reading of the input), the solve time (iterations only), the time per iteration, the throughput in MLUP/s (million lattice-point updates per second) and a checksum. It also prints one line starting with `CSV,` that you can collect directly into the result tables of Task 2.1:

```
CSV,variant,nx,ny,nz,iterations,final_residual,t_load,t_solve,MLUP/s,checksum
```

The programs do not report GFLOP/s. Counting the floating-point operations and the memory traffic of the reference and of your implementations is part of Task 2.4. The reference codes run sequentially and need several seconds per iteration on the largest datasets, so use a small iteration count (e.g. `-i 10`) for reference timings on large grids and report the time per iteration.

## Validating your implementations

Every implementation (OpenMP, MPI, CUDA, SYCL, …) must produce the same result as the reference within rounding error. Run the reference and your version on the **same input file with the same iteration count**, and compare the two output files:

```bash
tools/gen_input.py poisson7 --dims 128 128 128 -o in_128.dat
bin/x86_64-gnu/poisson7 -i 100 in_128.dat ref_128.bin        # reference
./my_poisson7          -i 100 in_128.dat mine_128.bin       # your version
bin/x86_64-gnu/compare ref_128.bin mine_128.bin             # PASS / FAIL (exit code 0 / 1)
```

`compare` checks every point, `|u − u_ref| ≤ 1e-3·|u_ref| + 1e-5·max|u_ref|`, and the relative L2 error (≤ 1e-4). It reports the worst point as `(i,j,k)`, which helps to locate halo and loop-bound bugs. Results are not expected to be bit-identical: vectorisation, FMA contraction, single-precision intermediates and parallel reductions all change the rounding. **Validate with `-i 100`** on sizes that are not multiples of your vector width or block size. Rounding differences grow slowly with the iteration count. Missing iterations, wrong halos or off-by-one loop bounds produce errors far above the tolerance.

`scripts/validate.sh <variant> <reference binary> <your command…>` automates this. It generates four small datasets (32³, 101³, 128³ and 40×33×57), runs both programs and compares them, e.g. `scripts/validate.sh varcoef bin/x86_64-gnu/varcoef srun -n 4 ./my_varcoef_mpi`.

The output file has a 40-byte header (magic, variant name, nx, ny, nz, iterations) followed by the `nx·ny·nz` interior values as `float32`, `k` fastest, then `j`, then `i`. Keep this format in all your versions; for GPU versions, copy the solution back before writing it.

## Visualising results

`tools/plot_result.py` renders a result file as three orthogonal slices through the grid centre (`i`, `j` and `k` constant), using the same `(i,j,k)` indices as `compare`. It reads only those three planes, so it also works on multi-GiB results. It needs numpy and matplotlib ≥ 3.2 (on Deucalion: `ml SciPy-bundle matplotlib`).

```bash
tools/plot_result.py result.bin                          # solution -> result.png
tools/plot_result.py result.bin --input varcoef.dat      # varcoef: draw material boundaries
tools/plot_result.py mine.bin --diff ref.bin             # error map against the reference
tools/plot_result.py result.bin --at 64 10 64            # slices through another point
```

The `--diff` mode is the most useful one while optimising. It shows `|u − u_ref| / max|u_ref|` on a log scale with the `compare` tolerance marked. A correct version shows only faint rounding noise (~1e-7), while a bug shows up as a bright structure where it happens: a missing plane, a wrong boundary, a bad halo between threads, ranks or GPUs. When `compare` fails, this image usually tells you where to look.

The images show the field after the given number of iterations, not the converged solution. Jacobi needs O(N²) iterations to converge, so after the 100 iterations used for validation the field is still far from its steady state. This is expected and does not affect validation, which compares iterate against iterate.

## Rules

You may change anything in the implementation: data structures, memory layout, loop order, precision of intermediate computations, the way the boundary is handled, the way the convergence measure is computed, parallelisation, blocking in space and in time, the way the input is read, compiler and flags. You must keep the mathematical problem and the method:

- the same equation, discretisation (stencil and coefficients) and source term, with the problem data (initial guess, boundary values, source, materials) taken from the input file;
- the **Jacobi** iteration: each new iterate is computed only from the previous one. Gauss–Seidel, red-black ordering, SOR, multigrid or other solvers change the sequence of iterates and are not allowed. Fusing several Jacobi iterations into one pass over the data (temporal blocking) is allowed, because it computes exactly the same iterates;
- the command line, the input and output file formats, and the final `||u_new − u||₂` of the last iteration. The `-t` option must keep working. You may evaluate it less often than every iteration if you state how often in the report;
- storage of the solution in single precision (`float`), so that the single-precision rooflines of Task 1 apply.

Report timings of the iteration loop (the `solve time`). Report the load time separately; it is not part of the performance analysis unless you study I/O. For GPU versions, also report the host↔device transfer times separately and an end-to-end time that includes them.

## Job scripts

`scripts/job_arm.sh`, `scripts/job_x86.sh` and `scripts/job_gpu.sh` are Slurm templates for the `normal-arm`, `normal-x86` and `normal-a100-40` partitions. Replace the account (accounts end in `a`, `x` or `g` for the ARM, x86 and GPU partitions; run `billing` to see yours), the module lines and the dataset. Use the `dev-*` partitions only for compilation, dataset generation and short tests. CPU nodes are allocated exclusively. GPU nodes are shared unless you request `--exclusive`, which you need for scalability and energy measurements.
