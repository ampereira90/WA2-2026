#!/usr/bin/env python3
"""
ACA - Work Assignment 2 - input dataset generator.

Creates the input file read by the poisson7, box27 and varcoef programs.
The grid size is chosen from a required dataset size in MiB, or given
explicitly.

examples
  tools/gen_input.py varcoef --size-mib 1024                  # ~1 GiB input file
  tools/gen_input.py box27   --size-mib 1024 --measure memory --round up
                                                              # >= 1 GiB in the reference's memory
  tools/gen_input.py poisson7 --size-mib 0.03 --round down    # fits a 32 KiB L1 (reference layout)
  tools/gen_input.py poisson7 --dims 40 33 57 -o test.dat     # explicit (non-cubic) grid
  tools/gen_input.py box27 --size-mib 64 --dry-run            # only print the chosen size

--measure file    (default) the size of the input file
--measure memory  the memory used by the data structures of the REFERENCE
                  code. If you change the data structures, compute the
                  footprint of your own version.

Cubic grids are chosen (N^3 interior points); --round selects the grid
whose size is nearest to the target (default), the largest one not above
it (down, e.g. "fits in L1"), or the smallest one above it (up, e.g.
"> 1 GiB").

File format (little-endian; see common/common.h):
  header (56 B) | field table (32 B per field) | material table (56 B each)
  | fields, each starting at a 4096-byte aligned offset.
Each field covers the grid INCLUDING the halo, (nx+2)*(ny+2)*(nz+2) values,
k fastest, then j, then i. Halo values are the Dirichlet boundary values.

Requires numpy (on Deucalion: `ml SciPy-bundle`, or `pip install --user numpy`).
"""
import argparse
import os
import struct
import sys

try:
    import numpy as np
except ImportError:
    sys.exit("gen_input.py needs numpy (on Deucalion: `ml SciPy-bundle`, "
             "or `python3 -m pip install --user numpy`)")

MIB = 1024 * 1024
ALIGN = 4096
DTYPE_FLOAT32, DTYPE_UINT8 = 0, 1

# ----------------------------------------------------------------------------
# Problem definitions
# ----------------------------------------------------------------------------
# fields stored in the file: (name, dtype code, bytes per value)
FIELDS = {
    "poisson7": [("u", DTYPE_FLOAT32, 4)],
    "box27":    [("u", DTYPE_FLOAT32, 4), ("f", DTYPE_FLOAT32, 4)],
    "varcoef":  [("u", DTYPE_FLOAT32, 4), ("material", DTYPE_UINT8, 1)],
}
# bytes per grid point used by the data structures of the reference codes
REFERENCE_MEMORY = {"poisson7": 8, "box27": 16, "varcoef": 12}

# box27: Gaussian charges (x, y, z, width, strength)
CHARGES = [
    (0.30, 0.30, 0.30, 0.08,  400.0),
    (0.70, 0.60, 0.40, 0.12, -250.0),
    (0.45, 0.75, 0.70, 0.05,  900.0),
]

# varcoef: materials (name, conductivity, density, heat capacity, heat source)
MATERIALS = [
    ("epoxy",   1.0,  1200.0, 1100.0,   0.0),
    ("carbon", 20.0,  1800.0,  700.0,   0.0),
    ("air",     0.05,    1.2, 1000.0,   0.0),
    ("heater",  5.0,  7800.0,  500.0, 500.0),
]
MAT_MATRIX, MAT_FIBRE, MAT_PORE, MAT_HEATER = 0, 1, 2, 3


def plane_coords(i, ny, nz, h):
    """physical coordinates of plane i (x scalar, y column, z row)"""
    x = i * h
    y = (np.arange(ny + 2, dtype=np.float64) * h)[:, None]
    z = (np.arange(nz + 2, dtype=np.float64) * h)[None, :]
    return x, y, z


def interior_mask(i, nx, ny, nz):
    m = np.zeros((ny + 2, nz + 2), dtype=bool)
    if 1 <= i <= nx:
        m[1:ny + 1, 1:nz + 1] = True
    return m


def plane_poisson7(i, nx, ny, nz, h):
    # initial guess 0, homogeneous Dirichlet boundary
    return {"u": np.zeros((ny + 2, nz + 2), dtype=np.float32)}


def plane_box27(i, nx, ny, nz, h):
    x, y, z = plane_coords(i, ny, nz, h)
    f = np.zeros((ny + 2, nz + 2), dtype=np.float64)
    for cx, cy, cz, s, a in CHARGES:
        dx, dy, dz = x - cx, y - cy, z - cz
        f = f + a * np.exp(-(dx * dx + dy * dy + dz * dz) / (2.0 * s * s))
    f = f.astype(np.float32)
    f[~interior_mask(i, nx, ny, nz)] = 0.0
    return {"u": np.zeros((ny + 2, nz + 2), dtype=np.float32), "f": f}


def plane_varcoef(i, nx, ny, nz, h):
    x, y, z = plane_coords(i, ny, nz, h)
    shape = (ny + 2, nz + 2)
    zz = np.broadcast_to(z, shape)
    # same tests and priority as a chain of ifs: heater > pore > fibre > matrix
    mat = np.full(shape, MAT_MATRIX, dtype=np.uint8)
    mat[np.fmod(zz * 8.0, 1.0) < 0.3] = MAT_FIBRE
    px, py, pz = np.fmod(x, 0.2) - 0.1, np.fmod(y, 0.2) - 0.1, np.fmod(z, 0.2) - 0.1
    mat[(px * px + py * py + pz * pz) < 0.035 * 0.035] = MAT_PORE
    dx, dy, dz = x - 0.5, y - 0.5, z - 0.5
    mat[(dx * dx + dy * dy + dz * dz) < 0.12 * 0.12] = MAT_HEATER
    # boundary condition: u = 1 on the face i = 0, 0 elsewhere; initial guess 0
    u = np.full(shape, 1.0 if i == 0 else 0.0, dtype=np.float32)
    return {"u": u, "material": mat}


PLANE = {"poisson7": plane_poisson7, "box27": plane_box27, "varcoef": plane_varcoef}

# ----------------------------------------------------------------------------
# Sizes
# ----------------------------------------------------------------------------
def align(n):
    return (n + ALIGN - 1) // ALIGN * ALIGN


def layout(variant, nx, ny, nz):
    """byte offsets of each field and total file size"""
    nmat = len(MATERIALS) if variant == "varcoef" else 0
    offset = align(56 + 32 * len(FIELDS[variant]) + 56 * nmat)
    npts = (nx + 2) * (ny + 2) * (nz + 2)
    offsets = []
    for name, code, size in FIELDS[variant]:
        offsets.append(offset)
        offset = align(offset + npts * size)
    return offsets, offset


def file_bytes(variant, n):
    return layout(variant, n, n, n)[1]


def memory_bytes_dims(variant, nx, ny, nz):
    b = REFERENCE_MEMORY[variant] * (nx + 2) * (ny + 2) * (nz + 2)
    if variant == "poisson7":
        # two jagged float*** arrays: one pointer + one malloc'd row
        # (~16 B allocator overhead) per (i,j), one pointer array per i
        b += 2 * ((nx + 2) * (ny + 2) * (8 + 16) + (nx + 2) * (8 + 16))
    return b


def memory_bytes(variant, n):
    return memory_bytes_dims(variant, n, n, n)


def choose_n(variant, target, measure, rounding):
    size = file_bytes if measure == "file" else memory_bytes
    n = 1
    while size(variant, n + 1) <= target:
        n += 1
    # now size(n) <= target < size(n+1), unless even n = 1 is too big
    if size(variant, n) > target:
        return 1
    if rounding == "down":
        return n
    if rounding == "up":
        return n + 1
    return n if target - size(variant, n) <= size(variant, n + 1) - target else n + 1


# ----------------------------------------------------------------------------
# Writer
# ----------------------------------------------------------------------------
def write_dataset(path, variant, nx, ny, nz):
    h = 1.0 / (nx + 1)
    fields = FIELDS[variant]
    offsets, total = layout(variant, nx, ny, nz)
    nmat = len(MATERIALS) if variant == "varcoef" else 0
    plane_bytes = {name: (ny + 2) * (nz + 2) * size for name, code, size in fields}

    with open(path, "wb") as fp:
        fp.write(struct.pack("<8s16s4id2i", b"ACAWA2IN", variant.encode(),
                             nx, ny, nz, len(fields), h, nmat, 0))
        for (name, code, size), off in zip(fields, offsets):
            fp.write(struct.pack("<16s2iq", name.encode(), code, 0, off))
        if nmat:
            for name, kappa, rho, cp, q in MATERIALS:
                fp.write(struct.pack("<24s4d", name.encode(), kappa, rho, cp, q))
        fp.truncate(total)
        for i in range(nx + 2):
            data = PLANE[variant](i, nx, ny, nz, h)
            for (name, code, size), off in zip(fields, offsets):
                arr = np.ascontiguousarray(data[name])
                assert arr.nbytes == plane_bytes[name]
                fp.seek(off + i * plane_bytes[name])
                fp.write(arr.tobytes())
    return total


def main():
    ap = argparse.ArgumentParser(
        description="Create an input dataset for the WA2 3D Jacobi codes.",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    ap.add_argument("variant", choices=sorted(FIELDS))
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--size-mib", type=float, help="required dataset size in MiB")
    g.add_argument("--dims", type=int, nargs=3, metavar=("NX", "NY", "NZ"),
                   help="explicit number of interior points")
    ap.add_argument("--measure", choices=["file", "memory"], default="file",
                    help="what --size-mib refers to (default: file)")
    ap.add_argument("--round", dest="rounding", choices=["nearest", "down", "up"],
                    default="nearest", help="how to round to a cubic grid (default: nearest)")
    ap.add_argument("-o", "--output", help="output file (default: <variant>_<nx>x<ny>x<nz>.dat)")
    ap.add_argument("--dry-run", action="store_true", help="only print the chosen grid")
    args = ap.parse_args()

    if args.dims:
        nx, ny, nz = args.dims
        if min(nx, ny, nz) < 1:
            ap.error("dimensions must be >= 1")
    else:
        if args.size_mib <= 0:
            ap.error("--size-mib must be positive")
        n = choose_n(args.variant, args.size_mib * MIB, args.measure, args.rounding)
        nx = ny = nz = n

    total = layout(args.variant, nx, ny, nz)[1]
    mem = memory_bytes_dims(args.variant, nx, ny, nz)
    path = args.output or "{}_{}x{}x{}.dat".format(args.variant, nx, ny, nz)
    print("variant              : {}".format(args.variant))
    print("interior grid        : {} x {} x {}  ({:.3e} points, {:.3e} with halo)".format(
        nx, ny, nz, float(nx * ny * nz), float((nx + 2) * (ny + 2) * (nz + 2))))
    if args.size_mib:
        print("requested            : {:.4f} MiB ({}, round {})".format(
            args.size_mib, args.measure, args.rounding))
    print("input file size      : {:.4f} MiB".format(total / MIB))
    print("reference memory use : {:.4f} MiB (approx. {} B/point; recompute for your data structures)".format(
        mem / MIB, REFERENCE_MEMORY[args.variant]))
    if args.dry_run:
        return
    write_dataset(path, args.variant, nx, ny, nz)
    print("written              : {}".format(path))


if __name__ == "__main__":
    main()
