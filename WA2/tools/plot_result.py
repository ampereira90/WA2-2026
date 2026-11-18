#!/usr/bin/env python3
"""
ACA - Work Assignment 2 - render a result file as an image.

Draws three orthogonal slices of the solution u through a point of the grid
(default: the centre), with one shared colour scale. Only the three planes
are read from disk, so multi-GiB result files are fine.

examples
  tools/plot_result.py result.bin                         # -> result.png
  tools/plot_result.py result.bin --input varcoef.dat     # + material boundaries
  tools/plot_result.py mine.bin --diff ref.bin            # error map vs. the reference
  tools/plot_result.py result.bin --at 64 10 64 -o s.png  # slices through (i,j,k)=(64,10,64)

The slices are:  i = const (a j-k plane), j = const (an i-k plane) and
k = const (an i-j plane). Coordinates are the interior indices 1..n used by
tools/compare, with x = i*h, y = j*h, z = k*h.

Signed solutions use a diverging colour map centred at 0; non-negative ones a
sequential map from 0. With --diff, |u - u_ref| / max|u_ref| is shown on a
logarithmic scale (1e-8 to 1e-2) with the absolute tolerance of tools/compare
marked on the colour bar: rounding differences appear as faint ~1e-7 noise,
bugs as bright structures (a missing plane, a wrong boundary, a shifted
stencil, a bad halo between threads or ranks).

Requires numpy and matplotlib (on Deucalion: `ml SciPy-bundle matplotlib`).
"""
import argparse
import struct
import sys

try:
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, TwoSlopeNorm
except ImportError:
    sys.exit("plot_result.py needs numpy and matplotlib "
             "(on Deucalion: `ml SciPy-bundle matplotlib`)")

RESULT_HEADER = 40


def open_result(path):
    with open(path, "rb") as fp:
        head = fp.read(RESULT_HEADER)
    if len(head) < RESULT_HEADER or head[:6] != b"ACAWA2" or head[6:8] == b"IN":
        sys.exit("{}: not a WA2 result file".format(path))
    variant = head[8:24].split(b"\0")[0].decode()
    nx, ny, nz, iters = struct.unpack("<4i", head[24:40])
    u = np.memmap(path, dtype="<f4", mode="r", offset=RESULT_HEADER, shape=(nx, ny, nz))
    return variant, (nx, ny, nz), iters, u


def open_material(path, dims):
    """material ids of the interior from a varcoef input file, or None"""
    with open(path, "rb") as fp:
        head = fp.read(56)
        if head[:8] != b"ACAWA2IN":
            sys.exit("{}: not a WA2 input file".format(path))
        variant = head[8:24].split(b"\0")[0].decode()
        nx, ny, nz, nfields = struct.unpack("<4i", head[24:40])
        if (nx, ny, nz) != tuple(dims):
            sys.exit("{}: grid {}x{}x{} does not match the result".format(path, nx, ny, nz))
        fields = {}
        for _ in range(nfields):
            name, dtype, _r, off = struct.unpack("<16s2iq", fp.read(32))
            fields[name.split(b"\0")[0].decode()] = (dtype, off)
    if "material" not in fields:
        print("note: {} ({}) has no material field; --input ignored".format(path, variant))
        return None
    full = np.memmap(path, dtype=np.uint8, mode="r", offset=fields["material"][1],
                     shape=(nx + 2, ny + 2, nz + 2))
    return full[1:-1, 1:-1, 1:-1]


def slices(a, at):
    """the three planes through index 'at' (0-based), each as a 2-D array"""
    i, j, k = at
    return [np.asarray(a[i, :, :]),      # j-k plane
            np.asarray(a[:, j, :]),      # i-k plane
            np.asarray(a[:, :, k])]      # i-j plane


def main():
    ap = argparse.ArgumentParser(description="Render a WA2 result file as an image.",
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 epilog=__doc__)
    ap.add_argument("result")
    ap.add_argument("--input", help="dataset used for the run (draws material boundaries, varcoef)")
    ap.add_argument("--diff", metavar="REF", help="show the error relative to this reference result")
    ap.add_argument("--at", type=int, nargs=3, metavar=("I", "J", "K"),
                    help="slice through this interior point (1-based; default: centre)")
    ap.add_argument("-o", "--output", help="image file (default: <result>.png)")
    ap.add_argument("--dpi", type=int, default=200)
    args = ap.parse_args()

    variant, dims, iters, u = open_result(args.result)
    nx, ny, nz = dims
    at = [d // 2 for d in dims] if args.at is None else [a - 1 for a in args.at]
    if any(not 0 <= a < d for a, d in zip(at, dims)):
        ap.error("--at must be within 1..nx, 1..ny, 1..nz")

    planes = slices(u, at)
    if args.diff:
        rvariant, rdims, riters, ref = open_result(args.diff)
        if (rvariant, rdims, riters) != (variant, dims, iters):
            sys.exit("reference and result differ in variant, grid or iterations")
        scale = float(np.max(np.abs(ref))) or 1.0       # full pass over the reference
        planes = [np.abs(p - q) / scale for p, q in zip(planes, slices(ref, at))]
        floor = 1e-8
        planes = [np.maximum(p, floor) for p in planes]
        norm = LogNorm(vmin=floor, vmax=max(1e-2, max(float(p.max()) for p in planes)))
        cmap, clabel = "magma", "|u − u_ref| / max|u_ref|"
    else:
        lo = min(float(p.min()) for p in planes)
        hi = max(float(p.max()) for p in planes)
        if lo < 0 < hi:
            m = max(-lo, hi)
            norm, cmap = TwoSlopeNorm(vcenter=0.0, vmin=-m, vmax=m), "RdBu_r"
        else:
            norm, cmap = plt.Normalize(vmin=min(lo, 0.0), vmax=hi if hi > 0 else 1.0), "inferno"
        clabel = "u"

    mats = None
    if args.input and not args.diff:
        m3 = open_material(args.input, dims)
        if m3 is not None:
            mats = slices(m3, at)

    # (plane, horizontal axis, vertical axis, sizes, fixed index text)
    spec = [(planes[0], "k (z)", "j (y)", (nz, ny), "i = {}".format(at[0] + 1)),
            (planes[1], "k (z)", "i (x)", (nz, nx), "j = {}".format(at[1] + 1)),
            (planes[2], "j (y)", "i (x)", (ny, nx), "k = {}".format(at[2] + 1))]

    plt.rcParams.update({"font.size": 9, "axes.titlesize": 9, "axes.labelsize": 9,
                         "xtick.labelsize": 7, "ytick.labelsize": 7,
                         "xtick.direction": "out", "ytick.direction": "out"})
    widths = [s[3][0] for s in spec]
    # figure height follows the tallest panel's aspect ratio
    panel_w = 8.4 * max(widths) / sum(widths)
    aspect = max(s[3][1] / float(s[3][0]) for s in spec)
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 1.3 + panel_w * aspect),
                             gridspec_kw={"width_ratios": widths, "wspace": 0.35})
    for n, (ax, (img, hx, vy, (w, hgt), fixed)) in enumerate(zip(axes, spec)):
        im = ax.imshow(img, origin="lower", cmap=cmap, norm=norm, interpolation="nearest",
                       extent=(0.5, w + 0.5, 0.5, hgt + 0.5), aspect="equal")
        if mats is not None:
            ax.contour(np.arange(1, w + 1), np.arange(1, hgt + 1), mats[n].astype(float),
                       levels=[0.5, 1.5, 2.5], colors="white", linewidths=0.4, alpha=0.7)
        ax.set_xlabel(hx)
        ax.set_ylabel(vy)
        ax.set_title(fixed, loc="left")
    cb = fig.colorbar(im, ax=list(axes), shrink=0.85, pad=0.02)
    cb.set_label(clabel)
    if args.diff:
        # tools/compare fails a point above 1e-5*max|ref| (+ 1e-3*|ref|)
        cb.ax.axhline(1e-5, color="white", linewidth=1.2)
        cb.set_label(clabel + "\n(white line: tools/compare absolute tolerance)", fontsize=8)
    what = "error vs. reference" if args.diff else "solution"
    fig.canvas.draw()
    top = axes[0].get_position().y1
    fig.text(axes[0].get_position().x0, top + 0.6 / fig.get_figheight(),
             "{} {}, {}×{}×{} grid, {} iterations".format(variant, what, nx, ny, nz, iters),
             ha="left", va="bottom", fontsize=9)
    out = args.output or (args.result.rsplit(".", 1)[0] + (".diff.png" if args.diff else ".png"))
    fig.savefig(out, dpi=args.dpi, bbox_inches="tight")
    print("written {}".format(out))


if __name__ == "__main__":
    main()
