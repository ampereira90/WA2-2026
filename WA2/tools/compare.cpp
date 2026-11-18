/*
 * ACA - Work Assignment 2 - result validation
 *
 *   compare REFERENCE.bin CANDIDATE.bin [rtol [atol]]
 *
 * Both files must come from the same variant, grid size and iteration
 * count. A candidate PASSES when
 *   (1) every point satisfies  |c - r| <= rtol * |r| + atol * max|r|
 *       (defaults rtol = 1e-3, atol = 1e-5), and
 *   (2) the relative L2 error  ||c - r||_2 / ||r||_2  is <= 1e-4, and
 *   (3) there are no NaN / Inf values.
 * The point-wise test (1) catches errors confined to small regions (e.g. a
 * missing plane next to the boundary, where the solution is small), which a
 * global norm alone would hide. The worst point is reported to help locate
 * halo and loop-bound bugs. Exit status 0 = PASS, 1 = FAIL, 2 = usage error.
 *
 * Floating-point results are not expected to be bit-identical: changing the
 * order of the operations (vectorisation, FMA contraction, single instead
 * of double precision intermediates, reductions on GPUs, ...) changes the
 * rounding, and the differences grow slowly with the number of iterations.
 * The tolerances allow for that up to a few thousand iterations; validate
 * with the iteration counts given in README.md (-i 100).
 */
#include "../common/common.h"

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

static FILE *open_result(const char *path, ResultHeader *h)
{
    FILE *fp = std::fopen(path, "rb");
    if (!fp) { std::perror(path); std::exit(2); }
    if (std::fread(h, sizeof *h, 1, fp) != 1 || std::memcmp(h->magic, "ACAWA2", 6) != 0) {
        std::fprintf(stderr, "%s: not a WA2 result file\n", path);
        std::exit(2);
    }
    h->variant[sizeof h->variant - 1] = '\0';
    return fp;
}

/* print the (i,j,k) interior coordinates (1-based) of linear index p */
static void print_ijk(long p, const ResultHeader &h)
{
    if (p < 0) return;
    long k = p % h.nz, j = (p / h.nz) % h.ny, i = p / ((long) h.nz * h.ny);
    std::printf("  at (i,j,k) = (%ld,%ld,%ld)", i + 1, j + 1, k + 1);
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        std::fprintf(stderr, "usage: %s REFERENCE.bin CANDIDATE.bin [rtol [atol]]\n", argv[0]);
        return 2;
    }
    const double rtol   = argc > 3 ? std::atof(argv[3]) : 1e-3;
    const double atol   = argc > 4 ? std::atof(argv[4]) : 1e-5;
    const double tol_l2 = 1e-4;

    ResultHeader ha, hb;
    FILE *fa = open_result(argv[1], &ha);
    FILE *fb = open_result(argv[2], &hb);

    std::printf("reference : %s  %s %dx%dx%d, %d iterations\n",
                argv[1], ha.variant, ha.nx, ha.ny, ha.nz, ha.iters);
    std::printf("candidate : %s  %s %dx%dx%d, %d iterations\n",
                argv[2], hb.variant, hb.nx, hb.ny, hb.nz, hb.iters);
    if (std::strcmp(ha.variant, hb.variant) || ha.nx != hb.nx || ha.ny != hb.ny ||
        ha.nz != hb.nz || ha.iters != hb.iters) {
        std::printf("FAIL: headers differ (variant, grid size or iteration count)\n");
        return 1;
    }

    const size_t chunk = 1 << 20;
    const long   n = (long) ha.nx * ha.ny * ha.nz;
    std::vector<float> a(chunk), b(chunk);

    /* pass 1: max |ref| */
    double max_ref = 0.0;
    for (long done = 0; done < n; ) {
        size_t m = (size_t)((n - done) < (long) chunk ? (n - done) : (long) chunk);
        if (std::fread(&a[0], sizeof(float), m, fa) != m) { std::fprintf(stderr, "unexpected end of file\n"); return 2; }
        for (size_t q = 0; q < m; q++) if (std::fabs((double) a[q]) > max_ref) max_ref = std::fabs((double) a[q]);
        done += (long) m;
    }
    std::fseek(fa, (long) sizeof(ResultHeader), SEEK_SET);
    const double abs_floor = atol * max_ref;

    /* pass 2: errors */
    long   nonfinite = 0, violations = 0, worst = -1, worst_abs_at = -1;
    double max_abs = 0.0, worst_ratio = 0.0, sum_d2 = 0.0, sum_r2 = 0.0;
    for (long done = 0; done < n; ) {
        size_t m = (size_t)((n - done) < (long) chunk ? (n - done) : (long) chunk);
        if (std::fread(&a[0], sizeof(float), m, fa) != m ||
            std::fread(&b[0], sizeof(float), m, fb) != m) {
            std::fprintf(stderr, "unexpected end of file\n");
            return 2;
        }
        for (size_t q = 0; q < m; q++) {
            double r = a[q], c = b[q];
            if (!std::isfinite(c)) { nonfinite++; continue; }
            double d = std::fabs(c - r);
            double allowed = rtol * std::fabs(r) + abs_floor;
            double ratio = allowed > 0 ? d / allowed : (d > 0 ? 1e300 : 0.0);
            if (ratio > 1.0) violations++;
            if (ratio > worst_ratio) { worst_ratio = ratio; worst = done + (long) q; }
            if (d > max_abs) { max_abs = d; worst_abs_at = done + (long) q; }
            sum_d2 += d * d;
            sum_r2 += r * r;
        }
        done += (long) m;
    }
    double max_rel = max_ref > 0 ? max_abs / max_ref : max_abs;
    double l2_rel  = sum_r2 > 0 ? std::sqrt(sum_d2 / sum_r2) : std::sqrt(sum_d2);

    std::printf("max |ref|            : %.6e\n", max_ref);
    std::printf("max abs error        : %.6e", max_abs);
    print_ijk(worst_abs_at, ha);
    std::printf("\nmax error / max|ref| : %.6e\n", max_rel);
    std::printf("relative L2 error    : %.6e  (tol %.1e)\n", l2_rel, tol_l2);
    std::printf("point-wise test      : %ld of %ld points outside rtol=%.1e, atol=%.1e*max|ref|\n",
                violations, n, rtol, atol);
    if (worst >= 0 && violations > 0) {
        std::printf("worst point          :");
        print_ijk(worst, ha);
        std::fseek(fa, (long) sizeof(ResultHeader) + worst * (long) sizeof(float), SEEK_SET);
        std::fseek(fb, (long) sizeof(ResultHeader) + worst * (long) sizeof(float), SEEK_SET);
        float ra = 0, cb = 0;
        if (std::fread(&ra, sizeof ra, 1, fa) == 1 && std::fread(&cb, sizeof cb, 1, fb) == 1)
            std::printf("  ref = %.7e  candidate = %.7e", ra, cb);
        std::printf("\n");
    }
    std::printf("non-finite values    : %ld\n", nonfinite);

    bool ok = nonfinite == 0 && violations == 0 && l2_rel <= tol_l2;
    std::printf("%s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
