/*
 * ACA - Work Assignment 2 - variant "varcoef": driver
 */
#include "varcoef.h"

#include <cstdio>

int main(int argc, char **argv)
{
    Options opt;
    parse_options(argc, argv, &opt, VARIANT_NAME, VARIANT_DESC);

    double t0 = wall_time();
    Grid g;
    {
        InputFile in(opt.input, VARIANT_NAME);
        opt.nx = in.nx(); opt.ny = in.ny(); opt.nz = in.nz();
        grid_alloc(&g, opt.nx, opt.ny, opt.nz, in.h());
        grid_load(&g, &in);
    }
    double t_load = wall_time() - t0;

    double residual = 0.0;
    int it;
    t0 = wall_time();
    for (it = 1; it <= opt.max_iters; it++) {
        jacobi_sweep(&g);
        residual = jacobi_residual(&g);
        jacobi_update(&g);
        if (opt.print_every > 0 && it % opt.print_every == 0)
            std::printf("iter %6d  ||du||_2 = %.6e\n", it, residual);
        if (residual < opt.tol) break;
    }
    if (it > opt.max_iters) it = opt.max_iters;
    double t_solve = wall_time() - t0;

    ResultSink sink(opt.output, VARIANT_NAME, opt.nx, opt.ny, opt.nz, it);
    grid_emit(&g, &sink);
    print_report(VARIANT_NAME, &opt, it, residual, t_load, t_solve, &sink);

    grid_free(&g);
    return 0;
}
