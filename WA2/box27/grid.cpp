/*
 * ACA - Work Assignment 2 - variant "box27": grid management
 */
#include "box27.h"

#include <cstdio>
#include <cstdlib>
#include <vector>

void grid_alloc(Grid *g, int nx, int ny, int nz, double h)
{
    g->nx = nx; g->ny = ny; g->nz = nz;
    g->h  = h;
    long n = (long)(nx + 2) * (ny + 2) * (nz + 2);
    g->cells = (Cell *) std::malloc(n * sizeof(Cell));
    if (!g->cells) { std::fprintf(stderr, "out of memory\n"); std::exit(1); }
}

void grid_load(Grid *g, InputFile *in)
{
    std::vector<float> u0(in->plane_size()), f(in->plane_size());
    for (int i = 0; i <= g->nx + 1; i++) {
        in->read_planes("u", i, 1, &u0[0]);
        in->read_planes("f", i, 1, &f[0]);
        for (int j = 0; j <= g->ny + 1; j++)
            for (int k = 0; k <= g->nz + 1; k++) {
                long q  = (long) j * (g->nz + 2) + k;
                Cell *c = cell_at(g, i, j, k);
                bool boundary = (i == 0 || i == g->nx + 1 ||
                                 j == 0 || j == g->ny + 1 ||
                                 k == 0 || k == g->nz + 1);
                c->type  = boundary ? BOUNDARY : INTERIOR;
                c->u     = u0[q];
                c->u_new = u0[q];
                c->f     = boundary ? 0.0f : f[q];
            }
    }
}

void grid_free(Grid *g)
{
    std::free(g->cells);
}

void grid_emit(const Grid *g, ResultSink *sink)
{
    std::vector<float> row(g->nz);
    for (int i = 1; i <= g->nx; i++)
        for (int j = 1; j <= g->ny; j++) {
            for (int k = 1; k <= g->nz; k++) row[k - 1] = cell_at(g, i, j, k)->u;
            sink->row(&row[0], g->nz);
        }
}
