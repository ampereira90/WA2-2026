/*
 * ACA - Work Assignment 2 - variant "poisson7": grid management
 */
#include "poisson7.h"

#include <cstdlib>
#include <cstdio>
#include <vector>

static float ***alloc_3d(int n1, int n2, int n3)
{
    float ***a = (float ***) std::malloc(n1 * sizeof(float **));
    if (!a) { std::fprintf(stderr, "out of memory\n"); std::exit(1); }
    for (int i = 0; i < n1; i++) {
        a[i] = (float **) std::malloc(n2 * sizeof(float *));
        if (!a[i]) { std::fprintf(stderr, "out of memory\n"); std::exit(1); }
        for (int j = 0; j < n2; j++) {
            a[i][j] = (float *) std::malloc(n3 * sizeof(float));
            if (!a[i][j]) { std::fprintf(stderr, "out of memory\n"); std::exit(1); }
        }
    }
    return a;
}

static void free_3d(float ***a, int n1, int n2)
{
    for (int i = 0; i < n1; i++) {
        for (int j = 0; j < n2; j++) std::free(a[i][j]);
        std::free(a[i]);
    }
    std::free(a);
}

void grid_alloc(Grid *g, int nx, int ny, int nz, double h)
{
    g->nx = nx; g->ny = ny; g->nz = nz;
    g->h  = h;
    g->u     = alloc_3d(nx + 2, ny + 2, nz + 2);
    g->u_new = alloc_3d(nx + 2, ny + 2, nz + 2);
}

void grid_load(Grid *g, InputFile *in)
{
    /* initial guess and boundary values, one i-plane at a time */
    std::vector<float> plane(in->plane_size());
    for (int i = 0; i <= g->nx + 1; i++) {
        in->read_planes("u", i, 1, &plane[0]);
        for (int j = 0; j <= g->ny + 1; j++)
            for (int k = 0; k <= g->nz + 1; k++) {
                float v = plane[(long) j * (g->nz + 2) + k];
                g->u[i][j][k]     = v;
                g->u_new[i][j][k] = v;
            }
    }
}

void grid_free(Grid *g)
{
    free_3d(g->u,     g->nx + 2, g->ny + 2);
    free_3d(g->u_new, g->nx + 2, g->ny + 2);
}

void grid_emit(const Grid *g, ResultSink *sink)
{
    for (int i = 1; i <= g->nx; i++)
        for (int j = 1; j <= g->ny; j++)
            sink->row(&g->u[i][j][1], g->nz);
}
