/*
 * ACA - Work Assignment 2 - variant "varcoef": grid and material setup
 */
#include "varcoef.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

void grid_alloc(Grid *g, int nx, int ny, int nz, double h)
{
    g->nx = nx; g->ny = ny; g->nz = nz;
    g->h  = h;
    long n = (long)(nx + 2) * (ny + 2) * (nz + 2);
    g->u        = (float *) std::malloc(n * sizeof(float));
    g->u_new    = (float *) std::malloc(n * sizeof(float));
    g->material = (int *)   std::malloc(n * sizeof(int));
    g->materials = (Material *) std::malloc(MAX_MATERIALS * sizeof(Material));
    g->nmaterials = 0;
    if (!g->u || !g->u_new || !g->material || !g->materials) {
        std::fprintf(stderr, "out of memory\n");
        std::exit(1);
    }
}

void grid_load(Grid *g, InputFile *in)
{
    g->nmaterials = in->nmaterials();
    for (int m = 0; m < g->nmaterials; m++) {
        const MaterialRecord &r = in->material(m);
        Material *mat = &g->materials[m];
        std::memset(mat, 0, sizeof *mat);
        std::memcpy(mat->name, r.name, sizeof mat->name - 1);   /* r.name is NUL-terminated */
        mat->conductivity  = r.conductivity;
        mat->density       = r.density;
        mat->heat_capacity = r.heat_capacity;
        mat->heat_source   = r.heat_source;
        mat->id            = m;
    }
    std::vector<float>   u0(in->plane_size());
    std::vector<uint8_t> id(in->plane_size());
    for (int i = 0; i <= g->nx + 1; i++) {
        in->read_planes("u", i, 1, &u0[0]);
        in->read_planes("material", i, 1, &id[0]);
        for (int j = 0; j <= g->ny + 1; j++)
            for (int k = 0; k <= g->nz + 1; k++) {
                long q = (long) j * (g->nz + 2) + k;
                long p = idx(g, i, j, k);
                if (id[q] >= g->nmaterials) {
                    std::fprintf(stderr, "invalid material id %d\n", id[q]);
                    std::exit(1);
                }
                g->u[p]        = u0[q];
                g->u_new[p]    = u0[q];
                g->material[p] = id[q];
            }
    }
}

void grid_free(Grid *g)
{
    std::free(g->u);
    std::free(g->u_new);
    std::free(g->material);
    std::free(g->materials);
}

void grid_emit(const Grid *g, ResultSink *sink)
{
    for (int i = 1; i <= g->nx; i++)
        for (int j = 1; j <= g->ny; j++)
            sink->row(&g->u[idx(g, i, j, 1)], g->nz);
}
