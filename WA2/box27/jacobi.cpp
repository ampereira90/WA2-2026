/*
 * ACA - Work Assignment 2 - variant "box27": Jacobi kernels
 */
#include "box27.h"

#include <cmath>

void jacobi_sweep(Grid *g)
{
    for (int i = 0; i <= g->nx + 1; i++)
        for (int j = 0; j <= g->ny + 1; j++)
            for (int k = 0; k <= g->nz + 1; k++) {
                Cell *c = cell_at(g, i, j, k);
                if (c->type == BOUNDARY) {
                    c->u_new = c->u;
                    continue;
                }
                double sum = 0.0;
                for (int di = -1; di <= 1; di++)
                    for (int dj = -1; dj <= 1; dj++)
                        for (int dk = -1; dk <= 1; dk++) {
                            if (di == 0 && dj == 0 && dk == 0) continue;
                            sum += stencil_weight(di, dj, dk)
                                 * cell_at(g, i + di, j + dj, k + dk)->u;
                        }
                c->u_new = (sum + g->h * g->h * c->f) / stencil_weight(0, 0, 0);
            }
}

double jacobi_residual(const Grid *g)
{
    double sum = 0.0;
    for (int i = 0; i <= g->nx + 1; i++)
        for (int j = 0; j <= g->ny + 1; j++)
            for (int k = 0; k <= g->nz + 1; k++) {
                const Cell *c = cell_at(g, i, j, k);
                if (c->type == INTERIOR)
                    sum += std::pow(c->u_new - c->u, 2.0);
            }
    return std::sqrt(sum);
}

void jacobi_update(Grid *g)
{
    long n = (long)(g->nx + 2) * (g->ny + 2) * (g->nz + 2);
    for (long p = 0; p < n; p++)
        g->cells[p].u = g->cells[p].u_new;
}
