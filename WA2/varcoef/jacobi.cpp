/*
 * ACA - Work Assignment 2 - variant "varcoef": Jacobi kernels
 */
#include "varcoef.h"

#include <cmath>

static double harmonic_mean(double a, double b)
{
    return 2.0 * a * b / (a + b);
}

void jacobi_sweep(Grid *g)
{
    const double h2 = g->h * g->h;
    for (int k = 1; k <= g->nz; k++)
        for (int j = 1; j <= g->ny; j++)
            for (int i = 1; i <= g->nx; i++) {
                long p = idx(g, i, j, k);
                long nb[6] = { idx(g, i - 1, j, k), idx(g, i + 1, j, k),
                               idx(g, i, j - 1, k), idx(g, i, j + 1, k),
                               idx(g, i, j, k - 1), idx(g, i, j, k + 1) };
                double kp = g->materials[g->material[p]].conductivity;
                double num = 0.0, den = 0.0;
                for (int n = 0; n < 6; n++) {
                    double kq  = g->materials[g->material[nb[n]]].conductivity;
                    double kpq = harmonic_mean(kp, kq);
                    num += kpq * g->u[nb[n]];
                    den += kpq;
                }
                num += h2 * g->materials[g->material[p]].heat_source;
                g->u_new[p] = num / den;
            }
}

double jacobi_residual(const Grid *g)
{
    double sum = 0.0;
    for (int k = 1; k <= g->nz; k++)
        for (int j = 1; j <= g->ny; j++)
            for (int i = 1; i <= g->nx; i++) {
                long p = idx(g, i, j, k);
                sum += std::pow(g->u_new[p] - g->u[p], 2.0);
            }
    return std::sqrt(sum);
}

void jacobi_update(Grid *g)
{
    for (int i = 0; i <= g->nx + 1; i++)
        for (int j = 0; j <= g->ny + 1; j++)
            for (int k = 0; k <= g->nz + 1; k++) {
                long p = idx(g, i, j, k);
                g->u[p] = g->u_new[p];
            }
}
