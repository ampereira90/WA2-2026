/*
 * ACA - Work Assignment 2 - variant "poisson7": Jacobi kernels
 */
#include "poisson7.h"

#include <cmath>

static const double PI = 3.14159265358979323846;

double source(double x, double y, double z)
{
    return 3.0 * PI * PI * std::sin(PI * x) * std::sin(PI * y) * std::sin(PI * z);
}

void jacobi_sweep(Grid *g)
{
    for (int i = 0; i <= g->nx + 1; i++)
        for (int k = 0; k <= g->nz + 1; k++)
            for (int j = 0; j <= g->ny + 1; j++) {
                if (i == 0 || i == g->nx + 1 ||
                    j == 0 || j == g->ny + 1 ||
                    k == 0 || k == g->nz + 1) {
                    /* boundary point: Dirichlet value is kept */
                    g->u_new[i][j][k] = g->u[i][j][k];
                    continue;
                }
                double x = i * g->h, y = j * g->h, z = k * g->h;
                g->u_new[i][j][k] =
                    ( g->u[i-1][j][k] + g->u[i+1][j][k]
                    + g->u[i][j-1][k] + g->u[i][j+1][k]
                    + g->u[i][j][k-1] + g->u[i][j][k+1]
                    + g->h * g->h * source(x, y, z) ) / 6.0;
            }
}

double jacobi_residual(const Grid *g)
{
    double sum = 0.0;
    for (int i = 1; i <= g->nx; i++)
        for (int k = 1; k <= g->nz; k++)
            for (int j = 1; j <= g->ny; j++)
                sum += std::pow(g->u_new[i][j][k] - g->u[i][j][k], 2.0);
    return std::sqrt(sum);
}

void jacobi_update(Grid *g)
{
    for (int i = 0; i <= g->nx + 1; i++)
        for (int k = 0; k <= g->nz + 1; k++)
            for (int j = 0; j <= g->ny + 1; j++)
                g->u[i][j][k] = g->u_new[i][j][k];
}
