/*
 * ACA - Work Assignment 2 - variant "poisson7"
 *
 * Jacobi iteration for the 3D Poisson equation  -lap(u) = f  on a regular
 * grid with spacing h, homogeneous Dirichlet boundary conditions (u = 0 on
 * the boundary) and the second-order 7-point discretisation of the Laplacian:
 *
 *   u_new(i,j,k) = ( u(i-1,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i,j+1,k)
 *                  + u(i,j,k-1) + u(i,j,k+1) + h^2 f(x_i,y_j,z_k) ) / 6
 *
 * with f(x,y,z) = 3 pi^2 sin(pi x) sin(pi y) sin(pi z),
 * x_i = i h, y_j = j h, z_k = k h.
 *
 * Input dataset (tools/gen_input.py poisson7): grid size, h and one field
 *   "u"  float32  initial guess, with the boundary values in the halo
 * The source term f is evaluated by the program.
 *
 * The grid stores nx*ny*nz interior points plus a one-point halo that holds
 * the boundary values, i.e. indices 0..nx+1, 0..ny+1, 0..nz+1.
 */
#ifndef POISSON7_H
#define POISSON7_H

#include "../common/common.h"

#define VARIANT_NAME "poisson7"
#define VARIANT_DESC "3D Poisson, 7-point stencil, Jacobi iteration"

struct Grid {
    int     nx, ny, nz;   /* interior points                      */
    double  h;            /* grid spacing                         */
    float ***u;           /* current iterate,  u[i][j][k]         */
    float ***u_new;       /* next iterate                         */
};

/* grid.cpp */
void grid_alloc(Grid *g, int nx, int ny, int nz, double h);
void grid_load(Grid *g, InputFile *in);
void grid_free(Grid *g);
void grid_emit(const Grid *g, ResultSink *sink);

/* jacobi.cpp */
double source(double x, double y, double z);
void   jacobi_sweep(Grid *g);            /* u_new <- J(u)                 */
double jacobi_residual(const Grid *g);   /* ||u_new - u||_2               */
void   jacobi_update(Grid *g);           /* u <- u_new                    */

#endif
