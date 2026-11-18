/*
 * ACA - Work Assignment 2 - variant "box27"
 *
 * Jacobi iteration for the 3D Poisson equation  -lap(u) = f  using the
 * compact 27-point ("box") discretisation of the Laplacian:
 *
 *   lap(u)(p) ~ 1/(30 h^2) * ( -128 u(p) + 14 * sum_{6 face neighbours}  u
 *                                        +  3 * sum_{12 edge neighbours} u
 *                                        +  1 * sum_{8 corner neighbours} u )
 *
 * A neighbour q = p + (di,dj,dk), di,dj,dk in {-1,0,1}, is a face, edge or
 * corner neighbour when |di|+|dj|+|dk| = 1, 2 or 3. One Jacobi iteration is
 *
 *   u_new(p) = ( sum_{q != p} w(q) u(q) + h^2 f(p) ) / w(p)
 *
 * with w(p) = 128/30 and w(q) = 14/30, 3/30 and 1/30 for face, edge and
 * corner neighbours. Homogeneous Dirichlet boundary conditions (u = 0).
 *
 * Input dataset (tools/gen_input.py box27): grid size, h and two fields
 *   "u"  float32  initial guess, with the boundary values in the halo
 *   "f"  float32  source term (the generator places three Gaussian charges)
 * The grid holds (nx+2) x (ny+2) x (nz+2) cells, including a one-cell halo
 * that stores the boundary values.
 */
#ifndef BOX27_H
#define BOX27_H

#include "../common/common.h"

#define VARIANT_NAME "box27"
#define VARIANT_DESC "3D Poisson, 27-point compact stencil, Jacobi iteration"

enum CellType { INTERIOR = 0, BOUNDARY = 1 };

struct Cell {
    float u;        /* current iterate   */
    float u_new;    /* next iterate      */
    float f;        /* source term       */
    int   type;     /* INTERIOR/BOUNDARY */
};

struct Grid {
    int    nx, ny, nz;   /* interior cells */
    double h;            /* grid spacing   */
    Cell  *cells;        /* (nx+2)*(ny+2)*(nz+2) cells, k fastest */
};

static inline Cell *cell_at(const Grid *g, int i, int j, int k)
{
    long idx = ((long) i * (g->ny + 2) + j) * (g->nz + 2) + k;
    return &g->cells[idx];
}

/* stencil.cpp */
double stencil_weight(int di, int dj, int dk);

/* grid.cpp */
void grid_alloc(Grid *g, int nx, int ny, int nz, double h);
void grid_load(Grid *g, InputFile *in);
void grid_free(Grid *g);
void grid_emit(const Grid *g, ResultSink *sink);

/* jacobi.cpp */
void   jacobi_sweep(Grid *g);
double jacobi_residual(const Grid *g);
void   jacobi_update(Grid *g);

#endif
