/*
 * ACA - Work Assignment 2 - variant "varcoef"
 *
 * Jacobi iteration for steady-state heat conduction in a heterogeneous
 * material,   -div( kappa(x) grad u ) = q(x),
 * discretised with a 7-point finite-volume stencil. The conductivity on the
 * face between cells p and q is the harmonic mean of the cell values:
 *
 *   kappa_pq = 2 kappa_p kappa_q / (kappa_p + kappa_q)
 *
 *   u_new(p) = ( sum_{q in 6 neighbours} kappa_pq u(q) + h^2 q(p) )
 *              / ( sum_{q in 6 neighbours} kappa_pq )
 *
 * Each cell belongs to one material; kappa and q are properties of the
 * material. Input dataset (tools/gen_input.py varcoef): grid size, h, the
 * material table and two fields
 *   "u"         float32  initial guess, with the Dirichlet boundary values
 *                        in the halo (the generator uses u = 1 on the face
 *                        i = 0 and u = 0 on the other faces)
 *   "material"  uint8    material id of every cell, halo included (the
 *                        generator builds a layered composite with pores
 *                        and a central heater)
 *
 * Arrays have (nx+2) x (ny+2) x (nz+2) entries, including a one-cell halo
 * that stores the boundary values; k is the fastest index.
 */
#ifndef VARCOEF_H
#define VARCOEF_H

#include "../common/common.h"

#define VARIANT_NAME "varcoef"
#define VARIANT_DESC "3D heterogeneous heat conduction, 7-point FV stencil, Jacobi iteration"

struct Material {
    char   name[24];
    double conductivity;     /* kappa  [W/(m K)]   */
    double density;          /* rho    [kg/m^3]    */
    double heat_capacity;    /* c_p    [J/(kg K)]  */
    double heat_source;      /* q      [W/m^3]     */
    int    id;
};

#define MAX_MATERIALS 256

struct Grid {
    int       nx, ny, nz;    /* interior cells                 */
    double    h;             /* grid spacing                   */
    float    *u;             /* current iterate                */
    float    *u_new;         /* next iterate                   */
    int      *material;      /* material id of each cell       */
    Material *materials;     /* material table                 */
    int       nmaterials;
};

static inline long idx(const Grid *g, int i, int j, int k)
{
    return ((long) i * (g->ny + 2) + j) * (g->nz + 2) + k;
}

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
