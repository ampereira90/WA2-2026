/*
 * ACA - Work Assignment 2 - 3D Jacobi
 * Common infrastructure: command line, timing, input dataset, result output.
 *
 * You normally do NOT need to change this file. If you do, keep the input
 * and output file formats unchanged so that tools/gen_input.py and
 * tools/compare keep working with your versions.
 */
#ifndef WA2_COMMON_H
#define WA2_COMMON_H

#include <cstdio>
#include <cstdint>
#include <vector>

/* ------------------------------------------------------------------ */
/* Command line:   <variant> [options] INPUT OUTPUT                    */
/* ------------------------------------------------------------------ */
struct Options {
    const char *input;           /* dataset produced by tools/gen_input.py */
    const char *output;          /* result file (for tools/compare)        */
    int         max_iters;       /* number of Jacobi iterations            */
    double      tol;             /* stop when ||u_new - u||_2 < tol (0=off) */
    int         print_every;     /* print convergence every N iters (0=off)*/
    int         nx, ny, nz;      /* filled in from the input file          */
};

void parse_options(int argc, char **argv, Options *opt,
                   const char *variant, const char *description);

/* ------------------------------------------------------------------ */
/* Timing                                                              */
/* ------------------------------------------------------------------ */
double wall_time(void);          /* seconds, monotonic clock */

/* ------------------------------------------------------------------ */
/* Input dataset                                                       */
/*                                                                     */
/* Little-endian binary file, written by tools/gen_input.py:           */
/*   InputHeader                                (56 bytes)             */
/*   InputField[nfields]                        (32 bytes each)        */
/*   MaterialRecord[nmaterials]                 (56 bytes each)        */
/*   field data, at the offsets given in InputField                    */
/* Every field covers the whole grid INCLUDING the halo:               */
/*   (nx+2) x (ny+2) x (nz+2) values, k fastest, then j, then i.       */
/* Halo values are the Dirichlet boundary values.                      */
/* ------------------------------------------------------------------ */
enum { DTYPE_FLOAT32 = 0, DTYPE_UINT8 = 1 };

struct InputHeader {
    char    magic[8];            /* "ACAWA2IN"                        */
    char    variant[16];         /* "poisson7", "box27" or "varcoef"  */
    int32_t nx, ny, nz;          /* interior points                   */
    int32_t nfields;
    double  h;                   /* grid spacing                      */
    int32_t nmaterials;
    int32_t reserved;
};

struct InputField {
    char    name[16];            /* e.g. "u", "f", "material"         */
    int32_t dtype;               /* DTYPE_FLOAT32 or DTYPE_UINT8      */
    int32_t reserved;
    int64_t offset;              /* byte offset from start of file    */
};

struct MaterialRecord {
    char   name[24];
    double conductivity;
    double density;
    double heat_capacity;
    double heat_source;
};

class InputFile {
public:
    InputFile(const char *path, const char *expected_variant);
    ~InputFile();
    int    nx() const { return hdr_.nx; }
    int    ny() const { return hdr_.ny; }
    int    nz() const { return hdr_.nz; }
    double h()  const { return hdr_.h; }
    int    nmaterials() const { return hdr_.nmaterials; }
    const MaterialRecord &material(int m) const { return materials_[m]; }
    /* number of values in one i-plane, (ny+2)*(nz+2) */
    long   plane_size() const { return (long)(hdr_.ny + 2) * (hdr_.nz + 2); }
    /* read 'count' consecutive i-planes of a field, starting at plane i
       (0 <= i, i+count <= nx+2), into buf (float* or uint8_t*) */
    void   read_planes(const char *field, int i, int count, void *buf);
private:
    FILE                       *fp_;
    const char                 *path_;
    InputHeader                 hdr_;
    std::vector<InputField>     fields_;
    std::vector<MaterialRecord> materials_;
};

/* ------------------------------------------------------------------ */
/* Result output                                                       */
/*                                                                     */
/* The result is the interior of the final solution u (halo excluded), */
/* written row by row with index k varying fastest, then j, then i:    */
/*      for i in 1..nx, for j in 1..ny: write u(i, j, 1..nz)           */
/* File layout: 40-byte header followed by nx*ny*nz float32 values.    */
/* ------------------------------------------------------------------ */
struct ResultHeader {
    char    magic[8];            /* "ACAWA2\0\0"                     */
    char    variant[16];         /* e.g. "poisson7"                  */
    int32_t nx, ny, nz, iters;
};

class ResultSink {
public:
    ResultSink(const char *path, const char *variant,
               int nx, int ny, int nz, int iters);
    ~ResultSink();
    /* append one row of nz consecutive values u(i, j, 1..nz) */
    void   row(const float *values, int n);
    double checksum() const { return sum_; }
    double checksum_abs() const { return sum_abs_; }
    long   count() const { return count_; }
private:
    FILE  *fp_;
    double sum_, sum_abs_;
    long   count_;
};

/* ------------------------------------------------------------------ */
/* Report                                                              */
/* ------------------------------------------------------------------ */
void print_report(const char *variant, const Options *opt, int iters_done,
                  double residual, double t_load, double t_solve,
                  const ResultSink *res);

#endif
