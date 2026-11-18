/*
 * ACA - Work Assignment 2 - common infrastructure (see common.h)
 */
#include "common.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>

/* ------------------------------------------------------------------ */
/* Command line                                                        */
/* ------------------------------------------------------------------ */
static void usage(const char *prog, const char *variant, const char *desc)
{
    std::fprintf(stderr,
        "%s - %s\n\n"
        "usage: %s [options] INPUT OUTPUT\n"
        "  INPUT       dataset created with tools/gen_input.py %s ...\n"
        "  OUTPUT      file for the final solution (check it with tools/compare)\n"
        "  -i ITERS    number of Jacobi iterations (default 100)\n"
        "  -t TOL      stop early when ||u_new-u||_2 < TOL (default 0 = never)\n"
        "  -p EVERY    print convergence every EVERY iterations (default 0)\n"
        "  -h          this help\n",
        variant, desc, prog, variant);
}

void parse_options(int argc, char **argv, Options *opt,
                   const char *variant, const char *desc)
{
    opt->input       = NULL;
    opt->output      = NULL;
    opt->max_iters   = 100;
    opt->tol         = 0.0;
    opt->print_every = 0;
    opt->nx = opt->ny = opt->nz = 0;

    int npos = 0;
    for (int a = 1; a < argc; a++) {
        const char *s = argv[a];
        bool has_arg = (a + 1 < argc);
        if      (!std::strcmp(s, "-i") && has_arg) { opt->max_iters = std::atoi(argv[++a]); }
        else if (!std::strcmp(s, "-t") && has_arg) { opt->tol = std::atof(argv[++a]); }
        else if (!std::strcmp(s, "-p") && has_arg) { opt->print_every = std::atoi(argv[++a]); }
        else if (s[0] == '-') { usage(argv[0], variant, desc); std::exit(!std::strcmp(s, "-h") ? 0 : 1); }
        else if (npos == 0) { opt->input  = s; npos++; }
        else if (npos == 1) { opt->output = s; npos++; }
        else { usage(argv[0], variant, desc); std::exit(1); }
    }
    if (npos != 2) { usage(argv[0], variant, desc); std::exit(1); }
    if (opt->max_iters < 1) { std::fprintf(stderr, "invalid iteration count\n"); std::exit(1); }
}

double wall_time(void)
{
    using namespace std::chrono;
    return duration_cast<duration<double> >(
               steady_clock::now().time_since_epoch()).count();
}

/* ------------------------------------------------------------------ */
/* Input dataset                                                       */
/* ------------------------------------------------------------------ */
static void input_error(const char *path, const char *msg)
{
    std::fprintf(stderr, "%s: %s\n", path, msg);
    std::exit(1);
}

InputFile::InputFile(const char *path, const char *expected_variant)
    : fp_(NULL), path_(path)
{
    fp_ = std::fopen(path, "rb");
    if (!fp_) { std::perror(path); std::exit(1); }
    if (std::fread(&hdr_, sizeof hdr_, 1, fp_) != 1 ||
        std::memcmp(hdr_.magic, "ACAWA2IN", 8) != 0)
        input_error(path, "not a WA2 input file (create it with tools/gen_input.py)");
    hdr_.variant[sizeof hdr_.variant - 1] = '\0';
    if (std::strcmp(hdr_.variant, expected_variant) != 0) {
        std::fprintf(stderr, "%s: dataset is for variant '%s', not '%s'\n",
                     path, hdr_.variant, expected_variant);
        std::exit(1);
    }
    if (hdr_.nx < 1 || hdr_.ny < 1 || hdr_.nz < 1 || hdr_.nfields < 1 ||
        hdr_.nfields > 64 || hdr_.nmaterials < 0 || hdr_.nmaterials > 256)
        input_error(path, "corrupted header");
    fields_.resize(hdr_.nfields);
    if (std::fread(&fields_[0], sizeof(InputField), fields_.size(), fp_) != fields_.size())
        input_error(path, "truncated field table");
    for (size_t f = 0; f < fields_.size(); f++)
        fields_[f].name[sizeof fields_[f].name - 1] = '\0';
    if (hdr_.nmaterials > 0) {
        materials_.resize(hdr_.nmaterials);
        if (std::fread(&materials_[0], sizeof(MaterialRecord), materials_.size(), fp_) != materials_.size())
            input_error(path, "truncated material table");
        for (size_t m = 0; m < materials_.size(); m++)
            materials_[m].name[sizeof materials_[m].name - 1] = '\0';
    }
}

InputFile::~InputFile()
{
    if (fp_) std::fclose(fp_);
}

void InputFile::read_planes(const char *field, int i, int count, void *buf)
{
    const InputField *fd = NULL;
    for (size_t f = 0; f < fields_.size(); f++)
        if (!std::strcmp(fields_[f].name, field)) fd = &fields_[f];
    if (!fd) {
        std::fprintf(stderr, "%s: field '%s' not found\n", path_, field);
        std::exit(1);
    }
    if (i < 0 || count < 0 || i + count > hdr_.nx + 2)
        input_error(path_, "plane index out of range");
    const size_t elem = (fd->dtype == DTYPE_UINT8) ? 1 : 4;
    const size_t n    = (size_t) plane_size() * (size_t) count;
    const long   off  = (long)(fd->offset + (int64_t) i * plane_size() * (int64_t) elem);
    if (std::fseek(fp_, off, SEEK_SET) != 0 || std::fread(buf, elem, n, fp_) != n)
        input_error(path_, "read error (truncated file?)");
}

/* ------------------------------------------------------------------ */
/* Result output                                                       */
/* ------------------------------------------------------------------ */
ResultSink::ResultSink(const char *path, const char *variant,
                       int nx, int ny, int nz, int iters)
    : fp_(NULL), sum_(0.0), sum_abs_(0.0), count_(0)
{
    if (!path) return;
    fp_ = std::fopen(path, "wb");
    if (!fp_) { std::perror(path); std::exit(1); }
    ResultHeader h;
    std::memset(&h, 0, sizeof h);
    std::memcpy(h.magic, "ACAWA2", 6);
    std::strncpy(h.variant, variant, sizeof h.variant - 1);
    h.nx = nx; h.ny = ny; h.nz = nz; h.iters = iters;
    if (std::fwrite(&h, sizeof h, 1, fp_) != 1) { std::perror(path); std::exit(1); }
}

ResultSink::~ResultSink()
{
    if (fp_) std::fclose(fp_);
}

void ResultSink::row(const float *v, int n)
{
    for (int k = 0; k < n; k++) { sum_ += v[k]; sum_abs_ += std::fabs(v[k]); }
    count_ += n;
    if (fp_ && std::fwrite(v, sizeof(float), (size_t)n, fp_) != (size_t)n) {
        std::perror("write"); std::exit(1);
    }
}

/* ------------------------------------------------------------------ */
/* Report                                                              */
/* ------------------------------------------------------------------ */
void print_report(const char *variant, const Options *opt, int iters_done,
                  double residual, double t_load, double t_solve,
                  const ResultSink *res)
{
    double points = (double)opt->nx * opt->ny * opt->nz;
    double mlups  = points * iters_done / t_solve / 1e6;
    std::printf("variant        : %s\n", variant);
    std::printf("input          : %s\n", opt->input);
    std::printf("output         : %s\n", opt->output);
    std::printf("grid (interior): %d x %d x %d  (%.3e points)\n",
                opt->nx, opt->ny, opt->nz, points);
    std::printf("iterations     : %d\n", iters_done);
    std::printf("final ||du||_2 : %.6e\n", residual);
    std::printf("load time   (s): %.6f  (allocation + reading the input)\n", t_load);
    std::printf("solve time  (s): %.6f\n", t_solve);
    std::printf("time/iter  (ms): %.6f\n", 1e3 * t_solve / iters_done);
    std::printf("MLUP/s         : %.3f\n", mlups);
    std::printf("checksum       : %.10e  (sum |u| = %.10e)\n",
                res->checksum(), res->checksum_abs());
    /* machine-readable line for building result tables */
    std::printf("CSV,%s,%d,%d,%d,%d,%.6e,%.6f,%.6f,%.3f,%.10e\n",
                variant, opt->nx, opt->ny, opt->nz, iters_done, residual,
                t_load, t_solve, mlups, res->checksum());
}
