#ifndef DIRECT_H
#define DIRECT_H

#include "Python.h"
#include <math.h>
#include <stdio.h>

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

typedef enum {
     DIRECT_ORIGINAL, DIRECT_GABLONSKY
} direct_algorithm;

typedef enum {
     DIRECT_INVALID_BOUNDS = -1,
     DIRECT_MAXFEVAL_TOOBIG = -2,
     DIRECT_INIT_FAILED = -3,
     DIRECT_SAMPLEPOINTS_FAILED = -4,
     DIRECT_SAMPLE_FAILED = -5,
     DIRECT_MAXLEVELS_REACHED = -6,
     DIRECT_MAXFEVAL_EXCEEDED = 1,
     DIRECT_MAXITER_EXCEEDED = 2,
     DIRECT_GLOBAL_FOUND = 3,
     DIRECT_VOLTOL = 4,
     DIRECT_SIGMATOL = 5,

     DIRECT_OUT_OF_MEMORY = -100,
     DIRECT_INVALID_ARGS = -101,
     DIRECT_FORCED_STOP = -102
} direct_return_code;

typedef struct {
     int numfunc;
     int numiter;
} direct_return_info;

#define DIRECT_UNKNOWN_FGLOBAL (-HUGE_VAL)
#define DIRECT_UNKNOWN_FGLOBAL_RELTOL (0.0)

extern PyObject* direct_optimize(
    PyObject* f, double *x, PyObject *x_seq, PyObject *args,
    int dimension,
    const double *lower_bounds, const double *upper_bounds,
    double *minf,
    int max_feval, int max_iter,
    double magic_eps, double magic_eps_abs,
    double volume_reltol, double sigma_reltol,
    int *force_stop,
    double fglobal,
    double fglobal_reltol,
    FILE *logfile,
    direct_algorithm algorithm,
    direct_return_info *info,
    direct_return_code *ret_code,
    PyObject* callback);

#ifdef __cplusplus
}  /* extern "C" */
#endif /* __cplusplus */

#endif /* DIRECT_H */
