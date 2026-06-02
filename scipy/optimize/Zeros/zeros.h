/* Written by Charles Harris charles.harris@sdl.usu.edu */

/* Modified to not depend on Python everywhere by Travis Oliphant.
 */

#ifndef ZEROS_H
#define ZEROS_H

typedef struct {
    int funcalls;
    int iterations;
    int error_num;
} scipy_zeros_info;


/* Must agree with _ECONVERGED, _ESIGNERR, _ECONVERR  in zeros.py */
#define CONVERGED 0
#define SIGNERR -1
#define CONVERR -2
#define EVALUEERR -3
#define INPROGRESS 1

typedef double (*callback_type)(double, void*);
typedef double (*solver_type)(callback_type, double, double, double, double,
                              int, void *, scipy_zeros_info*);

extern double bisect(callback_type f, double xa, double xb, double xtol,
                     double rtol, int iter, void *func_data_param,
                     scipy_zeros_info *solver_stats);
extern double ridder(callback_type f, double xa, double xb, double xtol,
                     double rtol, int iter, void *func_data_param,
                     scipy_zeros_info *solver_stats);
extern double brenth(callback_type f, double xa, double xb, double xtol,
                     double rtol, int iter, void *func_data_param,
                     scipy_zeros_info *solver_stats);
extern double brentq(callback_type f, double xa, double xb, double xtol,
                     double rtol, int iter, void *func_data_param,
                     scipy_zeros_info *solver_stats);

#endif
