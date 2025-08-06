#ifndef __SLSQPLIB_H
#define __SLSQPLIB_H

#include <math.h>
#include "blaslapack_declarations.h"
#include "nnls.h"

// The SLSQP_vars struct holds the state of the algorithm and passed to Python and back such that it is thread-safe.
struct SLSQP_vars {
    double acc, alpha, f0, gs, h1, h2, h3, h4, t, t0, tol;
    int exact, inconsistent, reset, iter, itermax, line, m, meq, mode, n;
};

void __slsqp_body(struct SLSQP_vars* S, double* funx, double* gradx, double* C, double* d, double* sol, double* mult, double* xl, double* xu, double* buffer, int* indices);

#endif // __SLSQPLIB_H
