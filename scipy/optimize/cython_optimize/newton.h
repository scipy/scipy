#ifndef NEWTON_H
#define NEWTON_H

#include "../Zeros/zeros.h"

#define SIGNERR -1
#define CONVERR -2

typedef double (*callback_type)(double, void*);

extern double newton(callback_type f, double p0, callback_type fprime, default_parameters *params, double tol, int maxiter);
extern double secant(callback_type f, double p0, default_parameters *params, double tol, int maxiter);
extern double halley(callback_type f, double p0, callback_type fprime, default_parameters *params, double tol, int maxiter, callback_type fprime2);

#endif
