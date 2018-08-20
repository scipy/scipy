/* Written by Charles Harris charles.harris@sdl.usu.edu */

/* Modified to not depend on Python everywhere by Travis Oliphant.
 */

#ifndef ZEROS_H
#define ZEROS_H

typedef struct {
    int funcalls;
    int iterations;
    int error_num;
} default_parameters;

#define SIGNERR -1
#define CONVERR -2

typedef double (*callback_type)(double,void*);

extern double newton(callback_type f, double p0, callback_type fprime, default_parameters *params, double tol, int maxiter);
extern double secant(callback_type f, double p0, default_parameters *params, double tol, int maxiter);
extern double halley(callback_type f, double p0, callback_type fprime, default_parameters *params, double tol, int maxiter, callback_type fprime2);

#endif
