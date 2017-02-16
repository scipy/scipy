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
typedef double (*solver_type)(callback_type, double, double, double, double, int,default_parameters*);

extern double bisect(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params);
extern double ridder(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params);
extern double brenth(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params);
extern double brentq(callback_type f, double xa, double xb, double xtol, double rtol, int iter, default_parameters *params);

#endif
