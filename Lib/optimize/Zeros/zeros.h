/* Written by Charles Harris charles.harris@sdl.usu.edu */

/* Modified to not depend on Python everywhere by Travis Oliphant.
 */


#ifndef ZEROS_H
#define ZEROS_H

#define ZEROS_PARAM_HEAD int funcalls; int iterations

typedef struct {
    ZEROS_PARAM_HEAD;
} default_parameters;

extern int zeros_error_num;

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define SIGN(a)   ((a) > 0.0 ? 1.0 : -1.0)
#define ERROR(num,val) zeros_error_num=(num); return (val)
#define SIGNERR -1
#define CONVERR -2

double bisect(double (*f)(double, void*), double xa, double xb, double xtol, double rtol, int iter, default_parameters *params);
double ridder(double (*f)(double, void*), double xa, double xb, double xtol, double rtol, int iter, default_parameters *params);
double brenth(double (*f)(double, void*), double xa, double xb, double xtol, double rtol, int iter, default_parameters *params);
double brentq(double (*f)(double, void*), double xa, double xb, double xtol, double rtol, int iter, default_parameters *params);

#endif
