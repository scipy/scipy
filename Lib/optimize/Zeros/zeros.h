/* Written by Charles Harris charles.harris@sdl.usu.edu */

/* Modified to not depend on Python everywhere by Travis Oliphant.
 */


#ifndef ZEROS_H
#define ZEROS_H

#define ZEROS_PARAM_HEAD int funcalls; int iterations; int error_num

typedef struct {
    ZEROS_PARAM_HEAD;
} default_parameters;

static double dminarg1,dminarg2;
#define DMIN(a,b) (dminarg1=(a),dminarg2=(b),(dminarg1) < (dminarg2) ?\
        (dminarg1) : (dminarg2))

#define SIGN(a)   ((a) > 0.0 ? 1.0 : -1.0)
#define ERROR(params,num,val) (params)->error_num=(num); return (val)
#define SIGNERR -1
#define CONVERR -2

extern double bisect(double (*f)(double, void*), double xa, double xb, double xtol, double rtol, int iter, default_parameters *params);
extern double ridder(double (*f)(double, void*), double xa, double xb, double xtol, double rtol, int iter, default_parameters *params);
extern double brenth(double (*f)(double, void*), double xa, double xb, double xtol, double rtol, int iter, default_parameters *params);
extern double brentq(double (*f)(double, void*), double xa, double xb, double xtol, double rtol, int iter, default_parameters *params);

extern double fabs(double);
extern double sqrt(double);

#endif
