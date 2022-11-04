/* This file is a collection of wrappers around the
 *  Amos Fortran library of functions that take complex
 *  variables (see www.netlib.org) so that they can be called from
 *  the cephes library of corresponding name but work with complex
 *  arguments.
 */

#ifndef _CDF_WRAPPERS_H
#define _CDF_WRAPPERS_H
#ifndef _AMOS_WRAPPERS_H
#include "Python.h"
#endif

#include "sf_error.h"

#include <numpy/npy_math.h>

extern double cdfbet3_wrap(double p, double x, double b);
extern double cdfbet4_wrap(double p, double x, double a);

extern double cdfbin2_wrap(double p, double xn, double pr);
extern double cdfbin3_wrap(double p, double s, double pr);

extern double cdfchi3_wrap(double p, double x);

extern double cdfchn1_wrap(double x, double df, double nc);
extern double cdfchn2_wrap(double p, double df, double nc);
extern double cdfchn3_wrap(double p, double x, double nc);
extern double cdfchn4_wrap(double p, double x, double df);

extern double cdff3_wrap(double p, double f, double dfd);
extern double cdff4_wrap(double p, double f, double dfn);

extern double cdffnc1_wrap(double f, double dfn, double dfd, double nc);
extern double cdffnc2_wrap(double p, double dfn, double dfd, double nc);
extern double cdffnc3_wrap(double p, double f, double dfd, double nc);
extern double cdffnc4_wrap(double p, double f, double dfn, double nc);
extern double cdffnc5_wrap(double p, double f, double dfn, double dfd);

extern double cdfgam1_wrap(double p, double x, double scl);
extern double cdfgam2_wrap(double p, double x, double shp);
extern double cdfgam3_wrap(double p, double x, double scl);
extern double cdfgam4_wrap(double p, double x, double shp);

extern double cdfnbn2_wrap(double p, double xn, double pr);
extern double cdfnbn3_wrap(double p, double s, double pr);

extern double cdfnor3_wrap(double p, double x, double std);
extern double cdfnor4_wrap(double p, double x, double mn);

extern double cdfpoi2_wrap(double p, double xlam);

extern double cdft1_wrap(double p, double t);
extern double cdft2_wrap(double p, double t);
extern double cdft3_wrap(double p, double t);

extern double cdftnc1_wrap(double df, double nc, double t);
extern double cdftnc2_wrap(double df, double nc, double p);
extern double cdftnc3_wrap(double p, double nc, double t);
extern double cdftnc4_wrap(double df, double p, double t);

extern double tukeylambdacdf(double x, double lambda);
#endif
