/*
 * This file is a collection (more can be added) of wrappers around some
 * CDFLIB Fortran algorithms.
 */

/*
 * Notice q and p are used in reverse from their meanings in distributions.py
 */

#include "cdf_wrappers.h"

#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif

/*
 * Call macros for different numbers of distribution parameters.
 */

#define CDFLIB_CALL2(func, name, a, b, result, return_bound)            \
    if (npy_isnan(p) || npy_isnan(q) || npy_isnan(a) || npy_isnan(b) || \
        npy_isnan(bound)) {                                             \
        return NPY_NAN;                                                 \
    }                                                                   \
    func(&which, &p, &q, &a, &b, &status, &bound);                      \
    return get_result(name, status, bound, result, return_bound)

#define CDFLIB_CALL3(func, name, a, b, c, result, return_bound)         \
    if (npy_isnan(p) || npy_isnan(q) || npy_isnan(a) || npy_isnan(b) || \
        npy_isnan(c) || npy_isnan(bound)) {                             \
        return NPY_NAN;                                                 \
    }                                                                   \
    func(&which, &p, &q, &a, &b, &c, &status, &bound);                  \
    return get_result(name, status, bound, result, return_bound)

#define CDFLIB_CALL4(func, name, a, b, c, d, result, return_bound)      \
    if (npy_isnan(p) || npy_isnan(q) || npy_isnan(a) || npy_isnan(b) || \
        npy_isnan(c) || npy_isnan(d) || npy_isnan(bound)) {             \
        return NPY_NAN;                                                 \
    }                                                                   \
    func(&which, &p, &q, &a, &b, &c, &d, &status, &bound);              \
    return get_result(name, status, bound, result, return_bound)


/* Return nan on status==1,2 */
#define NO_RETURN_BOUND 0
/* Return bound on status==1,2 */
#define RETURN_BOUND 1


/*
 * Return status checking function
 */

static double get_result(char *name, int status, double bound, double result, int return_bound) {
  if (status < 0) {
      sf_error(name, SF_ERROR_ARG, "(Fortran) input parameter %d is out of range", (-status));
  }
  else {
    switch (status) {
    case 0:
      /* no error */
      return result;
    case 1:
      sf_error(name, SF_ERROR_OTHER, "Answer appears to be lower than lowest search bound (%g)", bound);
      if (return_bound) {
          return bound;
      }
      break;
    case 2:
      sf_error(name, SF_ERROR_OTHER, "Answer appears to be higher than highest search bound (%g)", bound);
      if (return_bound) {
          return bound;
      }
      break;
    case 3:
    case 4:
      sf_error(name, SF_ERROR_OTHER, "Two parameters that should sum to 1.0 do not");
      break;
    case 10:
      sf_error(name, SF_ERROR_OTHER, "Computational error");
      break;
    default:
      sf_error(name, SF_ERROR_OTHER, "Unknown error");
    }
  }
  return NPY_NAN;
}


extern void F_FUNC(cdfbet,CDFBET)(int*,double*,double*,double*,double*,double*,double*,int*,double*);

double cdfbet3_wrap(double p, double b, double x) {
  int which=3;
  double q=1.0-p, y=1.0-x, a=0, bound=0;
  int status=10;

  CDFLIB_CALL4(F_FUNC(cdfbet,CDFBET), "btdtria",
               x, y, a, b,
               a, RETURN_BOUND);
}

double cdfbet4_wrap(double a, double p, double x) {
  int which=4;
  double q=1.0-p, y=1.0-x, b=0, bound=0;
  int status=10;

  CDFLIB_CALL4(F_FUNC(cdfbet,CDFBET), "btdtrib",
               x, y, a, b,
               b, RETURN_BOUND);
}


extern void F_FUNC(cdfbin,CDFBIN)(int*,double*,double*,double*,double*,double*,double*,int*,double*);

double cdfbin2_wrap(double p, double xn, double pr) {
  int which=2;
  double q=1.0-p, s=0, ompr=1.0-pr, bound=0;
  int status=10;

  CDFLIB_CALL4(F_FUNC(cdfbin,CDFBIN), "bdtrik",
               s, xn, pr, ompr,
               s, RETURN_BOUND);
}

double cdfbin3_wrap(double s, double p, double pr) {
  int which=3;
  double q=1.0-p, xn=0, ompr=1.0-pr, bound=0;
  int status=10;

  CDFLIB_CALL4(F_FUNC(cdfbin,CDFBIN), "bdtrin",
               s, xn, pr, ompr,
               xn, RETURN_BOUND);
}

extern void F_FUNC(cdfchi,CDFCHI)(int*,double*,double*,double*,double*,int*,double*);
double cdfchi3_wrap(double p, double x){
  int which=3;
  double q=1.0-p, df=0, bound=0;
  int status=10;

  CDFLIB_CALL2(F_FUNC(cdfchi,CDFCHI), "chdtriv",
               x, df,
               df, RETURN_BOUND);
}

extern void F_FUNC(cdfchn,CDFCHN)(int*,double*,double*,double*,double*,double*,int*,double*);
double cdfchn1_wrap(double x, double df, double nc) {
  int which=1;
  double q=0, p=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdfchn,CDFCHN), "chndtr",
               x, df, nc,
               p, RETURN_BOUND);
}

double cdfchn2_wrap(double p, double df, double nc) {
  int which=2;
  double q=1.0-p, x=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdfchn,CDFCHN), "chndtrix",
               x, df, nc,
               x, NO_RETURN_BOUND);
}

double cdfchn3_wrap(double x, double p, double nc) {
  int which=3;
  double q=1.0-p, df=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdfchn,CDFCHN), "chndtridf",
               x, df, nc,
               df, RETURN_BOUND);
}

double cdfchn4_wrap(double x, double df, double p) {
  int which=4;
  double q=1.0-p, nc=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdfchn,CDFCHN), "chndtrinc",
               x, df, nc,
               nc, RETURN_BOUND);
}

extern void F_FUNC(cdff,CDFF)(int*,double*,double*,double*,double*,double*,int*,double*);
/*
double cdff1_wrap(double dfn, double dfd, double f) {
  int which=1;
  double q, p, bound;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdff,CDFF), "fdtr",
               f, dfn, dfd,
               p, NO_RETURN_BOUND);
}

double cdff2_wrap(double dfn, double dfd, double p) {
  int which=2;
  double q=1.0-p, f, bound;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdff,CDFF), "fdtri",
               f, dfn, dfd,
               f, NO_RETURN_BOUND);
}
*/

/* This seem to give some trouble.  No idea why... */
double cdff3_wrap(double p, double dfd, double f) {
  int which=3;
  double q=1.0-p, dfn=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdff,CDFF), "fdtridfn",
               f, dfn, dfd,
               dfn, RETURN_BOUND);
}

double cdff4_wrap(double dfn, double p, double f) {
  int which=4;
  double q=1.0-p, dfd=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdff,CDFF), "fdtridfd",
               f, dfn, dfd,
               dfd, RETURN_BOUND);
}


extern void F_FUNC(cdffnc,CDFFNC)(int*,double*,double*,double*,double*,double*,double*,int*,double*);
double cdffnc1_wrap(double dfn, double dfd, double nc, double f) {
  int which=1;
  double q=0, p=0, bound=0;
  int status=10;

  CDFLIB_CALL4(F_FUNC(cdffnc,CDFFNC), "ncfdtr",
               f, dfn, dfd, nc,
               p, NO_RETURN_BOUND);
}

double cdffnc2_wrap(double dfn, double dfd, double nc, double p) {
  int which=2;
  double q=1.0-p, f=0, bound=0;
  int status=10;

  CDFLIB_CALL4(F_FUNC(cdffnc,CDFFNC), "ncfdtri",
               f, dfn, dfd, nc,
               f, RETURN_BOUND);
}


double cdffnc3_wrap(double p, double dfd, double nc, double f) {
  int which=3;
  double q=1.0-p, dfn=0, bound=0;
  int status=10;

  CDFLIB_CALL4(F_FUNC(cdffnc,CDFFNC), "ncfdtridfn",
               f, dfn, dfd, nc,
               dfn, RETURN_BOUND);
}

double cdffnc4_wrap(double dfn, double p, double nc, double f) {
  int which=4;
  double q=1.0-p, dfd=0, bound=0;
  int status=10;

  CDFLIB_CALL4(F_FUNC(cdffnc,CDFFNC), "ncfdtridfd",
               f, dfn, dfd, nc,
               dfd, RETURN_BOUND);
}

double cdffnc5_wrap(double dfn, double dfd, double p, double f) {
  int which=5;
  double q=1.0-p, nc=0, bound=0;
  int status=10;

  CDFLIB_CALL4(F_FUNC(cdffnc,CDFFNC), "ncfdtrinc",
               f, dfn, dfd, nc,
               nc, RETURN_BOUND);
}

/* scl == a in gdtr
   shp == b in gdtr
*/ 
extern void F_FUNC(cdfgam,CDFGAM)(int*,double*,double*,double*,double*,double*,int*,double*);
double cdfgam1_wrap(double scl, double shp, double x) {
  int which=1;
  double q=0, p=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdfgam,CDFGAM), "gdtr",
               x, shp, scl,
               p, NO_RETURN_BOUND);
}

double cdfgam2_wrap(double scl, double shp, double p) {
  int which=2;
  double q=1.0-p, x=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdfgam,CDFGAM), "gdtrix",
               x, shp, scl,
               x, RETURN_BOUND);
}

double cdfgam3_wrap(double scl, double p, double x) {
  int which=3;
  double q=1.0-p, shp=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdfgam,CDFGAM), "gdtrib",
               x, shp, scl,
               shp, RETURN_BOUND);
}

double cdfgam4_wrap(double p, double shp, double x) {
  int which=4;
  double q=1.0-p, scl=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdfgam,CDFGAM), "gdtria",
               x, shp, scl,
               scl, RETURN_BOUND);
}

extern void F_FUNC(cdfnbn,CDFNBN)(int*,double*,double*,double*,double*,double*,double*,int*,double*);
double cdfnbn2_wrap(double p, double xn, double pr) {
  int which=2;
  double q=1.0-p, s=0, ompr=1.0-pr, bound=0;
  int status=10;

  CDFLIB_CALL4(F_FUNC(cdfnbn,CDFNBN), "nbdtrik",
               s, xn, pr, ompr,
               s, RETURN_BOUND);
}

double cdfnbn3_wrap(double s, double p, double pr) {
  int which=3;
  double q=1.0-p, xn=0, ompr=1.0-pr, bound=0;
  int status=10;

  CDFLIB_CALL4(F_FUNC(cdfnbn,CDFNBN), "nbdtrin",
               s, xn, pr, ompr,
               xn, RETURN_BOUND);
}

extern void F_FUNC(cdfnor,CDFNOR)(int*,double*,double*,double*,double*,double*,int*,double*);
double cdfnor3_wrap(double p, double std, double x) {
  int which=3;
  double q=1.0-p, mn=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdfnor,CDFNOR), "nrdtrimn",
               x, mn, std,
               mn, RETURN_BOUND);
}

double cdfnor4_wrap(double mn, double p, double x) {
  int which=4;
  double q=1.0-p, std=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdfnor,CDFNOR), "nrdtrisd",
               x, mn, std,
               std, RETURN_BOUND);
}

extern void F_FUNC(cdfpoi,CDFPOI)(int*,double*,double*,double*,double*,int*,double*);
double cdfpoi2_wrap(double p, double xlam){
  int which=2;
  double q=1.0-p, s=0, bound=0;
  int status=10;

  CDFLIB_CALL2(F_FUNC(cdfpoi,CDFPOI), "pdtrik",
               s, xlam,
               s, RETURN_BOUND);
}

extern void F_FUNC(cdft,CDFT)(int*,double*,double*,double*,double*,int*,double*);
double cdft1_wrap(double df, double t){
  int which=1;
  double q=0, p=0, bound=0;
  int status=10;

  CDFLIB_CALL2(F_FUNC(cdft,CDFT), "stdtr",
               t, df,
               p, NO_RETURN_BOUND);
}

double cdft2_wrap(double df, double p){
  int which=2;
  double q=1.0-p, t=0, bound=0;
  int status=10;

  CDFLIB_CALL2(F_FUNC(cdft,CDFT), "stdtrit",
               t, df,
               t, RETURN_BOUND);
}

double cdft3_wrap(double p, double t){
  int which=3;
  double q=1.0-p, df=0, bound=0;
  int status=10;

  CDFLIB_CALL2(F_FUNC(cdft,CDFT), "stdtridf",
               t, df,
               df, RETURN_BOUND);
}

extern void F_FUNC(cdftnc,CDFTNC)(int*,double*,double*,double*,double*,double*,int*,double*);
double cdftnc1_wrap(double df, double nc, double t) {
  int which=1;
  double q=0, p=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdftnc,CDFTNC), "nctdtr",
               t, df, nc,
               p, RETURN_BOUND);
}

double cdftnc2_wrap(double df, double nc, double p) {
  int which=2;
  double q=1.0-p, t=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdftnc,CDFTNC), "nctdtrit",
               t, df, nc,
               t, RETURN_BOUND);
}

double cdftnc3_wrap(double p, double nc, double t) {
  int which=3;
  double q=1.0-p, df=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdftnc,CDFTNC), "nctdtridf",
               t, df, nc,
               df, RETURN_BOUND);
}

double cdftnc4_wrap(double df, double p, double t) {
  int which=4;
  double q=1.0-p, nc=0, bound=0;
  int status=10;

  CDFLIB_CALL3(F_FUNC(cdftnc,CDFTNC), "nctdtrinc",
               t, df, nc,
               nc, RETURN_BOUND);
}
