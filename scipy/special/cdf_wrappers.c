/* This file is a collection (more can be added) of wrappers around some
 *  CDF Fortran algorithms, so that they can be called from
 *  cephesmodule.so
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

/* This must be linked with fortran
 */

/* Notice q and p are used in reverse from their meanings in distributions.py
 */

static void show_error(char *func, int status, int bound) {
  /* show_error message */

  if (status < 0) {
      sf_error(func, SF_ERROR_ARG, "(Fortran) input parameter %d is out of range", (-status));
  }
  else {
    switch (status) {
    case 1:
      sf_error(func, SF_ERROR_OTHER, "Answer appears to be lower than lowest search bound (%d)", bound);
      break;
    case 2:
      sf_error(func, SF_ERROR_OTHER, "Answer appears to be higher than highest search bound (%d)", bound);
      break;
    case 3:
    case 4:
      sf_error(func, SF_ERROR_OTHER, "Two parameters that should sum to 1.0 do not");
      break;
    case 10:
      sf_error(func, SF_ERROR_OTHER, "Computational error");
      break;
    default:
      sf_error(func, SF_ERROR_OTHER, "Unknown error");
    }
  }
}

extern void F_FUNC(cdfbet,CDFBET)(int*,double*,double*,double*,double*,double*,double*,int*,double*);

double cdfbet3_wrap(double p, double b, double x) {
  int which=3;
  double q=1.0-p, y=1.0-x, a, bound;
  int status;  
  
  F_FUNC(cdfbet,CDFBET)(&which, &p, &q, &x, &y, &a, &b, &status, &bound);
  if (status) {
    show_error("cdfbet3", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return a;
}

double cdfbet4_wrap(double a, double p, double x) {
  int which=4;
  double q=1.0-p, y=1.0-x, b, bound;
  int status;  
  
  F_FUNC(cdfbet,CDFBET)(&which, &p, &q, &x, &y, &a, &b, &status, &bound);
  if (status) {
    show_error("cdfbet4", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return b;
}


extern void F_FUNC(cdfbin,CDFBIN)(int*,double*,double*,double*,double*,double*,double*,int*,double*);

double cdfbin2_wrap(double p, double xn, double pr) {
  int which=2;
  double q=1.0-p, s, ompr=1.0-pr, bound;
  int status;  
  
  F_FUNC(cdfbin,CDFBIN)(&which, &p, &q, &s, &xn, &pr, &ompr, &status, &bound);
  if (status) {
    show_error("cdfbin2", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return s;
}

double cdfbin3_wrap(double s, double p, double pr) {
  int which=3;
  double q=1.0-p, xn, ompr=1.0-pr, bound;
  int status;  

  F_FUNC(cdfbin,CDFBIN)(&which, &p, &q, &s, &xn, &pr, &ompr, &status, &bound); 
  if (status) {
    show_error("cdfbin3", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return xn;
}

extern void F_FUNC(cdfchi,CDFCHI)(int*,double*,double*,double*,double*,int*,double*);
double cdfchi3_wrap(double p, double x){
  int which=3;
  double q=1.0-p, df, bound;
  int status;  

  F_FUNC(cdfchi,CDFCHI)(&which, &p, &q, &x, &df, &status, &bound); 
  if (status) {
    show_error("cdfchi3", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return df;
}

extern void F_FUNC(cdfchn,CDFCHN)(int*,double*,double*,double*,double*,double*,int*,double*);
double cdfchn1_wrap(double x, double df, double nc) {
  int which=1;
  double q, p, bound;
  int status;  

  F_FUNC(cdfchn,CDFCHN)(&which, &p, &q, &x, &df, &nc, &status, &bound); 
  if (status) {
    show_error("cdfchn1", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return p;
}

double cdfchn2_wrap(double p, double df, double nc) {
  int which=2;
  double q=1.0-p, x, bound;
  int status;  

  F_FUNC(cdfchn,CDFCHN)(&which, &p, &q, &x, &df, &nc, &status, &bound); 
  if (status) {
    show_error("cdfchn2", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
  }
  return x;
}

double cdfchn3_wrap(double x, double p, double nc) {
  int which=3;
  double q=1.0-p, df, bound;
  int status;  

  F_FUNC(cdfchn,CDFCHN)(&which, &p, &q, &x, &df, &nc, &status, &bound); 
  if (status) {
    show_error("cdfchn3", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return df;
}

double cdfchn4_wrap(double x, double df, double p) {
  int which=4;
  double q=1.0-p, nc, bound;
  int status;  

  F_FUNC(cdfchn,CDFCHN)(&which, &p, &q, &x, &df, &nc, &status, &bound); 
  if (status) {
    show_error("cdfchn", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return nc;
}

extern void F_FUNC(cdff,CDFF)(int*,double*,double*,double*,double*,double*,int*,double*);
/*
double cdff1_wrap(double dfn, double dfd, double f) {
  int which=1;
  double q, p, bound;
  int status;

  F_FUNC(cdff,CDFF)(&which, &p, &q, &f, &dfn, &dfd, &status, &bound); 
  if (status) {
    show_error("cdff1", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
  }
  return p;
}

double cdff2_wrap(double dfn, double dfd, double p) {
  int which=2;
  double q=1.0-p, f, bound;
  int status;

  F_FUNC(cdff,CDFF)(&which, &p, &q, &f, &dfn, &dfd, &status, &bound); 
  if (status) {
    show_error("cdff2", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
  }
  return f;
}
*/

/* This seem to give some trouble.  No idea why... */
double cdff3_wrap(double p, double dfd, double f) {
  int which=3;
  double q=1.0-p, dfn, bound;
  int status;

  F_FUNC(cdff,CDFF)(&which, &p, &q, &f, &dfn, &dfd, &status, &bound); 
  if (status) {
    show_error("cdff3", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return dfn;
}

double cdff4_wrap(double dfn, double p, double f) {
  int which=4;
  double q=1.0-p, dfd, bound;
  int status;  

  F_FUNC(cdff,CDFF)(&which, &p, &q, &f, &dfn, &dfd, &status, &bound); 
  if (status) {
    show_error("cdff4", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return dfd;
}


extern void F_FUNC(cdffnc,CDFFNC)(int*,double*,double*,double*,double*,double*,double*,int*,double*);
double cdffnc1_wrap(double dfn, double dfd, double nc, double f) {
  int which=1;
  double q, p, bound;
  int status;

  F_FUNC(cdffnc,CDFFNC)(&which, &p, &q, &f, &dfn, &dfd, &nc, &status, &bound); 
  if (status) {
    show_error("cdffnc1", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
  }
  return p;
}

double cdffnc2_wrap(double dfn, double dfd, double nc, double p) {
  int which=2;
  double q=1.0-p, f, bound;
  int status;

  F_FUNC(cdffnc,CDFFNC)(&which, &p, &q, &f, &dfn, &dfd, &nc, &status, &bound); 
  if (status) {
    show_error("cdffnc2", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return f;
}


double cdffnc3_wrap(double p, double dfd, double nc, double f) {
  int which=3;
  double q=1.0-p, dfn, bound;
  int status;

  F_FUNC(cdffnc,CDFFNC)(&which, &p, &q, &f, &dfn, &dfd, &nc, &status, &bound); 
  if (status) {
    show_error("cdffnc3", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return dfn;
}
double cdffnc4_wrap(double dfn, double p, double nc, double f) {
  int which=4;
  double q=1.0-p, dfd, bound;
  int status;

  F_FUNC(cdffnc,CDFFNC)(&which, &p, &q, &f, &dfn, &dfd, &nc, &status, &bound); 
  if (status) {
    show_error("cdffnc4", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return dfd;
}

double cdffnc5_wrap(double dfn, double dfd, double p, double f) {
  int which=5;
  double q=1.0-p, nc, bound;
  int status;

  F_FUNC(cdffnc,CDFFNC)(&which, &p, &q, &f, &dfn, &dfd, &nc, &status, &bound); 
  if (status) {
    show_error("cdffnc5", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return nc;
}

/* scl == a in gdtr
   shp == b in gdtr
*/ 
extern void F_FUNC(cdfgam,CDFGAM)(int*,double*,double*,double*,double*,double*,int*,double*);
double cdfgam1_wrap(double scl, double shp, double x) {
  int which=1;
  double q, p, bound;
  int status;

  F_FUNC(cdfgam,CDFGAM)(&which, &p, &q, &x, &shp, &scl, &status, &bound); 
  if (status) {
    show_error("cdfgam1", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
  }
  return p;
}

double cdfgam2_wrap(double scl, double shp, double p) {
  int which=2;
  double q=1.0-p, x, bound;
  int status;

  F_FUNC(cdfgam,CDFGAM)(&which, &p, &q, &x, &shp, &scl,  &status, &bound); 
  if (status) {
    show_error("cdfgam2", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return x;
}

double cdfgam3_wrap(double scl, double p, double x) {
  int which=3;
  double q=1.0-p, shp, bound;
  int status;

  F_FUNC(cdfgam,CDFGAM)(&which, &p, &q, &x, &shp, &scl, &status, &bound); 
  if (status) {
    show_error("cdfgam3", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return shp;
}

double cdfgam4_wrap(double p, double shp, double x) {
  int which=4;
  double q=1.0-p, scl, bound;
  int status;

  F_FUNC(cdfgam,CDFGAM)(&which, &p, &q, &x, &shp, &scl, &status, &bound); 
  if (status) {
    show_error("cdfgam4", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return scl;
}

extern void F_FUNC(cdfnbn,CDFNBN)(int*,double*,double*,double*,double*,double*,double*,int*,double*);
double cdfnbn2_wrap(double p, double xn, double pr) {
  int which=2;
  double q=1.0-p, s, ompr=1.0-pr, bound;
  int status;  
  
  F_FUNC(cdfnbn,CDFNBN)(&which, &p, &q, &s, &xn, &pr, &ompr, &status, &bound);
  if (status) {
    show_error("cdfnbn2", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return s;
}

double cdfnbn3_wrap(double s, double p, double pr) {
  int which=3;
  double q=1.0-p, xn, ompr=1.0-pr, bound;
  int status;  

  F_FUNC(cdfnbn,CDFNBN)(&which, &p, &q, &s, &xn, &pr, &ompr, &status, &bound);
  if (status) {
    show_error("cdfnbn3", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return xn;
}

extern void F_FUNC(cdfnor,CDFNOR)(int*,double*,double*,double*,double*,double*,int*,double*);
double cdfnor3_wrap(double p, double std, double x) {
  int which=3;
  double q=1.0-p, mn, bound;
  int status;  

  F_FUNC(cdfnor,CDFNOR)(&which, &p, &q, &x, &mn, &std, &status, &bound); 
  if (status) {
    show_error("cdfnor3", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return mn;
}

double cdfnor4_wrap(double mn, double p, double x) {
  int which=4;
  double q=1.0-p, std, bound;
  int status;  

  F_FUNC(cdfnor,CDFNOR)(&which, &p, &q, &x, &mn, &std, &status, &bound); 
  if (status) {
    show_error("cdfnor4", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return std;
}

extern void F_FUNC(cdfpoi,CDFPOI)(int*,double*,double*,double*,double*,int*,double*);
double cdfpoi2_wrap(double p, double xlam){
  int which=2;
  double q=1.0-p, s, bound;
  int status;  

  F_FUNC(cdfpoi,CDFPOI)(&which, &p, &q, &s, &xlam, &status, &bound); 
  if (status) {
    show_error("cdfpoi2", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return s;
}

extern void F_FUNC(cdft,CDFT)(int*,double*,double*,double*,double*,int*,double*);
double cdft1_wrap(double df, double t){
  int which=1;
  double q, p, bound;
  int status;  

  F_FUNC(cdft,CDFT)(&which, &p, &q, &t, &df, &status, &bound); 
  if (status) {
    show_error("cdft1", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
  }
  return p;
}

double cdft2_wrap(double df, double p){
  int which=2;
  double q=1.0-p, t, bound;
  int status;  

  F_FUNC(cdft,CDFT)(&which, &p, &q, &t, &df, &status, &bound); 
  if (status) {
    show_error("cdft2", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return t;
}

double cdft3_wrap(double p, double t){
  int which=3;
  double q=1.0-p, df, bound;
  int status;  

  F_FUNC(cdft,CDFT)(&which, &p, &q, &t, &df, &status, &bound); 
  if (status) {
    show_error("cdft3", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return df;
}

extern void F_FUNC(cdftnc,CDFTNC)(int*,double*,double*,double*,double*,double*,int*,double*);
double cdftnc1_wrap(double df, double nc, double t) {
  int which=1;
  double q, p, bound;
  int status;  

  F_FUNC(cdftnc,CDFTNC)(&which, &p, &q, &t, &df, &nc, &status, &bound); 
  if (status) {
    show_error("cdftnc1", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return p;
}

double cdftnc2_wrap(double df, double nc, double p) {
  int which=2;
  double q=1.0-p, t, bound;
  int status;  

  F_FUNC(cdftnc,CDFTNC)(&which, &p, &q, &t, &df, &nc, &status, &bound); 
  if (status) {
    show_error("cdftnc2", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return t;
}

double cdftnc3_wrap(double p, double nc, double t) {
  int which=3;
  double q=1.0-p, df, bound;
  int status;  

  F_FUNC(cdftnc,CDFTNC)(&which, &p, &q, &t, &df, &nc, &status, &bound); 
  if (status) {
    show_error("cdftnc3", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;
  }
  return df;
}

double cdftnc4_wrap(double df, double p, double t) {
  int which=4;
  double q=1.0-p, nc, bound;
  int status;  

  F_FUNC(cdftnc,CDFTNC)(&which, &p, &q, &t, &df, &nc, &status, &bound); 
  if (status) {
    show_error("cdftnc4", status, bound);
    if ((status < 0) || (status==3) || (status==4)) return (NPY_NAN);
    if ((status == 1) || (status == 2)) return bound;

  }
  return nc;
}


