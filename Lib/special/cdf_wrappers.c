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


int status_to_mtherr( int status) {
     /* Return mtherr equivalents for ierr values */
  
  if (nz != 0) return UNDERFLOW;

  switch (ierr) {
  case 1:
    return DOMAIN;
  case 2:
    return OVERFLOW;
  case 3:
    return PLOSS;
  case 4:
    return TLOSS;
  case 5:   /* Algorithm termination condition not met */
    return TLOSS;    
  }
  return -1;
}

Py_complex cwofz_wrap( Py_complex z) {
  int errflag;
  Py_complex cy;

  F_FUNC(wofz,WOFZ)(CADDR(z), CADDR(cy), &errflag);
  if (errflag==1) mtherr("wofz:",3); /* wofz returns a single flag both
                                        for real overflows and for domain
                                        errors -- internal overflows from too
                                        large abs(z)*/
  return cy;
}


