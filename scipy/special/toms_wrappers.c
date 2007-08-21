/* This file is a collection (more can be added) of wrappers around some
 *  ToMS Fortran algorithm, so that they can be called from
 *  cephesmodule.so
 */

#include "toms_wrappers.h"
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
extern void F_FUNC(wofz,WOFZ)(double*,double*,double*,double*,int*);

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

