/* This file is a collection (more can be added) of wrappers around some
 *  ToMS Fortran algorithm, so that they can be called from
 *  cephesmodule.so
 */

#include "toms_wrappers.h"

/* This must be linked with g77
 */
Py_complex cwofz_wrap( Py_complex z) {
  int errflag;
  Py_complex cy;

  wofz_(CADDR(z), CADDR(cy), &errflag);
  if (errflag==1) mtherr("wofz:",3); /* wofz returns a single flag both
                                        for real overflows and for domain
                                        errors -- internal overflows from too
                                        large abs(z)*/
  return cy;
}

