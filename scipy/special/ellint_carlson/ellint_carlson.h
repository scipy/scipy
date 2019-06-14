#ifndef _ELLINT_CARLSON_H_INCLUDED
#define _ELLINT_CARLSON_H_INCLUDED
#include <complex.h>


#define ELLINT_STATUS_SUCCESS	(0)
#define ELLINT_STATUS_NITER	(4)
#define ELLINT_STATUS_BAD_ARGS	(7)
#define ELLINT_STATUS_BAD_RERR	(8)
#define ELLINT_STATUS_OTHER	(9)


#if ( __STDC_VERSION__ < 199901L )
#define restrict
#endif


extern int fellint_RF(double x, double y, double z, double rerr,
                      double * restrict res);

extern int cellint_RF(double complex x, double complex y, double complex z,
		      double rerr, double complex * restrict res);

extern int fellint_RD(double x, double y, double z, double rerr,
                      double * restrict res);

extern int cellint_RD(double complex x, double complex y, double complex z,
		      double rerr, double complex * restrict res);

extern int fellint_RJ(double x, double y, double z, double p,
                      double rerr, double * restrict res);

extern int cellint_RJ(double complex x, double complex y,
                      double complex z, double complex p,
                      double rerr, double complex * restrict res);

extern int fellint_RC(double x, double y, double rerr, double * restrict res);

extern int cellint_RC(double complex x, double complex y,
                      double rerr, double complex * restrict res);

extern int fellint_RG(double x, double y, double z, double rerr,
                      double * restrict res);

extern int cellint_RG(double complex x, double complex y, double complex z,
                      double rerr, double complex * restrict res);


#endif  /* _ELLINT_CARLSON_H_INCLUDED */
