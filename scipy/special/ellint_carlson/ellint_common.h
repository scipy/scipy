#ifndef _ELLINT_COMMON_H_INCLUDED
#define _ELLINT_COMMON_H_INCLUDED

#include <math.h>
#include <tgmath.h>

#if ( (!defined(NDEBUG)) || (__STDC_VERSION__ < 199901L) )
#define inline
#endif

#if ( __STDC_VERSION__ < 199901L )
#define restrict
#endif

#define ELLINT_TOO_SMALL	(1e-11)
#define ELLINT_MAXITER		(1000)

#define ELLINT_BAD_RERR(rerr, tol) ( ( (rerr) <= 0.0 ) || ( (rerr) >= (tol) ) )

/* Sexatic root of real number */
#define ELLINT_SXRT(r)	( sqrt(cbrt( (r) )) )
/* Octic root of real number */
#define ELLINT_OCRT(r)	( sqrt(sqrt(sqrt( (r) ))) )

#define ELLINT_FMAX3(a, b, c)		( fmax(fmax( (a) , (b) ), (c) ) )
#define ELLINT_FMAX4(a, b, c, d)	( fmax(fmax( (a) , (b) ), \
					  fmax( (c), (d) )) )
#define ELLINT_FAIL_WITH(code)	\
do {				\
    *res = ELLINT_NAN;		\
    status = (code);		\
    return status;		\
} while ( 0 )


#include "ellint_carlson.h"

#define ELLINT_RDJ_c02	(-0.2142857142857142857142857)	/* -3 / 14 */
#define ELLINT_RDJ_c03	(+0.1666666666666666666666667)	/*  1 / 6  */
#define ELLINT_RDJ_c22	(+0.1022727272727272727272727)	/*  9 / 88 */
#define ELLINT_RDJ_c04	(-0.1363636363636363636363636)	/* -3 / 22 */
#define ELLINT_RDJ_c23	(-0.1730769230769230769230769)	/* -9 / 52 */
#define ELLINT_RDJ_c05	(+0.1153846153846153846153846)	/*  3 / 26 */


#endif  /* _ELLINT_COMMON_H_INCLUDED */
