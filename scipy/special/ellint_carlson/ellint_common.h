#ifndef _ELLINT_COMMON_H_INCLUDED
#define _ELLINT_COMMON_H_INCLUDED


#include "ellint_poly.h"
#include "ellint_compat.h"


#define ELLINT_MAXITER		(1000)


/* Sexatic root of real number */
#define ELLINT_SXRT(r)	( sqrt(cbrt( (r) )) )
/* Octic root of real number */
#define ELLINT_OCRT(r)	( sqrt(sqrt(sqrt( (r) ))) )

#define ELLINT_FMAX3(a, b, c)		( fmax(fmax( (a) , (b) ), (c) ) )
#define ELLINT_FMAX4(a, b, c, d)	( fmax(fmax( (a) , (b) ), \
					  fmax( (c), (d) )) )
#define ELLINT_FMIN4(a, b, c, d)	( fmin(fmin( (a) , (b) ), \
					  fmin( (c), (d) )) )
#define ELLINT_FMIN5(a, b, c, d, e)	\
    ( fmin(ELLINT_FMIN4( (a), (b), (c), (d) ), (e)) )


#define CZERO	( MKCR(0.0) )
#define CPINF	( MKCR(HUGE_VAL) )
#define CP_I	( MKCMPLX(0.0, 1.0) )
#define CN_I	( MKCMPLX(0.0, -1.0) )


#define ELLINT_FAIL_WITH(code)	\
do {	\
    *res = ELLINT_NAN;	\
    status = ( code );	\
    return status;	\
} while ( 0 )


#define ELLINT_RETURN_WITH(code, value)	\
do {	\
    *res = ( value );	\
    status = ( code );	\
    return status;	\
} while ( 0 )


#define HORRIBLE_STATUS(s)	( ( (s) == ELLINT_STATUS_NORESULT ) || \
                                  ( (s) == ELLINT_STATUS_BAD_ARGS ) || \
				  ( (s) == ELLINT_STATUS_BAD_RERR ) || \
				  ( (s) == ELLINT_STATUS_OTHER ) )
#define TROUBLESOME_STATUS(s)	( ( (s) == ELLINT_STATUS_SINGULAR ) || \
                                  ( (s) == ELLINT_STATUS_UNDERFLOW ) || \
				  ( (s) == ELLINT_STATUS_OVERFLOW ) || \
				  ( (s) == ELLINT_STATUS_NITER ) || \
				  ( (s) == ELLINT_STATUS_PRECLOSS ) )


#define ELLINT_SWAP(x, y)	\
do {				\
    EllInt_Num_t _swap_tmp;	\
    _swap_tmp = ( x );		\
    ( x ) = ( y );		\
    ( y ) = _swap_tmp;	\
} while ( 0 )


static const double ELLINT_RDJ_C1[4] = {0, -875160, 417690, -255255};
static const double ELLINT_RDJ_C2[3] = {0, 680680, 306306};
static const double ELLINT_RDJ_C3[3] = {0, -706860, 675675};
static const double ELLINT_RDJ_C4[2] = {-556920, 612612};
static const double ELLINT_RDJ_C5[2] = {471240, -540540};
#define ELLINT_RDJ_DENOM	(4084080.0)


#include "ellint_argcheck.h"


#endif  /* _ELLINT_COMMON_H_INCLUDED */
