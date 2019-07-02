#ifndef _ELLINT_ARGCHECK_H_INCLUDED
#define _ELLINT_ARGCHECK_H_INCLUDED

#include <float.h>
#include "ellint_compat.h"


#define ELLINT_BAD_RERR(rerr, tol) ( ( (rerr) <= 0.0 ) || ( (rerr) >= (tol) ) )


static inline bool too_small(double t)
{
    int f;
    f = fpclassify(t);
    if ( (f == FP_ZERO) || (f == FP_SUBNORMAL) )
    {
	return true;
    }
    return false;
}


#define CEQZERO(z)	( (fpclassify(CREAL( (z) )) == FP_ZERO) && \
                          (fpclassify(CIMAG( (z) )) == FP_ZERO) )

#define Z_INFTY(z)	( fpclassify(FABS( (z) )) == FP_INFINITE )


/* principal value of complex phase is NOT +/- pi and not indetermined value */
static inline bool ph_is_not_pm_pi(EllInt_Num_t z)
{
    double r;
    int fp;
    r = CREAL(z);
    if ( ( fp = fpclassify(r) ) == FP_NAN )
    {
	return false;
    } else if ( fp == FP_INFINITE ) {
	if ( r < 0.0 )
	{
	    return false;
	}
	/* "bare" or true infinite point, the "north pole" of the Riemann
	 * sphere, without directional information about the approach to the
	 * point at infinity */
	if ( fpclassify(CIMAG(z)) == FP_INFINITE )
	{
	    return false;
	}
    } else if ( ( r < 0.0 ) && ( ( (fp = fpclassify(CIMAG(z))) == FP_ZERO ) ||
	                         ( fp == FP_NAN ) ) ) {
	return false;
    }
    return true;
}


#define PH_IS_PMPI_Z(z)		( !(ph_is_not_pm_pi( (z) )) )


/* The complex pair (a, b) are both non-zero. */
#define C_PAIR_NZ(a, b)		( !( too_small(FABS( (a) )) || \
                                     too_small(FABS( (b) )) ) )
/* The complex pair (a, b) are both non-zero and have phase less than pi. */
#define C_PAIR_NZ_NN(a, b)	( C_PAIR_NZ( (a), (b) ) && \
                                  ( !( PH_IS_PMPI_Z( (a) ) || \
				       PH_IS_PMPI_Z( (b) ) ) ) )
/* At most one in the complex triplet (a, b, c) is zero. */
#define C_ATMOST_1Z(a, b, c)	( C_PAIR_NZ( (a), (b) ) || \
                                  C_PAIR_NZ( (b), (c) ) || \
				  C_PAIR_NZ( (c), (a) ) )
/* At most one in the complex triplet (a, b, c) is zero and those non-zero ones
 * has phase less than pi. */
#define C_ATMOST_1Z_NN(a, b, c)	( C_PAIR_NZ_NN( (a), (b) ) || \
                                  C_PAIR_NZ_NN( (b), (c) ) || \
				  C_PAIR_NZ_NN( (c), (a) ) )

/* Three non-negative reals, at most only the first can be zero. To be used
 * with a, b, c permuted, or with a, b, c already in order. */
#define ELLINT_3NNEG_ATMOST1Z(a, b, c)	\
( ( (a) >= 0.0 ) && ( (b) > 0.0 ) && ( (c) > 0.0 ) )

/* The first number (a), given by real and imaginary parts, is real and
 * non-negative, while the other two are non-zero and complex conjugates and
 * have phases less than pi. */
#define ELLINT_1R_2CONJ(ar, ai, b, c)		\
( ( (ar) >= 0.0 ) && ( too_small( (ai) ) ) &&	\
  ( C_PAIR_NZ( (b), (c) ) ) &&			\
  ( !( PH_IS_PMPI_Z( (b) ) || PH_IS_PMPI_Z( (c) ) ) ) && \
  ( too_small( FABS(CONJ( (b) ) - (c) ) ) ) )


#endif  /* _ELLINT_ARGCHECK_H_INCLUDED */
