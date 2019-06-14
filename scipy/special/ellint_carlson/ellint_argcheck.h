#ifndef _ELLINT_ARGCHECK_H_INCLUDED
#define _ELLINT_ARGCHECK_H_INCLUDED

#include <stdbool.h>
#include <float.h>


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


/* The complex pair (a, b) are both non-zero. */
#define C_PAIR_NZ(a, b)		( !( too_small(fabs( (a) )) || \
                                     too_small(fabs( (b) )) ) )
/* At most one in the complex triplet (a, b, c) is zero. */
#define C_ATMOST_1Z(a, b, c)	( C_PAIR_NZ( (a), (b) ) || \
                                  C_PAIR_NZ( (b), (c) ) || \
				  C_PAIR_NZ( (c), (a) ) )

/* Three non-negative reals, at most only the first can be zero. To be used
 * with a, b, c permuted, or with a, b, c already in order. */
#define ELLINT_3NNEG_ATMOST1Z(a, b, c)	\
( ( (a) >= 0.0 ) && ( (b) > 0.0 ) && ( (c) > 0.0 ) )

/* The first number (a), given by real and imaginary parts, is real and
 * non-negative, while the other two are non-zero and complex conjugates. */
#define ELLINT_1R_2CONJ(ar, ai, b, c)		\
( ( (ar) >= 0.0 ) && ( too_small( (ai) ) ) &&	\
  ( C_PAIR_NZ( (b), (c) ) ) &&			\
  ( too_small( fabs(conj( (b) ) - (c) ) ) ) )


#define ELLINT_SWAP(x, y)	\
do {				\
    EllInt_Num_t _swap_tmp;	\
    _swap_tmp = ( x );		\
    ( x ) = ( y );		\
    ( y ) = _swap_tmp;	\
} while ( 0 )


#endif  /* _ELLINT_ARGCHECK_H_INCLUDED */
