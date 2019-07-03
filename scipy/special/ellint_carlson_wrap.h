#ifndef _ELLINT_CARLSON_WRAP_H_INCLUDED
#define _ELLINT_CARLSON_WRAP_H_INCLUDED

#include "Python.h"
#include "sf_error.h"
#include "ellint_carlson/ellint_carlson.h"
#include <numpy/npy_math.h>


#define SF_RERR		(2.5e-14)

#ifdef __cplusplus
extern "C" {
#endif

#ifndef _MSC_VER
    #ifndef CMPLX
	#if ( defined(__GNUC__) &&	\
	      ( __GNUC__ * 10000 + __GNUC_MINOR__ * 100 >= 40700 ) )
	    #define ELLINT_TO_C(z)	__builtin_complex( (z).real, (z).imag )
	#else	/* fallback */
	    #define ELLINT_TO_C(z)	\
		( (double_complex)((double)( (z).real ) + (double)( (z).imag ) * I) )
	#endif
    #else
	#define ELLINT_TO_C(z)		( CMPLX( (z).real, (z).imag ) )
    #endif
#else
    #define ELLINT_TO_C(z)	( _Cbuild( (z).real, (z).imag ) )
#endif

#define ELLINT_MAKE_N(t, z)	\
do { (t).real = creal( (z) ); (t).imag = cimag( (z) ); } while ( 0 )

extern npy_double fellint_RC_w(npy_double x, npy_double y);
extern npy_cdouble cellint_RC_w(npy_cdouble x, npy_cdouble y);

extern npy_double fellint_RD_w(npy_double x, npy_double y, npy_double z);
extern npy_cdouble cellint_RD_w(npy_cdouble x, npy_cdouble y, npy_cdouble z);

extern npy_double fellint_RF_w(npy_double x, npy_double y, npy_double z);
extern npy_cdouble cellint_RF_w(npy_cdouble x, npy_cdouble y, npy_cdouble z);

extern npy_double fellint_RG_w(npy_double x, npy_double y, npy_double z);
extern npy_cdouble cellint_RG_w(npy_cdouble x, npy_cdouble y, npy_cdouble z);

extern npy_double fellint_RJ_w(npy_double x, npy_double y,
                               npy_double z, npy_double p);
extern npy_cdouble cellint_RJ_w(npy_cdouble x, npy_cdouble y,
                                npy_cdouble z, npy_cdouble p);

#ifdef __cplusplus
}
#endif
#endif /* _ELLINT_CARLSON_WRAP_H_INCLUDED */
