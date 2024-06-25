#ifndef ELLINT_CARLSON_WRAP_HH_INCLUDED
#define ELLINT_CARLSON_WRAP_HH_INCLUDED


#undef _BEGIN_EXTERN_C
#undef _END_EXTERN_C
#ifdef __cplusplus
#define _BEGIN_EXTERN_C	extern "C" {
#define _END_EXTERN_C	}
#else
#define _BEGIN_EXTERN_C	/* nothing */
#define _END_EXTERN_C	/* nothing */
#endif

_BEGIN_EXTERN_C
#include <numpy/npy_math.h>
_END_EXTERN_C

#define ELLINT_NO_VALIDATE_RELATIVE_ERROR_BOUND
#include "ellint_carlson_cpp_lite/ellint_carlson.hh"


_BEGIN_EXTERN_C


extern double fellint_RC(double x, double y);
extern npy_cdouble cellint_RC(npy_cdouble x, npy_cdouble y);

extern double fellint_RD(double x, double y, double z);
extern npy_cdouble cellint_RD(npy_cdouble x, npy_cdouble y, npy_cdouble z);

extern double fellint_RF(double x, double y, double z);
extern npy_cdouble cellint_RF(npy_cdouble x, npy_cdouble y, npy_cdouble z);

extern double fellint_RG(double x, double y, double z);
extern npy_cdouble cellint_RG(npy_cdouble x, npy_cdouble y, npy_cdouble z);

extern double fellint_RJ(double x, double y, double z, double p);
extern npy_cdouble cellint_RJ(npy_cdouble x, npy_cdouble y, npy_cdouble z, npy_cdouble p);


_END_EXTERN_C


#undef _BEGIN_EXTERN_C
#undef _END_EXTERN_C


#endif /* ELLINT_CARLSON_WRAP_HH_INCLUDED */
