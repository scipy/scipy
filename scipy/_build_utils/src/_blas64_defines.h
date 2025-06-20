/*
 * A common include for fblas_64 and flapack_64 f2py sources.
 *
 * f2py accounts for the Fortran name mangling (upppercase/lowercase, trailing underscore),
 * via its hardcoded F_FUNC define.
 *
 * For ILP64 variants, we need a more flexible naming scheme, to potentially include
 * the _64 or 64_ suffixes. This is what the `BLAS_FUNC` macro from `npy_cblas.h` does.
 *
 * We therefore inject the define into the f2py-generated sources. 
 */

#ifndef HAVE_BLAS_ILP64
#error("HAVE_BLAS_ILP64 not defined.")
#endif
#ifndef BLAS_SYMBOL_SUFFIX
#error("BLAS_SYMBOL_SUFFIX  not defined")
#endif

#ifdef F_FUNC
#undef F_FUNC
#endif

#include "npy_cblas.h"
#define F_FUNC(f, F) BLAS_FUNC(f)

#ifdef FIX_MKL_2025_ILP64_MISSING_SYMBOL
#define cspr_64_ cspr_64
#endif

#define F_INT npy_int64

