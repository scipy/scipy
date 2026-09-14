#ifndef __SCIPY_BLAS_DEFINES__
#define __SCIPY_BLAS_DEFINES__


#include <stdint.h>

#ifdef __APPLE__
#include <AvailabilityMacros.h>
#endif

/* Allow the use in C++ code.  */
#ifdef __cplusplus
extern "C"
{
#endif


#ifdef ACCELERATE_NEW_LAPACK
    #if __MAC_OS_X_VERSION_MAX_ALLOWED < 130300
        #error "Accelerate support is only available with macOS 13.3 SDK or later"
    #else
        /* New Accelerate suffix is always $NEWLAPACK or $NEWLAPACK$ILP64 (no underscore appended) */
        #ifdef HAVE_BLAS_ILP64
            #define BLAS_SYMBOL_SUFFIX $NEWLAPACK$ILP64
        #else
            #define BLAS_SYMBOL_SUFFIX $NEWLAPACK
        #endif
    #endif
#endif

#ifdef NO_APPEND_FORTRAN
#define BLAS_FORTRAN_SUFFIX
#else
#define BLAS_FORTRAN_SUFFIX _
#endif

/* Accelerate doesn't use an underscore as suffix, so fix that up here */
#ifdef ACCELERATE_NEW_LAPACK
#undef BLAS_FORTRAN_SUFFIX
#define BLAS_FORTRAN_SUFFIX
#endif

#ifndef BLAS_SYMBOL_PREFIX
#define BLAS_SYMBOL_PREFIX
#endif

#ifndef BLAS_SYMBOL_SUFFIX
#define BLAS_SYMBOL_SUFFIX
#endif

#define BLAS_FUNC_CONCAT(name,prefix,suffix,suffix2) prefix ## name ## suffix ## suffix2
#define BLAS_FUNC_EXPAND(name,prefix,suffix,suffix2) BLAS_FUNC_CONCAT(name,prefix,suffix,suffix2)

/*
 * Use either the OpenBLAS scheme with the `64_` suffix behind the Fortran
 * compiler symbol mangling, or the MKL scheme (and upcoming
 * reference-lapack#666) which does it the other way around and uses `_64`.
 */
#ifdef OPENBLAS_ILP64_NAMING_SCHEME
#define BLAS_FUNC(name) BLAS_FUNC_EXPAND(name,BLAS_SYMBOL_PREFIX,BLAS_FORTRAN_SUFFIX,BLAS_SYMBOL_SUFFIX)
#else
#define BLAS_FUNC(name) BLAS_FUNC_EXPAND(name,BLAS_SYMBOL_PREFIX,BLAS_SYMBOL_SUFFIX,BLAS_FORTRAN_SUFFIX)
#endif

/*
 * Note that CBLAS doesn't include Fortran compiler symbol mangling, so ends up
 * being the same in both schemes
 */
#define CBLAS_FUNC(name) BLAS_FUNC_EXPAND(name,BLAS_SYMBOL_PREFIX,,BLAS_SYMBOL_SUFFIX)

#ifdef HAVE_BLAS_ILP64
#define CBLAS_INT int64_t
#define CBLAS_INT_MAX INT64_MAX
#else
#define CBLAS_INT int
#define CBLAS_INT_MAX INT_MAX
#endif


#ifdef __cplusplus
}
#endif

#endif  /* __SCIPY_BLAS_DEFINES__ */

