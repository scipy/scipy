/*
 * A common include for fblas_64 and flapack_64 f2py sources.
 *
 * f2py accounts for the Fortran name mangling (upppercase/lowercase, trailing underscore),
 * via its hardcoded F_FUNC define.
 *
 * For ILP64 variants, we need a more flexible naming scheme, to potentially include
 * the _64 or 64_ suffixes. This is what the `BLAS_FUNC` macro from
 * `scipy_blas_defines.h` does.
 *
 * We therefore inject the define into the f2py-generated sources.
 *
 * To be more specific:
 *
 * 32-bit LAPACK builds generate the `blas_lapack_wrapper.so` library
 * which contains trampolines from F_FUNC to BLAS_FUNC mangling
 * (i.e. create dgemm_ as a synonym for scipy_dgemm_ etc).
 * Then f2py uses F_FUNC mangled symbols, and which get picked up from the wrappers.
 *
 * With ILP64 builds this approach may lead to symbol clashes: at least the MKL's
 * _ilp64.so variants contain two kinds of symbols: both `dgemm_64_` and `dgemm_`.
 * We want to only use the former (is dgemm_ 32- or 64-bit?).
 *
 * Therefore, we instead force the BLAS_FUNC mangling at the f2py level:
 * inject a 'usercode' section into f2py-generated C files and replace
 * the F_FUNC define with BLAS_FUNC.
 * This relies on the fact that f2py includes the usercode section _after_
 * its F_FUNC define.
 *
 * As a net result, the f2py-generated C code of the sort
 *
 * extern void F_FUNC(srotg,SROTG)(float*,float*,float*,float* );
 *
 * will in fact be `BLAS_FUNC(srotg)(...)` --- which accounts for the ILP64 name mangling.
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

#include "scipy_blas_defines.h"
#define F_FUNC(f, F) BLAS_FUNC(f)

#include "_mkl_ilp64_fixes.h"

#define F_INT int64_t

/*
 * f2py translates abs() from pyf expressions directly into C abs(),
 * which takes int. With ILP64, F_INT is npy_int64 (long), so we
 * need labs() to avoid truncation warnings.
 */
#undef abs
#define abs labs
