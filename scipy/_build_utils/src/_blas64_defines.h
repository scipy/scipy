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

/*
 * Define F_INT to match the type that the f2cmap selects for `integer`:
 * prefer `long` if it's 64-bit (LP64: Linux, macOS), fall back to
 * `long long` (LLP64: Windows).  This avoids pointer-type mismatch
 * warnings between f2py-generated local variables and F_INT* prototypes.
 *
 * Note that the below looks a little awkward, that's because of f2py limitations.
 * We can't simply add `int64_t` in the f2cmap, because that mapping mechanism
 * only accepts a given set of types and int64_t isn't part of that set.
 * The `abs` redefinition then follows from the F_INT one.
 *
 * Note that this code will go away once we can free ourselves of f2py
 * completely (which is planned).
 */
#include <limits.h>
#if LONG_MAX >= 0x7FFFFFFFFFFFFFFF
#define F_INT long
#else
#define F_INT long long
#endif

/*
 * f2py translates abs() from pyf expressions directly into C abs(),
 * which takes int. With ILP64, F_INT is long or long long, so we
 * need the matching absolute-value function to avoid truncation.
 */
#undef abs
#if LONG_MAX >= 0x7FFFFFFFFFFFFFFF
#define abs labs
#else
#define abs llabs
#endif

