/*
 * Complex math wrappers for use from Cython.
 *
 * These replace the npy_cabs, npy_csqrt, etc. that used to come from
 * NumPy's npymath library. We implement them here using platform
 * <complex.h> to avoid linking against npymath (which drags in long double
 * complex symbols that break MinGW builds).
 *
 * Both GCC/Clang and MSVC (VS2019+) provide the underlying C99 complex
 * functions; only the complex type names differ.
 */

#ifndef SCIPY_SPECIAL_COMPLEXSTUFF_H
#define SCIPY_SPECIAL_COMPLEXSTUFF_H

#include <numpy/npy_common.h>
#include "scipy_complex_support.h"

#if defined(_MSC_VER)
    typedef _Dcomplex _scipy_dz;
#else
    typedef double complex _scipy_dz;
#endif

#define _SCIPY_TO_DZ(z) CMPLX(((double *)&(z))[0], ((double *)&(z))[1])

static inline npy_cdouble _scipy_from_dz(_scipy_dz z) {
    npy_cdouble r;
    ((double *)&r)[0] = creal(z);
    ((double *)&r)[1] = cimag(z);
    return r;
}

static inline double scipy_cabs(npy_cdouble z) {
    return cabs(_SCIPY_TO_DZ(z));
}

static inline double scipy_carg(npy_cdouble z) {
    return carg(_SCIPY_TO_DZ(z));
}

static inline npy_cdouble scipy_clog(npy_cdouble z) {
    return _scipy_from_dz(clog(_SCIPY_TO_DZ(z)));
}

static inline npy_cdouble scipy_cexp(npy_cdouble z) {
    return _scipy_from_dz(cexp(_SCIPY_TO_DZ(z)));
}

static inline npy_cdouble scipy_csin(npy_cdouble z) {
    return _scipy_from_dz(csin(_SCIPY_TO_DZ(z)));
}

static inline npy_cdouble scipy_ccos(npy_cdouble z) {
    return _scipy_from_dz(ccos(_SCIPY_TO_DZ(z)));
}

static inline npy_cdouble scipy_csqrt(npy_cdouble z) {
    return _scipy_from_dz(csqrt(_SCIPY_TO_DZ(z)));
}

static inline npy_cdouble scipy_cpow(npy_cdouble x, npy_cdouble y) {
    return _scipy_from_dz(cpow(_SCIPY_TO_DZ(x), _SCIPY_TO_DZ(y)));
}

#endif /* SCIPY_SPECIAL_COMPLEXSTUFF_H */
