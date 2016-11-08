/*
 * Portable versions of some C99 functions. In cmath they are only
 * available in C++11 and the npy_math versions don't work well (they
 * get undef'd by cmath; see gh-5689). Implementations based on
 * npy_math.h and friends.
 */
#ifndef C99COMPAT_H
#define C99COMPAT_H


#define sc_isnan(x) ((x) != (x))
#ifdef _MSC_VER
    #define sc_isfinite(x) _finite((x))
#else
    #define sc_isfinite(x) !sc_isnan((x) + (-x))
#endif
#define sc_isinf(x) (!sc_isfinite(x) && !sc_isnan(x))


#endif /* _c99compat.h */
