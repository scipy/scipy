/*
 * Portable versions of some C99 functions. In cmath they are only
 * available in C++11 and the npy_math versions either don't work well
 * (they get undef'd by cmath; see gh-5689) or aren't
 * available. Implementations based on npy_math.h and friends.
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


#if __STDC_VERSION__ >= 199901L
/* We have C99 */
#define sc_fma(x, y, z) fma(x, y, z)
#else
/* No C99, define fma as x*y + z and maybe the compiler picks it up */
#define sc_fma(x, y, z) ((x)*(y) + (z))
#endif


#endif /* _c99compat.h */
