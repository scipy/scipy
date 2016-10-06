/*
 * Portable versions of some C99 functions. In cmath they are only
 * available in C++11 and the npy_math versions don't work well (they
 * get undef'd by cmath; see gh-5689). Implementations based on
 * npy_math.h and friends.
 */
#ifndef C99COMPAT_H
#define C99COMPAT_H

#include <Python.h>
#include <numpy/numpyconfig.h>


int
sc_signbit(double x)
{
    union
    {
        double d;
        short s[4];
        int i[2];
    } u;

    u.d = x;

#if NPY_SIZEOF_INT == 4

#ifdef WORDS_BIGENDIAN /* defined in pyconfig.h */
    return u.i[0] < 0;
#else
    return u.i[1] < 0;
#endif

#else  /* NPY_SIZEOF_INT != 4 */

#ifdef WORDS_BIGENDIAN
    return u.s[0] < 0;
#else
    return u.s[3] < 0;
#endif

#endif  /* NPY_SIZEOF_INT */
}


#define sc_isnan(x) ((x) != (x))
#ifdef _MSC_VER
    #define sc_isfinite(x) _finite((x))
#else
    #define sc_isfinite(x) !sc_isnan((x) + (-x))
#endif
#define sc_isinf(x) (!sc_isfinite(x) && !sc_isnan(x))


#endif /* _c99compat.h */
