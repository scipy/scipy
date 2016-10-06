#ifndef NONFINITE_H
#define NONFINITE_H

/*
 * Portable isnan/isinf; in cmath only available in C++11 and npy_math
 * versions don't work well (they get undef'd by cmath, see gh-5689).
 * Implementation based on npy_math.h.
 */
#define sc_isnan(x) ((x) != (x))
#ifdef _MSC_VER
    #define sc_isfinite(x) _finite((x))
#else
    #define sc_isfinite(x) !sc_isnan((x) + (-x))
#endif
#define sc_isinf(x) (!sc_isfinite(x) && !sc_isnan(x))

#endif /* _nonfinite.h */
