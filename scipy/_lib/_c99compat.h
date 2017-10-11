/*
    C99 compatibility for SciPy
    
    The Rules:
      - Every function definition must have compiler guard so that it is not defined twice.
      - This file should only affect compilers that do not support C99 natively.
      - All functions should be defined as "static" to limit linker conflicts.
*/
#ifndef SCIPY_C99_COMPAT
    #define SCIPY_C99_COMPAT
    #include <float.h>
    
    #if defined(_MSC_VER) && _MSC_VER <= 1600
            /* MSVC 2008 and MSVC 2010 */
            #ifndef isnan
                #define isnan(x) _isnan((x))
            #endif

            #ifndef signbit
                static int signbit(double x)
                {
                    return x > 0;
                }
            #endif

            #ifndef copysign
                static double copysign(double x, double y)
                {
                    if (x >= 0) {
                        if (y >= 0) {
                            return x;
                        }
                        else {
                            return -x;
                        }
                    }
                    else {
                        if (y >= 0) {
                            return -x;
                        }
                        else {
                            return x;
                        }
                    }
                }
            #endif

            #ifndef fmax
                static double fmax(double x, double y)
                {
                    /* z > nan for z != nan is required by C the standard */
                    int xnan = isnan(x), ynan = isnan(y);

                    if (xnan || ynan)
                    {
                        if (xnan && !ynan)
                            return y;
                        if (!xnan && ynan)
                            return x;
                        return x;
                    }

                    /* +0 > -0 is preferred by C the standard */
                    if (x == 0 && y == 0)
                    {
                        int xs = signbit(x), ys = signbit(y);
                        if (xs && !ys)
                            return y;
                        if (!xs && ys)
                            return x;
                        return x;
                    }

                    if (x > y) {
                        return x;
                    }
                    else {
                        return y;
                    }
                }
            #endif

            #ifndef fmin
                static double fmin(double x, double y)
                {
                    /* z > nan for z != nan is required by C the standard */
                    int xnan = isnan(x), ynan = isnan(y);

                    if (xnan || ynan)
                    {
                        if (xnan && !ynan)
                            return y;
                        if (!xnan && ynan)
                            return x;
                        return x;
                    }

                    /* +0 > -0 is preferred by C the standard */
                    if (x == 0 && y == 0)
                    {
                        int xs = signbit(x), ys = signbit(y);
                        if (xs && !ys)
                            return x;
                        if (!xs && ys)
                            return y;
                        return y;
                    }

                    if (x > y) {
                        return y;
                    }
                    else {
                        return x;
                    }
                }
            #endif
    #endif

    #if (__STDC_VERSION__ < 199901L)
        /* Hopefully fail in less cases */

        /* For compilers which aren't MSVC and haven't defined isnan */
        #ifndef isnan
            #define isnan(x) ((x) != (x))
        #endif

        #ifndef isfinite
            #ifdef _MSC_VER
                /* MSVC 2015 and newer still don't have everything */
                #define isfinite(x) _finite((x))
            #else
                #define isfinite(x) !isnan((x) + (-x))
            #endif
        #endif

        #ifndef isinf
            #define isinf(x) (!isfinite(x) && !isnan(x))
        #endif

        #ifndef fma
            #define fma(x, y, z) ((x)*(y) + (z))
        #endif
    #endif

    /*
     * portable isnan/isinf; in cmath only available in C++11 and npy_math
     * versions don't work well (they get undef'd by cmath, see gh-5689)
     * Implementation based on npy_math.h
     */
    #ifndef sc_isnan
        #define sc_isnan isnan
    #endif
    #ifndef sc_isinf
        #define sc_isinf isinf
    #endif
    #ifndef sc_isfinite
        #define sc_isfinite isfinite
    #endif
    #ifndef sc_fma
        #define sc_fma fma
    #endif
#endif
