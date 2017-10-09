/*
    C99 compatibility for SciPy
    
    The Rules:
      - Every function definition must have compiler guard so that it is not defined twice.
      - This file should only affect compilers that do not support C99 natively.
      - All functions should be defined as "static" to limit linker conflicts.
*/
#if !defined(SCIPY_C99_COMPAT) && defined(_MSC_VER) && _MSC_VER <= 1600
#define SCIPY_C99_COMPAT
    #include <float.h>

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
