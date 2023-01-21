#include <math.h>

//
// We use exp1_wrap() from specfun_wrappers.c
//
extern double exp1_wrap(double x);

//
// Compute a factor of the result of the exponential integral E1.
// This is used in logexp1(x) for 1 < x <= 500.
//
// The function uses the continued fraction expansion given in equation 5.1.22
// of Abramowitz & Stegun, "Handbook of Mathematical Functions".
// For n=1, this is
//    E1(x) = exp(-x)*F(x)
// where F(x) is expressed as a continued fraction:
//    F(x) =                 1
//           ---------------------------------------
//                              1
//           x + ------------------------------------
//                                 1
//               1 + ---------------------------------
//                                    2
//                   x + ------------------------------
//                                       2
//                       1 + ---------------------------
//                                          3
//                           x + ------------------------
//                                             3
//                               1 + ---------------------
//                                                4
//                                   x + ------------------
//                                       1 +     [...]
//
static double
expint1_cont_frac_factor(double x)
{
    // The number of terms to use in the truncated continued fraction
    // depends on x.  Larger values of x require fewer terms.
    int m = 20 + (int) (80.0 / x);
    double t0 = 0.0;
    for (int k = m; k > 0; --k) {
        t0 = k/(1 + k/(x + t0));
    }
    return 1/(x + t0);
}

//
// Log of the exponential integral function E1 (for real x only).
//
double
logexp1(double x)
{
    if (x < 0) {
        return NAN;
    }
    if (x == 0) {
        return INFINITY;
    }
    if (x <= 1) {
        // For small x, the naive implementation is sufficiently accurate.
        return log(exp1_wrap(x));
    }
    if (x <= 500) {
        // For moderate x, use the continued fraction expansion.
        double t = expint1_cont_frac_factor(x);
        return -x + log(t);
    }
    // For large x, use the asymptotic expansion.  This is equation 5.1.51
    // from Abramowitz & Stegun, "Handbook of Mathematical Functions".
    double s = (-1 + (2 + (-6 + (24 - 120/x)/x)/x)/x)/x;
    return -x - log(x) + log1p(s);
}
