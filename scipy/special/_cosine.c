/*
 *  The functions
 *
 *    cosine_cdf
 *    cosine_invcdf
 *
 *  defined here are the kernels for the ufuncs
 *
 *    _cosine_cdf
 *    _cosine_invcdf
 *
 *  defined in scipy.special._ufuncs.
 *
 *  The ufuncs are used by the class scipy.stats.cosine_gen.
 */

#include <math.h>
#include "cephes/polevl.h"

// M_PI64 is the 64 bit floating point representation of π, e.g.
//   >>> math.pi.hex()
//   '0x1.921fb54442d18p+1'
// It is used in the function cosine_cdf_pade_approx_at_neg_pi,
// which depends on this value being the 64 bit representation.
// Do not replace this with M_PI from math.h or NPY_PI from the
// numpy header files.
#define M_PI64 0x1.921fb54442d18p+1

//
// p and q (below) are the coefficients in the numerator and denominator
// polynomials (resp.) of the Pade approximation of
//     f(x) = (pi + x + sin(x))/(2*pi)
// at x=-pi.  The coefficients are ordered from lowest degree to highest.
// These values are used in the function cosine_cdf_pade_approx_at_neg_pi(x).
//
// These coefficients can be derived by using mpmath as follows:
//
//    import mpmath
//
//    def f(x):
//        return (mpmath.pi + x + mpmath.sin(x)) / (2*mpmath.pi)
//
//    # Note: 40 digits might be overkill; a few more digits than the default
//    # might be sufficient.
//    mpmath.mp.dps = 40
//    ts = mpmath.taylor(f, -mpmath.pi, 20)
//    p, q = mpmath.pade(ts, 9, 10)
//
// (A python script with that code is in special/_precompute/cosine_cdf.py.)
//
// The following are the values after converting to 64 bit floating point:
// p = [0.0,
//      0.0,
//      0.0,
//      0.026525823848649224,
//      0.0,
//      -0.0007883197097740538,
//      0.0,
//      1.0235408442872927e-05,
//      0.0,
//      -3.8360369451359084e-08]
// q = [1.0,
//      0.0,
//      0.020281047093125535,
//      0.0,
//      0.00020944197182753272,
//      0.0,
//      1.4162345851873058e-06,
//      0.0,
//      6.498171564823105e-09,
//      0.0,
//      1.6955280904096042e-11]
//


//
// Compute the CDF of the standard cosine distribution for x close to but
// not less than -π.  A Pade approximant is used to avoid the loss of
// precision that occurs in the formula 1/2 + (x + sin(x))/(2*pi) when
// x is near -π.
//
static
double cosine_cdf_pade_approx_at_neg_pi(double x)
{
    double h, h2, h3;
    double numer, denom;
    double numer_coeffs[] = {-3.8360369451359084e-08,
                             1.0235408442872927e-05,
                             -0.0007883197097740538,
                             0.026525823848649224};
    double denom_coeffs[] = {1.6955280904096042e-11,
                             6.498171564823105e-09,
                             1.4162345851873058e-06,
                             0.00020944197182753272,
                             0.020281047093125535,
                             1.0};

    // M_PI64 is not exactly π.  In fact, float64(π - M_PI64) is
    // 1.2246467991473532e-16.  h is supposed to be x + π, so to compute
    // h accurately, we write the calculation as:
    h = (x + M_PI64) + 1.2246467991473532e-16;
    h2 = h*h;
    h3 = h2*h;
    numer = h3*polevl(h2, numer_coeffs,
                      sizeof(numer_coeffs)/sizeof(numer_coeffs[0]) - 1);
    denom = polevl(h2, denom_coeffs,
                   sizeof(denom_coeffs)/sizeof(denom_coeffs[0]) - 1);
    return numer / denom;
}


//
// cosine distribution cumulative distribution function (CDF).
//
double cosine_cdf(double x)
{
    if (x >= M_PI64) {
        return 1;
    }
    if (x < -M_PI64) {
        return 0;
    }
    if (x < -1.6) {
        return cosine_cdf_pade_approx_at_neg_pi(x);
    }
    return 0.5 + (x + sin(x))/(2*M_PI64);
}


// The CDF of the cosine distribution is
//     p = (pi + x + sin(x)) / (2*pi),
// We want the inverse of this.
//
// Move the factor 2*pi and the constant pi to express this as
//     pi*(2*p - 1) = x + sin(x)
// Then if f(x) = x + sin(x), and g(x) is the inverse of f, we have
//     x = g(pi*(2*p - 1)).

// The coefficients in the functions _p2 and _q2 are the coefficients in the
// Pade approximation at p=0.5 to the inverse of x + sin(x).
// The following steps were used to derive these:
// 1. Find the coefficients of the Taylor polynomial of x + sin(x) at x = 0.
//    A Taylor polynomial of order 22 was used.
// 2. "Revert" the Taylor coefficients to find the Taylor polynomial of the
//    inverse.  The method for series reversion is described at
//        https://en.wikipedia.org/wiki/Bell_polynomials#Reversion_of_series
// 3. Convert the Taylor coefficients of the inverse to the (11, 10) Pade
//    approximant. The coefficients of the Pade approximant (converted to 64
//    bit floating point) in increasing order of degree are:
//      p = [0.0, 0.5,
//           0.0, -0.11602142940208726,
//           0.0, 0.009350454384541677,
//           0.0, -0.00030539712907115167,
//           0.0, 3.4900934227012284e-06,
//           0.0, -6.8448463845552725e-09]
//      q = [1.0,
//           0.0, -0.25287619213750784,
//           0.0, 0.022927496105281435,
//           0.0, -0.0008916919927321117,
//           0.0, 1.3728570152788793e-05,
//           0.0, -5.579679571562129e-08]
//
// The nonzero values in p and q are used in the functions _p2 and _q2 below.
// The functions assume that the square of the variable is passed in, and
// _p2 does not include the outermost multiplicative factor.
// So to evaluate x = invcdf(p) for a given p, the following is used:
//        y = pi*(2*p - 1)
//        x = y * _p2(y**2) / _q2(y**2)

static
double _p2(double t)
{
    double coeffs[] = {-6.8448463845552725e-09,
                       3.4900934227012284e-06,
                       -0.00030539712907115167,
                       0.009350454384541677,
                       -0.11602142940208726,
                       0.5};
    double v;

    v = polevl(t, coeffs, sizeof(coeffs) / sizeof(coeffs[0]) - 1);
    return v;
}

static
double _q2(double t)
{
    double coeffs[] = {-5.579679571562129e-08,
                       1.3728570152788793e-05,
                       -0.0008916919927321117,
                       0.022927496105281435,
                       -0.25287619213750784,
                       1.0};
    double v;

    v = polevl(t, coeffs, sizeof(coeffs) / sizeof(coeffs[0]) - 1);
    return v;
}


//
// Part of the asymptotic expansion of the inverse function at p=0.
//
// See, for example, the wikipedia article "Kepler's equation"
// (https://en.wikipedia.org/wiki/Kepler%27s_equation).  In particular, see the
// series expansion for the inverse Kepler equation when the eccentricity e is 1.
//
static
double _poly_approx(double s)
{
    double s2;
    double p;
    double coeffs[] = {1.1911667949082915e-08,
                       1.683039183039183e-07,
                       43.0/17248000,
                       1.0/25200,
                       1.0/1400,
                       1.0/60,
                       1.0};
    //
    // p(s) = s + (1/60) * s**3 + (1/1400) * s**5 + (1/25200) * s**7 +
    //        (43/17248000) * s**9 + (1213/7207200000) * s**11 +
    //        (151439/12713500800000) * s**13 + ...
    //
    // Here we include terms up to s**13.
    //
    s2 = s*s;
    p = s*polevl(s2, coeffs, sizeof(coeffs)/sizeof(coeffs[0]) - 1);
    return p;
}


//
// cosine distribution inverse CDF (aka percent point function).
//
double cosine_invcdf(double p)
{
    double x;
    int sgn = 1;

    if ((p < 0) || (p > 1)) {
        return NAN;
    }
    if (p <= 1e-48) {
        return -M_PI64;
    }
    if (p == 1) {
        return M_PI64;
    }

    if (p > 0.5) {
        p = 1.0 - p;
        sgn = -1;
    }

    if (p < 0.0925) {
        x = _poly_approx(cbrt(12*M_PI64*p)) - M_PI64;
    }
    else {
        double y, y2;
        y = M_PI64*(2*p - 1);
        y2 = y*y;
        x = y * _p2(y2) / _q2(y2);
    }

    // For p < 0.0018, the asymptotic expansion at p=0 is sufficently
    // accurate that no more work is needed.  Similarly, for p > 0.42,
    // the Pade approximant is sufficiently accurate.  In between these
    // bounds, we refine the estimate with Halley's method.
    if ((0.0018 < p) && (p < 0.42)) {
        // Apply one iteration of Halley's method, with
        //    f(x)   = pi + x + sin(x) - y,
        //    f'(x)  = 1 + cos(x),
        //    f''(x) = -sin(x)
        // where y = 2*pi*p.
        double f0, f1, f2;
        f0 = M_PI64 + x + sin(x) - 2*M_PI64*p;
        f1 = 1 + cos(x);
        f2 = -sin(x);
        x = x - 2*f0*f1/(2*f1*f1 - f0*f2);
    }

    return sgn*x;
}
