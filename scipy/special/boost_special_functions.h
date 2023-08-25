#ifndef BOOST_SPECIAL_FUNCTIONS_H
#define BOOST_SPECIAL_FUNCTIONS_H

#include <cmath>
#include <stdexcept>
#include "sf_error.h"

#include "boost/math/special_functions/erf.hpp"
#include "boost/math/special_functions/powm1.hpp"
#include "boost/math/special_functions/hypergeometric_1F1.hpp"
#include "boost/math/special_functions/hypergeometric_pFq.hpp"

template<typename Real>
static inline
Real erfinv_wrap(Real x)
{
    Real y;

    if (x == -1) {
        return -INFINITY;
    }
    if (x == 1) {
        return INFINITY;
    }

    try {
        y = boost::math::erf_inv(x);
    } catch (const std::domain_error& e) {
        sf_error("erfinv", SF_ERROR_DOMAIN, NULL);
        y = NAN;
    } catch (const std::overflow_error& e) {
        sf_error("erfinv", SF_ERROR_OVERFLOW, NULL);
        y = INFINITY;
    } catch (const std::underflow_error& e) {
        sf_error("erfinv", SF_ERROR_UNDERFLOW, NULL);
        y = 0;
    } catch (...) {
        sf_error("erfinv", SF_ERROR_OTHER, NULL);
        y = NAN;
    }
    return y;
}

float
erfinv_float(float x)
{
    return erfinv_wrap(x);
}

double
erfinv_double(double x)
{
    return erfinv_wrap(x);
}


template<typename Real>
static inline
Real powm1_wrap(Real x, Real y)
{
    Real z;

    // Handle edge cases here instead of relying on boost.  This gives
    // better control of how we call sf_error().
    if (y == 0 || x == 1) {
        // (anything)**0 is 1
        // 1**(anything) is 1
        // This includes the case 0**0, and 'anything' includes inf and nan.
        return 0;
    }
    if (x == 0) {
        if (y < 0) {
            sf_error("powm1", SF_ERROR_DOMAIN, NULL);
            return INFINITY;
        }
        else if (y > 0) {
            return -1;
        }
    }
    if (x < 0 && std::trunc(y) != y) {
        // To compute x**y with x < 0, y must be an integer.
        sf_error("powm1", SF_ERROR_DOMAIN, NULL);
        return NAN;
    }

    try {
        z = boost::math::powm1(x, y);
    } catch (const std::domain_error& e) {
        sf_error("powm1", SF_ERROR_DOMAIN, NULL);
        z = NAN;
    } catch (const std::overflow_error& e) {
        sf_error("powm1", SF_ERROR_OVERFLOW, NULL);
        
        // See: https://en.cppreference.com/w/cpp/numeric/math/pow
        if (x > 0) {
            if (y < 0) {
                z = 0;
            }
            else if (y == 0) {
                z = 1;
            }
            else {
                z = INFINITY;
            }
        }
        else if (x == 0) {
            z = INFINITY;
        }
        else {
            if (y < 0) {
                if (std::fmod(y, 2) == 0) {
                    z = 0;
                }
                else {
                    z = -0;
                }
            }
            else if (y == 0) {
                z = 1;
            }
            else {
                if (std::fmod(y, 2) == 0) {
                    z = INFINITY;
                }
                else {
                    z = -INFINITY;
                }
            }
        }
    } catch (const std::underflow_error& e) {
        sf_error("powm1", SF_ERROR_UNDERFLOW, NULL);
        z = 0;
    } catch (...) {
        sf_error("powm1", SF_ERROR_OTHER, NULL);
        z = NAN;
    }
    return z;
}

float
powm1_float(float x, float y)
{
    return powm1_wrap(x, y);
}

double
powm1_double(double x, double y)
{
    return powm1_wrap(x, y);
}


//
// This wrapper of hypergeometric_pFq is here because there are a couple
// edge cases where hypergeometric_1F1 in Boost version 1.80 and earlier
// has a either bug or an inconsistent behavior.  It turns out that
// hypergeometric_pFq does the right thing in those cases, so we'll use
// it until our copy of Boost is updated.
//
template<typename Real>
static inline
Real call_hypergeometric_pFq(Real a, Real b, Real x)
{
    Real y;

    try {
        y = boost::math::hypergeometric_pFq({a}, {b}, x);
    } catch (const std::domain_error& e) {
        // The string "hyp1f1" is correct--don't change it to something like
        // "hypqfq".  The string refers to the SciPy special function, not the
        // underlying Boost function that is being called.
        sf_error("hyp1f1", SF_ERROR_DOMAIN, NULL);
        y = INFINITY;
    } catch (const std::overflow_error& e) {
        sf_error("hyp1f1", SF_ERROR_OVERFLOW, NULL);
        y = INFINITY;
    } catch (const std::underflow_error& e) {
        sf_error("hyp1f1", SF_ERROR_UNDERFLOW, NULL);
        y = 0;
    } catch (...) {
        sf_error("hyp1f1", SF_ERROR_OTHER, NULL);
        y = NAN;
    }
    return y;
}

template<typename Real>
static inline
Real hyp1f1_wrap(Real a, Real b, Real x)
{
    Real y;

    if (isnan(a) || isnan(b) || isnan(x)) {
        return NAN;
    }
    if (b <= 0 && std::trunc(b) == b) {
        // b is a non-positive integer.
        // Note: The code here is designed to preserve the historical behavior
        // of hyp1f1 in this edge case.  Other software, such as Boost, mpmath and
        // Mathematica, use a different convention for some of the subcases.
        if (b != 0  && a == b) {
            // In this case, use the Boost function hypergeometric_pFq
            // instead of hypergeometric_1F1.  This avoids an inconsistency
            // in Boost 1.80 and earlier; for details, see
            // https://github.com/boostorg/math/issues/829.
            return call_hypergeometric_pFq(a, b, x);
        }
        if (!(a < 0 && std::trunc(a) == a && a >= b)) {
            return INFINITY;
        }
        // Fall through and let the Boost function handle the remaining
        // cases.
    }
    if (a < 0 && std::trunc(a) == a && b > 0 && b == x) {
        // Avoid a bug in hypergeometric_1F1 in Boost 1.80 and earlier
        // that occurs when `a` is a negative integer, `b` is positive
        // and `b` == `x`.  hypergeometric_1F1 incorrectly sets a
        // floating point exception flag in this case; see
        // https://github.com/boostorg/math/issues/833.
        return call_hypergeometric_pFq(a, b, x);
    }

    // Use Boost's hypergeometric_1F1 for the basic calculation.  It should
    // also handle correctly any other special cases not covered above.
    // Catch all exceptions and handle them using the `special` error
    // handling function.
    try {
        y = boost::math::hypergeometric_1F1(a, b, x);
    } catch (const std::domain_error& e) {
        sf_error("hyp1f1", SF_ERROR_DOMAIN, NULL);
        y = INFINITY;
    } catch (const std::overflow_error& e) {
        sf_error("hyp1f1", SF_ERROR_OVERFLOW, NULL);
        y = INFINITY;
    } catch (const std::underflow_error& e) {
        sf_error("hyp1f1", SF_ERROR_UNDERFLOW, NULL);
        y = 0;
    } catch (...) {
        sf_error("hyp1f1", SF_ERROR_OTHER, NULL);
        y = NAN;
    }
    return y;
}

double
hyp1f1_double(double a, double b, double x)
{
    return hyp1f1_wrap(a, b, x);
}

//
// NOTE: It would be easy to also provide hyp1f1_float(...), but with the
// current ufunc generation code, it would not be used.  The ufunc
// generation code will implement the float types for the ufunc by
// casting the floats to double and using the double implementation.
// This is because we also have a complex version that is implemented
// in a different file, and the ufunc generation code requires just one
// kernel function per header when there are multiple headers (see the
// comments in _generate_pyx.py).
//

#endif
