#ifndef BOOST_SPECIAL_FUNCTIONS_H
#define BOOST_SPECIAL_FUNCTIONS_H

#include <cmath>
#include <stdexcept>
#include "sf_error.h"

#include "boost/math/special_functions/erf.hpp"
#include "boost/math/special_functions/powm1.hpp"


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
        z = INFINITY;
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

#endif
