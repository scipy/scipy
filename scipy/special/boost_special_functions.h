#ifndef BOOST_SPECIAL_FUNCTIONS_H
#define BOOST_SPECIAL_FUNCTIONS_H

#include <cmath>
#include <stdexcept>
#include "sf_error.h"


#include "boost/math/special_functions/beta.hpp"
#include "boost/math/special_functions/erf.hpp"
#include "boost/math/special_functions/powm1.hpp"
#include "boost/math/special_functions/hypergeometric_1F1.hpp"
#include "boost/math/special_functions/hypergeometric_pFq.hpp"

#include "boost/math/distributions.hpp"
#include <boost/math/distributions/inverse_gaussian.hpp>

typedef boost::math::policies::policy<
    boost::math::policies::promote_float<false >,
    boost::math::policies::promote_double<false >,
    boost::math::policies::max_root_iterations<400 > > SpecialPolicy;

// Round up to achieve correct ppf(cdf) round-trips for discrete distributions
typedef boost::math::policies::policy<
    boost::math::policies::domain_error<boost::math::policies::ignore_error >,
    boost::math::policies::overflow_error<boost::math::policies::user_error >,
    boost::math::policies::evaluation_error<boost::math::policies::user_error >,
    boost::math::policies::promote_float<false >,
    boost::math::policies::promote_double<false >,
    boost::math::policies::discrete_quantile<
        boost::math::policies::integer_round_up > > StatsPolicy;

// Raise a RuntimeWarning making users aware that something went wrong during
// evaluation of the function, but return the best guess
template <class RealType>
RealType
boost::math::policies::user_evaluation_error(const char* function, const char* message, const RealType& val) {
    std::string msg("Error in function ");
    std::string haystack {function};
    const std::string needle {"%1%"};
    msg += haystack.replace(haystack.find(needle), needle.length(), typeid(RealType).name()) + ": ";
    // "message" may have %1%, but arguments don't always contain all
    // required information, so don't call boost::format for now
    msg += message;
    PyGILState_STATE save = PyGILState_Ensure();
    PyErr_WarnEx(PyExc_RuntimeWarning, msg.c_str(), 1);
    PyGILState_Release(save);
    return val;
}


template <class RealType>
RealType
boost::math::policies::user_overflow_error(const char* function, const char* message, const RealType& val) {
    std::string msg("Error in function ");
    std::string haystack {function};
    const std::string needle {"%1%"};
    msg += haystack.replace(haystack.find(needle), needle.length(), typeid(RealType).name()) + ": ";
    // From Boost docs: "overflow and underflow messages do not contain this %1% specifier
    //                   (since the value of value is immaterial in these cases)."
    if (message) {
        msg += message;
    }
    PyGILState_STATE save = PyGILState_Ensure();
    PyErr_SetString(PyExc_OverflowError, msg.c_str());
    PyGILState_Release(save);
    return 0;
}

template<typename Real>
static inline
Real ibeta_wrap(Real a, Real b, Real x)
{
    Real y;

    if (isnan(a) || isnan(b) || isnan(x)) {
        return NAN;
    }
    if ((a <= 0) || (b <= 0) || (x < 0) || (x > 1)) {
        sf_error("betainc", SF_ERROR_DOMAIN, NULL);
        return NAN;
    }
    try {
        y = boost::math::ibeta(a, b, x, SpecialPolicy());
    } catch (const std::domain_error& e) {
        sf_error("betainc", SF_ERROR_DOMAIN, NULL);
        y = NAN;
    } catch (const std::overflow_error& e) {
        sf_error("betainc", SF_ERROR_OVERFLOW, NULL);
        y = INFINITY;
    } catch (const std::underflow_error& e) {
        sf_error("betainc", SF_ERROR_UNDERFLOW, NULL);
        y = 0;
    } catch (...) {
        sf_error("betainc", SF_ERROR_OTHER, NULL);
        y = NAN;
    }
    return y;
}

float
ibeta_float(float a, float b, float x)
{
    return ibeta_wrap(a, b, x);
}

double
ibeta_double(double a, double b, double x)
{
    return ibeta_wrap(a, b, x);
}


template<typename Real>
static inline
Real ibetac_wrap(Real a, Real b, Real x)
{
    Real y;

    if (isnan(a) || isnan(b) || isnan(x)) {
        return NAN;
    }
    if ((a <= 0) || (b <= 0) || (x < 0) || (x > 1)) {
        sf_error("betaincc", SF_ERROR_DOMAIN, NULL);
        return NAN;
    }
    try {
        y = boost::math::ibetac(a, b, x);
    } catch (const std::domain_error& e) {
        sf_error("betaincc", SF_ERROR_DOMAIN, NULL);
        y = NAN;
    } catch (const std::overflow_error& e) {
        sf_error("betaincc", SF_ERROR_OVERFLOW, NULL);
        y = INFINITY;
    } catch (const std::underflow_error& e) {
        sf_error("betaincc", SF_ERROR_UNDERFLOW, NULL);
        y = 0;
    } catch (...) {
        sf_error("betaincc", SF_ERROR_OTHER, NULL);
        y = NAN;
    }
    return y;
}

float
ibetac_float(float a, float b, float x)
{
    return ibetac_wrap(a, b, x);
}

double
ibetac_double(double a, double b, double x)
{
    return ibetac_wrap(a, b, x);
}


template<typename Real, typename Policy>
static inline
Real ibeta_inv_wrap(Real a, Real b, Real p, const Policy& policy_)
{
    Real y;

    if (isnan(a) || isnan(b) || isnan(p)) {
        return NAN;
    }
    if ((a <= 0) || (b <= 0) || (p < 0) || (p > 1)) {
        sf_error("betaincinv", SF_ERROR_DOMAIN, NULL);
        return NAN;
    }
    try {
        y = boost::math::ibeta_inv(a, b, p, policy_);
    } catch (const std::domain_error& e) {
        sf_error("betaincinv", SF_ERROR_DOMAIN, NULL);
        y = NAN;
    } catch (const std::overflow_error& e) {
        sf_error("betaincinv", SF_ERROR_OVERFLOW, NULL);
        y = INFINITY;
    } catch (const std::underflow_error& e) {
        sf_error("betaincinv", SF_ERROR_UNDERFLOW, NULL);
        y = 0;
    } catch (...) {
        sf_error("betaincinv", SF_ERROR_OTHER, NULL);
        y = NAN;
    }
    return y;
}

float
ibeta_inv_float(float a, float b, float p)
{
    return ibeta_inv_wrap(a, b, p, SpecialPolicy());
}

double
ibeta_inv_double(double a, double b, double p)
{
    return ibeta_inv_wrap(a, b, p, SpecialPolicy());
}


template<typename Real>
static inline
Real ibetac_inv_wrap(Real a, Real b, Real p)
{
    Real y;

    if (isnan(a) || isnan(b) || isnan(p)) {
        return NAN;
    }
    if ((a <= 0) || (b <= 0) || (p < 0) || (p > 1)) {
        sf_error("betainccinv", SF_ERROR_DOMAIN, NULL);
        return NAN;
    }
    try {
        y = boost::math::ibetac_inv(a, b, p, SpecialPolicy());
    } catch (const std::domain_error& e) {
        sf_error("betainccinv", SF_ERROR_DOMAIN, NULL);
        y = NAN;
    } catch (const std::overflow_error& e) {
        sf_error("betainccinv", SF_ERROR_OVERFLOW, NULL);
        y = INFINITY;
    } catch (const std::underflow_error& e) {
        sf_error("betainccinv", SF_ERROR_UNDERFLOW, NULL);
        y = 0;
    } catch (...) {
        sf_error("betainccinv", SF_ERROR_OTHER, NULL);
        y = NAN;
    }
    return y;
}

float
ibetac_inv_float(float a, float b, float p)
{
    return ibetac_inv_wrap(a, b, p);
}

double
ibetac_inv_double(double a, double b, double p)
{
    return ibetac_inv_wrap(a, b, p);
}

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
        y = boost::math::erf_inv(x, SpecialPolicy());
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
        z = boost::math::powm1(x, y, SpecialPolicy());
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
        Real p_abs_error = 0;
        y = boost::math::hypergeometric_pFq({a}, {b}, x, &p_abs_error, SpecialPolicy());
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
        // Real p_abs_error = 0;
        y = boost::math::hypergeometric_1F1(a, b, x, SpecialPolicy());
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

// patch for boost::math::beta_distribution throwing exception for
// x = 1, beta < 1 as well as x = 0, alpha < 1
template<typename Real>
Real beta_pdf_wrap(const Real x, const Real a, const Real b)
{
    if (std::isfinite(x)) {
        if ((x >= 1) && (b < 1)) {
            // x>1 should really be 0, but rv_continuous will do that for us
            return INFINITY;
        }
        else if ((x <= 0) && (a < 1)) {
            return INFINITY;
        }
        return boost::math::pdf(
            boost::math::beta_distribution<Real, StatsPolicy>(a, b), x);
    }
    return NAN;
}

float
beta_pdf_float(float x, float a, float b)
{
    return beta_pdf_wrap(x, a, b);
}

double
beta_pdf_double(double x, double a, double b)
{
    return beta_pdf_wrap(x, a, b);
}

template<typename Real>
Real beta_ppf_wrap(const Real x, const Real a, const Real b)
{
    typedef boost::math::policies::policy<
        boost::math::policies::domain_error<boost::math::policies::ignore_error >,
        boost::math::policies::overflow_error<boost::math::policies::user_error >,
        boost::math::policies::evaluation_error<boost::math::policies::user_error >,
        boost::math::policies::promote_double<false > > BetaPolicyForStats;

    return ibeta_inv_wrap(a, b, x, BetaPolicyForStats());
}

float
beta_ppf_float(float x, float a, float b)
{
    return beta_ppf_wrap(x, a, b);
}

double
beta_ppf_double(double x, double a, double b)
{
    return beta_ppf_wrap(x, a, b);
}

template<typename Real>
Real
invgauss_ppf_wrap(const Real x, const Real mu, const Real s)
{
    return boost::math::quantile(
        boost::math::inverse_gaussian_distribution<Real, StatsPolicy>(mu, s), x);
}

float
invgauss_ppf_float(float x, float mu, float s)
{
    return invgauss_ppf_wrap(x, mu, s);
}

double
invgauss_ppf_double(double x, double mu, double s)
{
    return invgauss_ppf_wrap(x, mu, s);
}

template<typename Real>
Real
invgauss_isf_wrap(const Real x, const Real mu, const Real s)
{
    return boost::math::quantile(boost::math::complement(
        boost::math::inverse_gaussian_distribution<Real, StatsPolicy>(mu, s), x));
}

float
invgauss_isf_float(float x, float mu, float s)
{
    return invgauss_isf_wrap(x, mu, s);
}

double
invgauss_isf_double(double x, double mu, double s)
{
    return invgauss_isf_wrap(x, mu, s);
}

template<typename Real>
Real
cauchy_ppf_wrap(const Real p, const Real loc, const Real scale)
{
    return boost::math::quantile(
        boost::math::cauchy_distribution<Real, StatsPolicy>(loc, scale), p);
}

float
cauchy_ppf_float(float p, float loc, float scale)
{
    return cauchy_ppf_wrap(p, loc, scale);
}

double
cauchy_ppf_double(double p, double loc, double scale)
{
    return cauchy_ppf_wrap(p, loc, scale);
}

template<typename Real>
Real
cauchy_isf_wrap(const Real p, const Real loc, const Real scale)
{
    return boost::math::quantile(boost::math::complement(
        boost::math::cauchy_distribution<Real, StatsPolicy>(loc, scale), p));
}

float
cauchy_isf_float(float p, float loc, float scale)
{
    return cauchy_isf_wrap(p, loc, scale);
}

double
cauchy_isf_double(double p, double loc, double scale)
{
    return cauchy_isf_wrap(p, loc, scale);
}

template<typename Real>
Real
ncx2_pdf_wrap(const Real x, const Real k, const Real l)
{
    if (std::isfinite(x)) {
        return boost::math::pdf(
            boost::math::non_central_chi_squared_distribution<Real, StatsPolicy>(k, l), x);
    }
    return NAN; // inf or -inf returns NAN
}

float
ncx2_pdf_float(float x, float k, float l)
{
    return ncx2_pdf_wrap(x, k, l);
}

double
ncx2_pdf_double(double x, double k, double l)
{
    return ncx2_pdf_wrap(x, k, l);
}

template<typename Real>
Real
ncx2_cdf_wrap(const Real x, const Real k, const Real l)
{
    if (std::isfinite(x)) {
        return boost::math::cdf(
            boost::math::non_central_chi_squared_distribution<Real, StatsPolicy>(k, l), x);
    }
    // -inf => 0, inf => 1
    return 1 - std::signbit(x);
}

float
ncx2_cdf_float(float x, float k, float l)
{
    return ncx2_cdf_wrap(x, k, l);
}

double
ncx2_cdf_double(double x, double k, double l)
{
    return ncx2_cdf_wrap(x, k, l);
}

template<typename Real>
Real
ncx2_ppf_wrap(const Real x, const Real k, const Real l)
{
    return boost::math::quantile(
        boost::math::non_central_chi_squared_distribution<Real, StatsPolicy>(k, l), x);
}

float
ncx2_ppf_float(float x, float k, float l)
{
    return ncx2_ppf_wrap(x, k, l);
}

double
ncx2_ppf_double(double x, double k, double l)
{
    return ncx2_ppf_wrap(x, k, l);
}

template<typename Real>
Real
ncx2_sf_wrap(const Real x, const Real k, const Real l)
{
    return boost::math::cdf(boost::math::complement(
        boost::math::non_central_chi_squared_distribution<Real, StatsPolicy>(k, l), x));
}

float
ncx2_sf_float(float x, float k, float l)
{
    return ncx2_sf_wrap(x, k, l);
}

double
ncx2_sf_double(double x, double k, double l)
{
    return ncx2_sf_wrap(x, k, l);
}

template<typename Real>
Real
ncx2_isf_wrap(const Real x, const Real k, const Real l)
{
    return boost::math::quantile(boost::math::complement(
        boost::math::non_central_chi_squared_distribution<Real, StatsPolicy>(k, l), x));
}

float
ncx2_isf_float(float x, float k, float l)
{
    return ncx2_isf_wrap(x, k, l);
}

double
ncx2_isf_double(double x, double k, double l)
{
    return ncx2_isf_wrap(x, k, l);
}

template<typename Real>
Real
ncf_pdf_wrap(const Real x, const Real v1, const Real v2, const Real l)
{
    if (std::isfinite(x)) {
        return boost::math::pdf(
            boost::math::non_central_f_distribution<Real, StatsPolicy>(v1, v2, l), x);
    }
    return NAN; // inf or -inf returns NAN
}

float
ncf_pdf_float(float x, float v1, float v2, float l)
{
    return ncf_pdf_wrap(x, v1, v2, l);
}

double
ncf_pdf_double(double x, double v1, double v2, double l)
{
    return ncf_pdf_wrap(x, v1, v2, l);
}

template<typename Real>
Real
ncf_cdf_wrap(const Real x, const Real v1, const Real v2, const Real l)
{
    if (std::isfinite(x)) {
        return boost::math::cdf(
            boost::math::non_central_f_distribution<Real, StatsPolicy>(v1, v2, l), x);
    }
    // -inf => 0, inf => 1
    return 1.0 - std::signbit(x);
}

float
ncf_cdf_float(float x, float v1, float v2, float l)
{
    return ncf_cdf_wrap(x, v1, v2, l);
}

double
ncf_cdf_double(double x, double v1, double v2, double l)
{
    return ncf_cdf_wrap(x, v1, v2, l);
}

template<typename Real>
Real
ncf_ppf_wrap(const Real x, const Real v1, const Real v2, const Real l)
{
    return boost::math::quantile<Real, StatsPolicy>(
        boost::math::non_central_f_distribution<Real, StatsPolicy>(v1, v2, l), x);
}

float
ncf_ppf_float(float x, float v1, float v2, float l)
{
    return ncf_ppf_wrap(x, v1, v2, l);
}

double
ncf_ppf_double(double x, double v1, double v2, double l)
{
    return ncf_ppf_wrap(x, v1, v2, l);
}

template<typename Real>
Real
ncf_sf_wrap(const Real x, const Real v1, const Real v2, const Real l)
{
    return boost::math::cdf(boost::math::complement(
        boost::math::non_central_f_distribution<Real, StatsPolicy>(v1, v2, l), x));
}

float
ncf_sf_float(float x, float v1, float v2, float l)
{
    return ncf_sf_wrap(x, v1, v2, l);
}

double
ncf_sf_double(double x, double v1, double v2, double l)
{
    return ncf_sf_wrap(x, v1, v2, l);
}

template<typename Real>
Real
ncf_isf_wrap(const Real x, const Real v1, const Real v2, const Real l)
{
    return boost::math::quantile(boost::math::complement(
        boost::math::non_central_f_distribution<Real, StatsPolicy>(v1, v2, l), x));
}

float
ncf_isf_float(float x, float v1, float v2, float l)
{
    return ncf_isf_wrap(x, v1, v2, l);
}

double
ncf_isf_double(double x, double v1, double v2, double l)
{
    return ncf_isf_wrap(x, v1, v2, l);
}

#define RETURN_NAN(v, c) if( v <= c ) { \
        return NAN; \
    } \

float
ncf_mean_float(float v1, float v2, float l)
{
    RETURN_NAN(v2, 2.0);
    return boost::math::mean(boost::math::non_central_f_distribution<float, StatsPolicy>(v1, v2, l));
}

double
ncf_mean_double(double v1, double v2, double l)
{
    RETURN_NAN(v2, 2.0);
    return boost::math::mean(boost::math::non_central_f_distribution<double, StatsPolicy>(v1, v2, l));
}

float
ncf_variance_float(float v1, float v2, float l)
{
    RETURN_NAN(v2, 4.0);
    return boost::math::variance(boost::math::non_central_f_distribution<float, StatsPolicy>(v1, v2, l));
}

double
ncf_variance_double(double v1, double v2, double l)
{
    RETURN_NAN(v2, 4.0);
    return boost::math::variance(boost::math::non_central_f_distribution<double, StatsPolicy>(v1, v2, l));
}

float
ncf_skewness_float(float v1, float v2, float l)
{
    RETURN_NAN(v2, 6.0);
    return boost::math::skewness(boost::math::non_central_f_distribution<float, StatsPolicy>(v1, v2, l));
}

double
ncf_skewness_double(double v1, double v2, double l)
{
    RETURN_NAN(v2, 6.0);
    return boost::math::skewness(boost::math::non_central_f_distribution<double, StatsPolicy>(v1, v2, l));
}

float
ncf_kurtosis_excess_float(float v1, float v2, float l)
{
    RETURN_NAN(v2, 8.0);
    return boost::math::kurtosis_excess(boost::math::non_central_f_distribution<float, StatsPolicy>(v1, v2, l));
}

double
ncf_kurtosis_excess_double(double v1, double v2, double l)
{
    RETURN_NAN(v2, 8.0);
    return boost::math::kurtosis_excess(boost::math::non_central_f_distribution<double, StatsPolicy>(v1, v2, l));
}

template<typename Real>
Real
nct_cdf_wrap(const Real x, const Real v, const Real l)
{
    if (std::isfinite(x)) {
        return boost::math::cdf(
            boost::math::non_central_t_distribution<Real, StatsPolicy>(v, l), x);
    }
    // -inf => 0, inf => 1
    return 1.0 - std::signbit(x);
}

float
nct_cdf_float(float x, float v, float l)
{
    return nct_cdf_wrap(x, v, l);
}

double
nct_cdf_double(double x, double v, double l)
{
    return nct_cdf_wrap(x, v, l);
}

template<typename Real>
Real
nct_pdf_wrap(const Real x, const Real v, const Real l)
{
    if (std::isfinite(x)) {
        return boost::math::pdf(
            boost::math::non_central_t_distribution<Real, StatsPolicy>(v, l), x);
    }
    return NAN;
}

float
nct_pdf_float(float x, float v, float l)
{
    return nct_pdf_wrap(x, v, l);
}

double
nct_pdf_double(double x, double v, double l)
{
    return nct_pdf_wrap(x, v, l);
}

template<typename Real>
Real
nct_ppf_wrap(const Real x, const Real v, const Real l)
{
    return boost::math::quantile(
        boost::math::non_central_t_distribution<Real, StatsPolicy>(v, l), x);
}

float
nct_ppf_float(float x, float v, float l)
{
    return nct_ppf_wrap(x, v, l);
}

double
nct_ppf_double(double x, double v, double l)
{
    return nct_ppf_wrap(x, v, l);
}

template<typename Real>
Real
nct_sf_wrap(const Real x, const Real v, const Real l)
{
    return boost::math::cdf(boost::math::complement(
        boost::math::non_central_t_distribution<Real, StatsPolicy>(v, l), x));
}

float
nct_sf_float(float x, float v, float l)
{
    return nct_sf_wrap(x, v, l);
}

double
nct_sf_double(double x, double v, double l)
{
    return nct_sf_wrap(x, v, l);
}

template<typename Real>
Real
nct_isf_wrap(const Real x, const Real v, const Real l)
{
    return boost::math::quantile(boost::math::complement(
        boost::math::non_central_t_distribution<Real, StatsPolicy>(v, l), x));
}

float
nct_isf_float(float x, float v, float l)
{
    return nct_isf_wrap(x, v, l);
}

double
nct_isf_double(double x, double v, double l)
{
    return nct_isf_wrap(x, v, l);
}

float
nct_mean_float(float v, float l)
{
    RETURN_NAN(v, 1.0);
    return boost::math::mean(boost::math::non_central_t_distribution<float, StatsPolicy>(v, l));
}

double
nct_mean_double(double v, double l)
{
    RETURN_NAN(v, 1.0);
    return boost::math::mean(boost::math::non_central_t_distribution<double, StatsPolicy>(v, l));
}

float
nct_variance_float(float v, float l)
{
    RETURN_NAN(v, 2.0);
    return boost::math::variance(boost::math::non_central_t_distribution<float, StatsPolicy>(v, l));
}

double
nct_variance_double(double v, double l)
{
    RETURN_NAN(v, 2.0);
    return boost::math::variance(boost::math::non_central_t_distribution<double, StatsPolicy>(v, l));
}

float
nct_skewness_float(float v, float l)
{
    RETURN_NAN(v, 3.0);
    return boost::math::skewness(boost::math::non_central_t_distribution<float, StatsPolicy>(v, l));
}

double
nct_skewness_double(double v, double l)
{
    RETURN_NAN(v, 3.0);
    return boost::math::skewness(boost::math::non_central_t_distribution<double, StatsPolicy>(v, l));
}

float
nct_kurtosis_excess_float(float v, float l)
{
    RETURN_NAN(v, 4.0);
    return boost::math::kurtosis_excess(boost::math::non_central_t_distribution<float, StatsPolicy>(v, l));
}

double
nct_kurtosis_excess_double(double v, double l)
{
    RETURN_NAN(v, 4.0);
    return boost::math::kurtosis_excess(boost::math::non_central_t_distribution<double, StatsPolicy>(v, l));
}

template<typename Real>
Real
skewnorm_cdf_wrap(const Real x, const Real l, const Real sc, const Real sh)
{
    if (std::isfinite(x)) {
        return boost::math::cdf(
            boost::math::skew_normal_distribution<Real, StatsPolicy>(l, sc, sh), x);
    }
    // -inf => 0, inf => 1
    return 1 - std::signbit(x);
}

float
skewnorm_cdf_float(float x, float l, float sc, float sh)
{
    return skewnorm_cdf_wrap(x, l, sc, sh);
}

double
skewnorm_cdf_double(double x, double l, double sc, double sh)
{
    return skewnorm_cdf_wrap(x, l, sc, sh);
}

template<typename Real>
Real
skewnorm_ppf_wrap(const Real x, const Real l, const Real sc, const Real sh)
{
    return boost::math::quantile(
        boost::math::skew_normal_distribution<Real, StatsPolicy>(l, sc, sh), x);
}

float
skewnorm_ppf_float(float x, float l, float sc, float sh)
{
    return skewnorm_ppf_wrap(x, l, sc, sh);
}

double
skewnorm_ppf_double(double x, double l, double sc, double sh)
{
    return skewnorm_ppf_wrap(x, l, sc, sh);
}

template<typename Real>
Real
skewnorm_isf_wrap(const Real x, const Real l, const Real sc, const Real sh)
{
    return boost::math::quantile(boost::math::complement(
        boost::math::skew_normal_distribution<Real, StatsPolicy>(l, sc, sh), x));
}

float
skewnorm_isf_float(float x, float l, float sc, float sh)
{
    return skewnorm_isf_wrap(x, l, sc, sh);
}

double
skewnorm_isf_double(double x, double l, double sc, double sh)
{
    return skewnorm_isf_wrap(x, l, sc, sh);
}

template<typename Real>
Real
binom_pmf_wrap(const Real x, const Real n, const Real p)
{
    if (std::isfinite(x)) {
        return boost::math::pdf(
            boost::math::binomial_distribution<Real, StatsPolicy>(n, p), x);
    }
    return NAN; // inf or -inf returns NAN
}

float
binom_pmf_float(float x, float n, float p)
{
    return binom_pmf_wrap(x, n, p);
}

double
binom_pmf_double(double x, double n, double p)
{
    return binom_pmf_wrap(x, n, p);
}

template<typename Real>
Real
binom_cdf_wrap(const Real x, const Real n, const Real p)
{
    if (std::isfinite(x)) {
        return boost::math::cdf(
            boost::math::binomial_distribution<Real, StatsPolicy>(n, p), x);
    }
    // -inf => 0, inf => 1
    return 1 - std::signbit(x);
}

float
binom_cdf_float(float x, float n, float p)
{
    return binom_cdf_wrap(x, n, p);
}

double
binom_cdf_double(double x, double n, double p)
{
    return binom_cdf_wrap(x, n, p);
}

template<typename Real>
Real
binom_ppf_wrap(const Real x, const Real n, const Real p)
{
    return boost::math::quantile(
        boost::math::binomial_distribution<Real, StatsPolicy>(n, p), x);
}

float
binom_ppf_float(float x, float n, float p)
{
    return binom_ppf_wrap(x, n, p);
}

double
binom_ppf_double(double x, double n, double p)
{
    return binom_ppf_wrap(x, n, p);
}

template<typename Real>
Real
binom_sf_wrap(const Real x, const Real n, const Real p)
{
    return boost::math::cdf(boost::math::complement(
        boost::math::binomial_distribution<Real, StatsPolicy>(n, p), x));
}

float
binom_sf_float(float x, float n, float p)
{
    return binom_sf_wrap(x, n, p);
}

double
binom_sf_double(double x, double n, double p)
{
    return binom_sf_wrap(x, n, p);
}

template<typename Real>
Real
binom_isf_wrap(const Real x, const Real n, const Real p)
{
    return boost::math::quantile(boost::math::complement(
        boost::math::binomial_distribution<Real, StatsPolicy>(n, p), x));
}

float
binom_isf_float(float x, float n, float p)
{
    return binom_isf_wrap(x, n, p);
}

double
binom_isf_double(double x, double n, double p)
{
    return binom_isf_wrap(x, n, p);
}

template<typename Real>
Real
nbinom_pmf_wrap(const Real x, const Real r, const Real p)
{
    if (std::isfinite(x)) {
        return boost::math::pdf(
            boost::math::negative_binomial_distribution<Real, StatsPolicy>(r, p), x);
    }
    return NAN; // inf or -inf returns NAN
}

float
nbinom_pmf_float(float x, float r, float p)
{
    return nbinom_pmf_wrap(x, r, p);
}

double
nbinom_pmf_double(double x, double r, double p)
{
    return nbinom_pmf_wrap(x, r, p);
}

template<typename Real>
Real
nbinom_cdf_wrap(const Real x, const Real r, const Real p)
{
    if (std::isfinite(x)) {
        return boost::math::cdf(
            boost::math::negative_binomial_distribution<Real, StatsPolicy>(r, p), x);
    }
    // -inf => 0, inf => 1
    return 1 - std::signbit(x);
}

float
nbinom_cdf_float(float x, float r, float p)
{
    return nbinom_cdf_wrap(x, r, p);
}

double
nbinom_cdf_double(double x, double r, double p)
{
    return nbinom_cdf_wrap(x, r, p);
}

template<typename Real>
Real
nbinom_ppf_wrap(const Real x, const Real r, const Real p)
{
    return boost::math::quantile(
        boost::math::negative_binomial_distribution<Real, StatsPolicy>(r, p), x);
}

float
nbinom_ppf_float(float x, float r, float p)
{
    return nbinom_ppf_wrap(x, r, p);
}

double
nbinom_ppf_double(double x, double r, double p)
{
    return nbinom_ppf_wrap(x, r, p);
}

template<typename Real>
Real
nbinom_sf_wrap(const Real x, const Real r, const Real p)
{
    return boost::math::cdf(boost::math::complement(
        boost::math::negative_binomial_distribution<Real, StatsPolicy>(r, p), x));
}

float
nbinom_sf_float(float x, float r, float p)
{
    return nbinom_sf_wrap(x, r, p);
}

double
nbinom_sf_double(double x, double r, double p)
{
    return nbinom_sf_wrap(x, r, p);
}

template<typename Real>
Real
nbinom_isf_wrap(const Real x, const Real r, const Real p)
{
    return boost::math::quantile(boost::math::complement(
        boost::math::negative_binomial_distribution<Real, StatsPolicy>(r, p), x));
}

float
nbinom_isf_float(float x, float r, float p)
{
    return nbinom_isf_wrap(x, r, p);
}

double
nbinom_isf_double(double x, double r, double p)
{
    return nbinom_isf_wrap(x, r, p);
}

float
nbinom_mean_float(float r, float p)
{
    return boost::math::mean(boost::math::negative_binomial_distribution<float, StatsPolicy>(r, p));
}

double
nbinom_mean_double(double r, double p)
{
    return boost::math::mean(boost::math::negative_binomial_distribution<double, StatsPolicy>(r, p));
}

float
nbinom_variance_float(float r, float p)
{
    return boost::math::variance(boost::math::negative_binomial_distribution<float, StatsPolicy>(r, p));
}

double
nbinom_variance_double(double r, double p)
{
    return boost::math::variance(boost::math::negative_binomial_distribution<double, StatsPolicy>(r, p));
}

float
nbinom_skewness_float(float r, float p)
{
    return boost::math::skewness(boost::math::negative_binomial_distribution<float, StatsPolicy>(r, p));
}

double
nbinom_skewness_double(double r, double p)
{
    return boost::math::skewness(boost::math::negative_binomial_distribution<double, StatsPolicy>(r, p));
}

float
nbinom_kurtosis_excess_float(float r, float p)
{
    return boost::math::kurtosis_excess(boost::math::negative_binomial_distribution<float, StatsPolicy>(r, p));
}

double
nbinom_kurtosis_excess_double(double r, double p)
{
    return boost::math::kurtosis_excess(boost::math::negative_binomial_distribution<double, StatsPolicy>(r, p));
}

template<typename Real>
Real
hypergeom_pmf_wrap(const Real k, const Real n, const Real N, const Real M)
{
    if (std::isfinite(k)) {
        return boost::math::pdf(
            boost::math::hypergeometric_distribution<Real, StatsPolicy>(n, N, M), k);
    }
    return NAN; // inf or -inf returns NAN
}

float
hypergeom_pmf_float(float k, float n, float N, float M)
{
    return hypergeom_pmf_wrap(k, n, N, M);
}

double
hypergeom_pmf_double(double k, double n, double N, double M)
{
    return hypergeom_pmf_wrap(k, n, N, M);
}

template<typename Real>
Real
hypergeom_cdf_wrap(const Real k, const Real n, const Real N, const Real M)
{
    if (std::isfinite(k)) {
        return boost::math::cdf(
            boost::math::hypergeometric_distribution<Real, StatsPolicy>(n, N, M), k);
    }
    // -inf => 0, inf => 1
    return 1 - std::signbit(k);
}

float
hypergeom_cdf_float(float k, float n, float N, float M)
{
    return hypergeom_cdf_wrap(k, n, N, M);
}

double
hypergeom_cdf_double(double k, double n, double N, double M)
{
    return hypergeom_cdf_wrap(k, n, N, M);
}

template<typename Real>
Real
hypergeom_sf_wrap(const Real k, const Real n, const Real N, const Real M)
{
    return boost::math::cdf(boost::math::complement(
        boost::math::hypergeometric_distribution<Real, StatsPolicy>(n, N, M), k));
}

float
hypergeom_sf_float(float k, float n, float N, float M)
{
    return hypergeom_sf_wrap(k, n, N, M);
}

double
hypergeom_sf_double(double k, double n, double N, double M)
{
    return hypergeom_sf_wrap(k, n, N, M);
}

float
hypergeom_mean_float(float n, float N, float M)
{
    return boost::math::mean(boost::math::hypergeometric_distribution<float, StatsPolicy>(n, N, M));
}

double
hypergeom_mean_double(double n, double N, double M)
{
    return boost::math::mean(boost::math::hypergeometric_distribution<double, StatsPolicy>(n, N, M));
}

float
hypergeom_variance_float(float n, float N, float M)
{
    return boost::math::variance(boost::math::hypergeometric_distribution<float, StatsPolicy>(n, N, M));
}

double
hypergeom_variance_double(double n, double N, double M)
{
    return boost::math::variance(boost::math::hypergeometric_distribution<double, StatsPolicy>(n, N, M));
}

float
hypergeom_skewness_float(float n, float N, float M)
{
    return boost::math::skewness(boost::math::hypergeometric_distribution<float, StatsPolicy>(n, N, M));
}

double
hypergeom_skewness_double(double n, double N, double M)
{
    return boost::math::skewness(boost::math::hypergeometric_distribution<double, StatsPolicy>(n, N, M));
}

#endif
