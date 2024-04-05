#ifndef BOOST_STATS_FUNCTIONS_H
#define BOOST_STATS_FUNCTIONS_H


// Override some default BOOST policies.
// These are required to ensure that the Boost function ibeta_inv
// handles extremely small p values with precision comparable to the
// Cephes incbi function.
#define BOOST_MATH_DOMAIN_ERROR_POLICY        ignore_error
#define BOOST_MATH_OVERFLOW_ERROR_POLICY      user_error
#define BOOST_MATH_EVALUATION_ERROR_POLICY    user_error
#define BOOST_MATH_PROMOTE_DOUBLE_POLICY      false

#include "boost/math/distributions.hpp"
#include <boost/math/distributions/inverse_gaussian.hpp>

// Round up to achieve correct ppf(cdf) round-trips for discrete distributions
typedef boost::math::policies::policy<
    boost::math::policies::discrete_quantile<
        boost::math::policies::integer_round_up > > Policy;

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
            boost::math::beta_distribution<Real>(a, b), x);
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
Real
invgauss_ppf_wrap(const Real x, const Real mu, const Real s)
{
    return boost::math::quantile(
        boost::math::inverse_gaussian_distribution<Real>(mu, s), x);
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
        boost::math::inverse_gaussian_distribution<Real>(mu, s), x));
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
ncx2_pdf_wrap(const Real x, const Real k, const Real l)
{
    if (std::isfinite(x)) {
        return boost::math::pdf(
            boost::math::non_central_chi_squared_distribution<Real>(k, l), x);
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
            boost::math::non_central_chi_squared_distribution<Real>(k, l), x);
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
        boost::math::non_central_chi_squared_distribution<Real>(k, l), x);
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
        boost::math::non_central_chi_squared_distribution<Real>(k, l), x));
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
        boost::math::non_central_chi_squared_distribution<Real>(k, l), x));
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
            boost::math::non_central_f_distribution<Real>(v1, v2, l), x);
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
            boost::math::non_central_f_distribution<Real>(v1, v2, l), x);
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
    return boost::math::quantile(
        boost::math::non_central_f_distribution<Real>(v1, v2, l), x);
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
        boost::math::non_central_f_distribution<Real>(v1, v2, l), x));
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
        boost::math::non_central_f_distribution<Real>(v1, v2, l), x));
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
    return boost::math::mean(boost::math::non_central_f_distribution<float>(v1, v2, l));
}

double
ncf_mean_double(double v1, double v2, double l)
{
    RETURN_NAN(v2, 2.0);
    return boost::math::mean(boost::math::non_central_f_distribution<double>(v1, v2, l));
}

float
ncf_variance_float(float v1, float v2, float l)
{
    RETURN_NAN(v2, 4.0);
    return boost::math::variance(boost::math::non_central_f_distribution<float>(v1, v2, l));
}

double
ncf_variance_double(double v1, double v2, double l)
{
    RETURN_NAN(v2, 4.0);
    return boost::math::variance(boost::math::non_central_f_distribution<double>(v1, v2, l));
}

float
ncf_skewness_float(float v1, float v2, float l)
{
    RETURN_NAN(v2, 6.0);
    return boost::math::skewness(boost::math::non_central_f_distribution<float>(v1, v2, l));
}

double
ncf_skewness_double(double v1, double v2, double l)
{
    RETURN_NAN(v2, 6.0);
    return boost::math::skewness(boost::math::non_central_f_distribution<double>(v1, v2, l));
}

float
ncf_kurtosis_excess_float(float v1, float v2, float l)
{
    RETURN_NAN(v2, 8.0);
    return boost::math::kurtosis_excess(boost::math::non_central_f_distribution<float>(v1, v2, l));
}

double
ncf_kurtosis_excess_double(double v1, double v2, double l)
{
    RETURN_NAN(v2, 8.0);
    return boost::math::kurtosis_excess(boost::math::non_central_f_distribution<double>(v1, v2, l));
}

template<typename Real>
Real
nct_cdf_wrap(const Real x, const Real v, const Real l)
{
    if (std::isfinite(x)) {
        return boost::math::cdf(
            boost::math::non_central_t_distribution<Real>(v, l), x);
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
nct_ppf_wrap(const Real x, const Real v, const Real l)
{
    return boost::math::quantile(
        boost::math::non_central_t_distribution<Real>(v, l), x);
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
        boost::math::non_central_t_distribution<Real>(v, l), x));
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
        boost::math::non_central_t_distribution<Real>(v, l), x));
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
    return boost::math::mean(boost::math::non_central_t_distribution<float>(v, l));
}

double
nct_mean_double(double v, double l)
{
    RETURN_NAN(v, 1.0);
    return boost::math::mean(boost::math::non_central_t_distribution<double>(v, l));
}

float
nct_variance_float(float v, float l)
{
    RETURN_NAN(v, 2.0);
    return boost::math::variance(boost::math::non_central_t_distribution<float>(v, l));
}

double
nct_variance_double(double v, double l)
{
    RETURN_NAN(v, 2.0);
    return boost::math::variance(boost::math::non_central_t_distribution<double>(v, l));
}

float
nct_skewness_float(float v, float l)
{
    RETURN_NAN(v, 3.0);
    return boost::math::skewness(boost::math::non_central_t_distribution<float>(v, l));
}

double
nct_skewness_double(double v, double l)
{
    RETURN_NAN(v, 3.0);
    return boost::math::skewness(boost::math::non_central_t_distribution<double>(v, l));
}

float
nct_kurtosis_excess_float(float v, float l)
{
    RETURN_NAN(v, 4.0);
    return boost::math::kurtosis_excess(boost::math::non_central_t_distribution<float>(v, l));
}

double
nct_kurtosis_excess_double(double v, double l)
{
    RETURN_NAN(v, 4.0);
    return boost::math::kurtosis_excess(boost::math::non_central_t_distribution<double>(v, l));
}

template<typename Real>
Real
skewnorm_cdf_wrap(const Real x, const Real l, const Real sc, const Real sh)
{
    if (std::isfinite(x)) {
        return boost::math::cdf(
            boost::math::skew_normal_distribution<Real>(l, sc, sh), x);
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
        boost::math::skew_normal_distribution<Real>(l, sc, sh), x);
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
        boost::math::skew_normal_distribution<Real>(l, sc, sh), x));
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
            boost::math::binomial_distribution<Real, Policy>(n, p), x);
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
            boost::math::binomial_distribution<Real, Policy>(n, p), x);
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
        boost::math::binomial_distribution<Real, Policy>(n, p), x);
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
        boost::math::binomial_distribution<Real, Policy>(n, p), x));
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
        boost::math::binomial_distribution<Real, Policy>(n, p), x));
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
            boost::math::negative_binomial_distribution<Real, Policy>(r, p), x);
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
            boost::math::negative_binomial_distribution<Real, Policy>(r, p), x);
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
        boost::math::negative_binomial_distribution<Real, Policy>(r, p), x);
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
        boost::math::negative_binomial_distribution<Real, Policy>(r, p), x));
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
        boost::math::negative_binomial_distribution<Real, Policy>(r, p), x));
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
    return boost::math::mean(boost::math::negative_binomial_distribution<float, Policy>(r, p));
}

double
nbinom_mean_double(double r, double p)
{
    return boost::math::mean(boost::math::negative_binomial_distribution<double, Policy>(r, p));
}

float
nbinom_variance_float(float r, float p)
{
    return boost::math::variance(boost::math::negative_binomial_distribution<float, Policy>(r, p));
}

double
nbinom_variance_double(double r, double p)
{
    return boost::math::variance(boost::math::negative_binomial_distribution<double, Policy>(r, p));
}

float
nbinom_skewness_float(float r, float p)
{
    return boost::math::skewness(boost::math::negative_binomial_distribution<float, Policy>(r, p));
}

double
nbinom_skewness_double(double r, double p)
{
    return boost::math::skewness(boost::math::negative_binomial_distribution<double, Policy>(r, p));
}

float
nbinom_kurtosis_excess_float(float r, float p)
{
    return boost::math::kurtosis_excess(boost::math::negative_binomial_distribution<float, Policy>(r, p));
}

double
nbinom_kurtosis_excess_double(double r, double p)
{
    return boost::math::kurtosis_excess(boost::math::negative_binomial_distribution<double, Policy>(r, p));
}

template<typename Real>
Real
hypergeom_pmf_wrap(const Real k, const Real n, const Real N, const Real M)
{
    if (std::isfinite(k)) {
        return boost::math::pdf(
            boost::math::hypergeometric_distribution<Real, Policy>(n, N, M), k);
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
            boost::math::hypergeometric_distribution<Real, Policy>(n, N, M), k);
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
        boost::math::hypergeometric_distribution<Real, Policy>(n, N, M), k));
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
    return boost::math::mean(boost::math::hypergeometric_distribution<float, Policy>(n, N, M));
}

double
hypergeom_mean_double(double n, double N, double M)
{
    return boost::math::mean(boost::math::hypergeometric_distribution<double, Policy>(n, N, M));
}

float
hypergeom_variance_float(float n, float N, float M)
{
    return boost::math::variance(boost::math::hypergeometric_distribution<float, Policy>(n, N, M));
}

double
hypergeom_variance_double(double n, double N, double M)
{
    return boost::math::variance(boost::math::hypergeometric_distribution<double, Policy>(n, N, M));
}

float
hypergeom_skewness_float(float n, float N, float M)
{
    return boost::math::skewness(boost::math::hypergeometric_distribution<float, Policy>(n, N, M));
}

double
hypergeom_skewness_double(double n, double N, double M)
{
    return boost::math::skewness(boost::math::hypergeometric_distribution<double, Policy>(n, N, M));
}

#endif
