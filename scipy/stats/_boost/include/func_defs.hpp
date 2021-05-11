#ifndef CLASS_DEF_HPP
#define CLASS_DEF_HPP

#include <utility>
#include <cmath>
#include "boost/math/distributions.hpp"

// Round up to achieve correct ppf(cdf) round-trips for discrete distributions
typedef boost::math::policies::policy<
    boost::math::policies::discrete_quantile<
        boost::math::policies::integer_round_up > > Policy;

template<template <typename, typename> class Dst, class RealType, class...Args>
RealType
boost_pdf(const RealType x, const Args ... args)
{
    if (std::isfinite(x)) {
        return boost::math::pdf(Dst<RealType, Policy>(args...), x);
    }
    return NAN; // inf or -inf returns NAN
}

// patch for boost::math::beta_distribution throwing exception for
// x = 1, beta < 1 as well as x = 0, alpha < 1
template<template <typename, typename> class Dst, class RealType, class...Args>
RealType
boost_pdf_beta(const RealType x, const RealType a, const RealType b)
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
	    boost::math::beta_distribution<RealType, Policy>(a, b), x);
    }
    return NAN;
}

template<template <typename, typename> class Dst, class RealType, class...Args>
RealType
boost_cdf(const RealType x, const Args ... args)
{
    if (std::isfinite(x)) {
        return boost::math::cdf(Dst<RealType, Policy>(args...), x);
    }
    // -inf => 0, inf => 1
    return 1 - std::signbit(x);
}

template<template <typename, typename> class Dst, class RealType, class...Args>
RealType
boost_sf(const RealType x, const Args ... args)
{
    return boost::math::cdf(
        boost::math::complement(Dst<RealType, Policy>(args...), x));
}

template<template <typename, typename> class Dst, class RealType, class...Args>
RealType
boost_ppf(const RealType q, const Args ... args)
{
    return boost::math::quantile(Dst<RealType, Policy>(args...), q);
}

template<template <typename, typename> class Dst, class RealType, class...Args>
RealType
boost_isf(const RealType q, const Args ... args)
{
    return boost::math::quantile(
        boost::math::complement(Dst<RealType, Policy>(args...), q));
}

template<template <typename, typename> class Dst, class RealType, class...Args>
RealType
boost_mean(const Args ... args)
{
    return boost::math::mean(Dst<RealType, Policy>(args...));
}

template<template <typename, typename> class Dst, class RealType, class...Args>
RealType
boost_variance(const Args ... args) {
    return boost::math::variance(Dst<RealType, Policy>(args...));
}

template<template <typename, typename> class Dst, class RealType, class...Args>
RealType
boost_skewness(const Args ... args) {
    return boost::math::skewness(Dst<RealType, Policy>(args...));
}

template<template <typename, typename> class Dst, class RealType, class...Args>
RealType
boost_kurtosis_excess(const Args ... args) {
    return boost::math::kurtosis_excess(Dst<RealType, Policy>(args...));
}

#endif // CLASS_DEF_HPP
