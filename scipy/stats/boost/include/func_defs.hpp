#ifndef CLASS_DEF_HPP
#define CLASS_DEF_HPP

#include <utility>
#include <cmath>
#include "boost/math/distributions.hpp"

typedef boost::math::policies::policy<> Policy;

template<template <typename, typename> class Dist, class RealType, class ... Args>
RealType boost_pdf(const RealType x, const Args ... args) {
  if (std::isfinite(x)) {
    return boost::math::pdf(Dist<RealType, Policy>(args...), x);
  }
  return NAN; // inf or -inf returns NAN
}

template<template <typename, typename> class Dist, class RealType, class ... Args>
RealType boost_cdf(const RealType x, const Args ... args) {
  if (std::isfinite(x)) {
    return boost::math::cdf(Dist<RealType, Policy>(args...), x);
  }
  if (std::signbit(x)) {
    // -inf
    return 0;
  }
  // inf
  return 1;
}

template<template <typename, typename> class Dist, class RealType, class ... Args>
RealType boost_icdf(const RealType x, const Args ... args) {
  return boost::math::cdf(boost::math::complement(Dist<RealType, Policy>(args...), x));
}

template<template <typename, typename> class Dist, class RealType, class ... Args>
RealType boost_quantile(const RealType q, const Args ... args) {
  return boost::math::quantile(Dist<RealType, Policy>(args...), q);
}

template<template <typename, typename> class Dist, class RealType, class ... Args >
RealType boost_iquantile(const RealType q, const Args ... args) {
  return boost::math::quantile(boost::math::complement(Dist<RealType, Policy>(args...), q));
}

template<template <typename, typename> class Dist, class RealType, class ... Args>
RealType boost_mean(const Args ... args) {
  return boost::math::mean(Dist<RealType, Policy>(args...));
}

template<template <typename, typename> class Dist, class RealType, class ... Args>
RealType boost_variance(const Args ... args) {
  return boost::math::variance(Dist<RealType, Policy>(args...));
}

template<template <typename, typename> class Dist, class RealType, class ... Args>
RealType boost_skewness(const Args ... args) {
  return boost::math::skewness(Dist<RealType, Policy>(args...));
}

template<template <typename, typename> class Dist, class RealType, class ... Args>
RealType boost_kurtosis_excess(const Args ... args) {
  return boost::math::kurtosis_excess(Dist<RealType, Policy>(args...));
}

#endif
