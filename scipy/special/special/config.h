#pragma once

#ifdef __CUDACC__
#define SPECFUN_HOST_DEVICE __host__ __device__
#else
#define SPECFUN_HOST_DEVICE
#endif

#if __has_include(<math.h>)
#include <math.h>
#else

#define M_E 2.71828182845904523536
#define M_LOG2E 1.44269504088896340736
#define M_LOG10E 0.434294481903251827651
#define M_LN2 0.693147180559945309417
#define M_LN10 2.30258509299404568402
#define M_PI 3.14159265358979323846
#define M_PI_2 1.57079632679489661923
#define M_PI_4 0.785398163397448309616
#define M_1_PI 0.318309886183790671538
#define M_2_PI 0.636619772367581343076
#define M_2_SQRTPI 1.12837916709551257390
#define M_SQRT2 1.41421356237309504880
#define M_SQRT1_2 0.707106781186547524401

#endif

#if __has_include(<cmath>)
#include <cmath>
#else

namespace std {

SPECFUN_HOST_DEVICE inline double abs(double num) { return thrust::abs(num); }

SPECFUN_HOST_DEVICE inline double exp(double num) { return thrust::exp(num); }

SPECFUN_HOST_DEVICE inline double fma(double x, double y, double z) { return x * y + z; }

SPECFUN_HOST_DEVICE inline double log(double num) { return thrust::log(num); }

SPECFUN_HOST_DEVICE inline double sqrt(double num) { return thrust::sqrt(num); }

SPECFUN_HOST_DEVICE inline bool isnan(double num) { return thrust::isnan(num); }

} // namespace std

#endif

#if __has_include(<limits>)

#include <limits>

#else

namespace std {

template <typename T>
class numeric_limits;

template <>
class numeric_limits<double> {
  public:
    SPECFUN_HOST_DEVICE static constexpr double infinity() { return CUDART_INF; }

    SPECFUN_HOST_DEVICE static constexpr double quiet_NaN() { return CUDART_NAN; }
};

} // namespace std

#endif

#if __has_include(<complex>)

#include <complex>

#else

namespace std {

template <typename T>
using complex = thrust::complex<T>;

template <typename T>
SPECFUN_HOST_DEVICE T abs(const complex<T> &z) {
    return thrust::abs(z);
}

template <typename T>
SPECFUN_HOST_DEVICE complex<T> exp(const complex<T> &z) {
    return thrust::exp(z);
}

template <typename T>
SPECFUN_HOST_DEVICE complex<T> log(const complex<T> &z) {
    return thrust::log(z);
}

template <typename T>
SPECFUN_HOST_DEVICE T norm(const complex<T> &z) {
    return thrust::norm(z);
}

template <typename T>
SPECFUN_HOST_DEVICE complex<T> sqrt(const complex<T> &z) {
    return thrust::sqrt(z);
}

} // namespace std

#endif
