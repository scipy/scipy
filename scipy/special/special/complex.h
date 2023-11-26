#pragma once

#include "config.h"

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
