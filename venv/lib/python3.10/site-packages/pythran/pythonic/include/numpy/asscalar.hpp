#ifndef PYTHONIC_INCLUDE_NUMPY_ASSCALAR_HPP
#define PYTHONIC_INCLUDE_NUMPY_ASSCALAR_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>
  using asscalar_result_type = typename std::conditional<
      std::is_integral<T>::value, long,
      typename std::conditional<std::is_floating_point<T>::value, double,
                                std::complex<double>>::type>::type;

  template <class E>
  asscalar_result_type<typename E::dtype> asscalar(E const &expr);

  DEFINE_FUNCTOR(pythonic::numpy, asscalar);
}
PYTHONIC_NS_END

#endif
