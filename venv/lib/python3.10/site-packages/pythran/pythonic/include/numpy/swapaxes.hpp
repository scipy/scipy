#ifndef PYTHONIC_INCLUDE_NUMPY_SWAPAXES_HPP
#define PYTHONIC_INCLUDE_NUMPY_SWAPAXES_HPP

#include "pythonic/include/numpy/transpose.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>
  auto swapaxes(T &&a, int axis1, int axis2) -> decltype(functor::transpose{}(
      std::forward<T>(a),
      std::declval<types::array<long, std::decay<T>::type::value>>()));

  DEFINE_FUNCTOR(pythonic::numpy, swapaxes);
}
PYTHONIC_NS_END

#endif
