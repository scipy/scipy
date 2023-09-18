#ifndef PYTHONIC_INCLUDE_NUMPY_ARGMIN_HPP
#define PYTHONIC_INCLUDE_NUMPY_ARGMIN_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  long argmin(E const &expr);

  template <class E>
  types::ndarray<long, types::array<long, E::value - 1>> argmin(E const &expr,
                                                                long axis);

  DEFINE_FUNCTOR(pythonic::numpy, argmin);
}
PYTHONIC_NS_END

#endif
