#ifndef PYTHONIC_INCLUDE_NUMPY_ARGWHERE_HPP
#define PYTHONIC_INCLUDE_NUMPY_ARGWHERE_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  typename types::ndarray<long, types::array<long, 2>> argwhere(E const &expr);

  DEFINE_FUNCTOR(pythonic::numpy, argwhere);
}
PYTHONIC_NS_END

#endif
