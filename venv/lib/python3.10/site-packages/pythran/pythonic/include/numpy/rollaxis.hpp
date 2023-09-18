#ifndef PYTHONIC_INCLUDE_NUMPY_ROLLAXIS_HPP
#define PYTHONIC_INCLUDE_NUMPY_ROLLAXIS_HPP

#include "pythonic/include/numpy/transpose.hpp"
#include "pythonic/include/numpy/copy.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  types::ndarray<T, types::array<long, std::tuple_size<pS>::value>>
  rollaxis(types::ndarray<T, pS> const &a, long axis, long start = 0);

  NUMPY_EXPR_TO_NDARRAY0_DECL(rollaxis);
  DEFINE_FUNCTOR(pythonic::numpy, rollaxis);
}
PYTHONIC_NS_END

#endif
