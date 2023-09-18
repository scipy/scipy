#ifndef PYTHONIC_INCLUDE_NUMPY_EYE_HPP
#define PYTHONIC_INCLUDE_NUMPY_EYE_HPP

#include "pythonic/include/numpy/zeros.hpp"
#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/builtins/None.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class dtype = functor::float64>
  types::ndarray<typename dtype::type, types::array<long, 2>>
  eye(long N, long M, long k = 0, dtype d = dtype());

  template <class dtype = functor::float64>
  types::ndarray<typename dtype::type, types::array<long, 2>>
  eye(long N, types::none_type M = builtins::None, long k = 0,
      dtype d = dtype());

  DEFINE_FUNCTOR(pythonic::numpy, eye);
}
PYTHONIC_NS_END

#endif
