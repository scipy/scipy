#ifndef PYTHONIC_INCLUDE_NUMPY_NDIM_HPP
#define PYTHONIC_INCLUDE_NUMPY_NDIM_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/numpy/shape.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class E>
  long ndim(E const &e);

  DEFINE_FUNCTOR(pythonic::numpy, ndim)
}
PYTHONIC_NS_END

#endif
