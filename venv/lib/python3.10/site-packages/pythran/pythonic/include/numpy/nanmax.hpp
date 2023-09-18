#ifndef PYTHONIC_INCLUDE_NUMPY_NANMAX_HPP
#define PYTHONIC_INCLUDE_NUMPY_NANMAX_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/numpy/isnan.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  typename E::dtype nanmax(E const &expr);

  DEFINE_FUNCTOR(pythonic::numpy, nanmax);
}
PYTHONIC_NS_END

#endif
