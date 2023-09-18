#ifndef PYTHONIC_INCLUDE_NUMPY_NANARGMIN_HPP
#define PYTHONIC_INCLUDE_NUMPY_NANARGMIN_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  long nanargmin(E const &expr);

  DEFINE_FUNCTOR(pythonic::numpy, nanargmin);
}
PYTHONIC_NS_END

#endif
