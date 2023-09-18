#ifndef PYTHONIC_INCLUDE_NUMPY_ALEN_HPP
#define PYTHONIC_INCLUDE_NUMPY_ALEN_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>
  long alen(T &&expr);

  DEFINE_FUNCTOR(pythonic::numpy, alen);
}
PYTHONIC_NS_END

#endif
