#ifndef PYTHONIC_INCLUDE_NUMPY_ARRAYEQUAL_HPP
#define PYTHONIC_INCLUDE_NUMPY_ARRAYEQUAL_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class U, class V>
  bool array_equal(U const &u, V const &v);

  DEFINE_FUNCTOR(pythonic::numpy, array_equal);
}
PYTHONIC_NS_END

#endif
