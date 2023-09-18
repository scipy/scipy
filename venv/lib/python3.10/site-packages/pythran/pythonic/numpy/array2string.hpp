#ifndef PYTHONIC_NUMPY_ARRAY2STRING_HPP
#define PYTHONIC_NUMPY_ARRAY2STRING_HPP

#include "pythonic/include/numpy/array2string.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/str.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  types::str array2string(E &&a)
  {
    std::ostringstream oss;
    oss << std::forward<E>(a);
    return oss.str();
  }
}
PYTHONIC_NS_END

#endif
