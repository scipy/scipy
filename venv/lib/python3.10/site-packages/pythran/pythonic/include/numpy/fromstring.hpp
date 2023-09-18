#ifndef PYTHONIC_INCLUDE_NUMPY_FROMSTRING_HPP
#define PYTHONIC_INCLUDE_NUMPY_FROMSTRING_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/numpy/float64.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/list.hpp"
#include "pythonic/include/types/str.hpp"

#include <limits>
#include <sstream>

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class dtype = functor::float64>
  types::ndarray<typename dtype::type, types::pshape<long>>
  fromstring(types::str const &string, dtype d = dtype(), long count = -1,
             types::str const &sep = {});

  DEFINE_FUNCTOR(pythonic::numpy, fromstring);
}
PYTHONIC_NS_END

#endif
