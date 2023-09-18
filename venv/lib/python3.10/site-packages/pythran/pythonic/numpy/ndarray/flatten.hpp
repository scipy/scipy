#ifndef PYTHONIC_NUMPY_NDARRAY_FLATTEN_HPP
#define PYTHONIC_NUMPY_NDARRAY_FLATTEN_HPP

#include "pythonic/include/numpy/ndarray/flatten.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace ndarray
  {
    template <class T, class pS>
    types::ndarray<T, types::pshape<long>>
    flatten(types::ndarray<T, pS> const &a)
    {
      return {a.mem, types::pshape<long>{a.flat_size()}};
    }

    NUMPY_EXPR_TO_NDARRAY0_IMPL(flatten);
  }
}
PYTHONIC_NS_END

#endif
