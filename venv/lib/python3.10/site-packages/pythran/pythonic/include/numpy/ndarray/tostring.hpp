#ifndef PYTHONIC_INCLUDE_NUMPY_NDARRAY_TOSTRING_HPP
#define PYTHONIC_INCLUDE_NUMPY_NDARRAY_TOSTRING_HPP

#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/types/str.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace ndarray
  {
    template <class T, class pS>
    types::str tostring(types::ndarray<T, pS> const &expr);

    NUMPY_EXPR_TO_NDARRAY0_DECL(tostring);
    DEFINE_FUNCTOR(pythonic::numpy::ndarray, tostring);
  }
}
PYTHONIC_NS_END
#endif
