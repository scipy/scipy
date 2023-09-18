#ifndef PYTHONIC_INCLUDE_NUMPY_BROADCAST_TO_HPP
#define PYTHONIC_INCLUDE_NUMPY_BROADCAST_TO_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/numpy/empty.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class pS>
  auto broadcast_to(E const &expr, pS shape)
      -> decltype(numpy::functor::empty{}(
          shape, typename types::dtype_t<typename types::dtype_of<E>::type>{}));

  DEFINE_FUNCTOR(pythonic::numpy, broadcast_to);
}
PYTHONIC_NS_END

#endif
