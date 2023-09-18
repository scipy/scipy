#ifndef PYTHONIC_INCLUDE_NUMPY_DIAG_HPP
#define PYTHONIC_INCLUDE_NUMPY_DIAG_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"
#include "pythonic/include/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS>
  typename std::enable_if<std::tuple_size<pS>::value == 2,
                          types::ndarray<T, types::pshape<long>>>::type
  diag(types::ndarray<T, pS> const &a, long k = 0);

  template <class T, class pS>
  typename std::enable_if<std::tuple_size<pS>::value == 1,
                          types::ndarray<T, types::array<long, 2>>>::type
  diag(types::ndarray<T, pS> const &a, long k = 0);

  template <class T>
  auto diag(types::list<T> const &a, long k = 0)
      -> decltype(diag(asarray(a), k));

  NUMPY_EXPR_TO_NDARRAY0_DECL(diag);
  DEFINE_FUNCTOR(pythonic::numpy, diag);
}
PYTHONIC_NS_END

#endif
