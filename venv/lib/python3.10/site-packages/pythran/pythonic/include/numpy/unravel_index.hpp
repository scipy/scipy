#ifndef PYTHONIC_INCLUDE_NUMPY_UNRAVEL_INDEX_HPP
#define PYTHONIC_INCLUDE_NUMPY_UNRAVEL_INDEX_HPP

#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/str.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class S>
  typename std::enable_if<std::is_scalar<E>::value,
                          types::array<long, std::tuple_size<S>::value>>::type
  unravel_index(E const &expr, S const &shape, types::str const &order = "C");

  DEFINE_FUNCTOR(pythonic::numpy, unravel_index);
}
PYTHONIC_NS_END

#endif
