#ifndef PYTHONIC_NUMPY_BROADCAST_TO_HPP
#define PYTHONIC_NUMPY_BROADCAST_TO_HPP

#include "pythonic/include/numpy/broadcast_to.hpp"
#include "pythonic/utils/functor.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/numpy/empty.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E, class pS>
  auto broadcast_to(E const &expr, pS shape)
      -> decltype(numpy::functor::empty{}(
          shape, typename types::dtype_t<typename types::dtype_of<E>::type>{}))
  {
    using dtype = typename types::dtype_of<E>::type;
    using BExpr =
        typename std::conditional<std::is_scalar<E>::value,
                                  types::broadcast<E, dtype>, E const &>::type;
    auto out = numpy::functor::empty{}(shape, typename types::dtype_t<dtype>{});
    using array_type = decltype(out);
    BExpr bexpr = expr;
    utils::broadcast_copy<array_type, E, array_type::value,
                          array_type::value -
                              utils::nested_container_depth<E>::value,
                          std::remove_reference<BExpr>::type::is_vectorizable>(
        out, bexpr);
    return out;
  }
}
PYTHONIC_NS_END

#endif
