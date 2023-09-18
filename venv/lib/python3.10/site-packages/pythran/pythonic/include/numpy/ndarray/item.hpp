#ifndef PYTHONIC_INCLUDE_NUMPY_NDARRAY_ITEM_HPP
#define PYTHONIC_INCLUDE_NUMPY_NDARRAY_ITEM_HPP

#include "pythonic/utils/functor.hpp"
#include "pythonic/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace ndarray
  {

    template <class T, class pS>
    T item(types::ndarray<T, pS> const &expr, long i);

    template <class E, size_t N>
    auto item(E &&expr, types::array<long, N> const &i) -> decltype(expr[i]);

    // only for compatibility purpose, very bad impl
    template <class E>
    typename std::decay<E>::type::dtype item(E &&expr, long i);

    DEFINE_FUNCTOR(pythonic::numpy::ndarray, item);
  }
}
PYTHONIC_NS_END

#endif
