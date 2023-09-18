#ifndef PYTHONIC_NUMPY_NDARRAY_ITEM_HPP
#define PYTHONIC_NUMPY_NDARRAY_ITEM_HPP

#include "pythonic/include/numpy/ndarray/item.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/types/ndarray.hpp"
#include "pythonic/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  namespace ndarray
  {

    template <class T, class pS>
    T item(types::ndarray<T, pS> const &expr, long i)
    {
      if (i < 0)
        i += expr.flat_size();
      return *(expr.fbegin() + i);
    }

    template <class E, size_t N>
    auto item(E &&expr, types::array<long, N> const &i) -> decltype(expr[i])
    {
      return expr[i];
    }

    // only for compatibility purpose, very bad impl
    template <class E>
    typename std::decay<E>::type::dtype item(E &&expr, long i)
    {
      if (i < 0)
        i += expr.flat_size();
      return asarray(std::forward<E>(expr)).flat()[i];
    }
  }
}
PYTHONIC_NS_END
#endif
