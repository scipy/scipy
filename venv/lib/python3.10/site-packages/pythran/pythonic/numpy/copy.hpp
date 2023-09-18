#ifndef PYTHONIC_NUMPY_COPY_HPP
#define PYTHONIC_NUMPY_COPY_HPP

#include "pythonic/include/numpy/copy.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/numpy_conversion.hpp"
#include "pythonic/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  // list case
  template <class E>
  typename std::enable_if<
      !types::is_array<E>::value && !types::is_dtype<E>::value,
      types::ndarray<typename E::dtype, types::array<long, E::value>>>::type
  copy(E const &v)
  {
    return {v};
  }

  // scalar / complex case
  template <class E>
  auto copy(E const &v) ->
      typename std::enable_if<types::is_dtype<E>::value, E>::type
  {
    return v;
  }

  // No copy is required for numpy_expr
  template <class E>
  auto copy(E &&v) ->
      typename std::enable_if<types::is_array<E>::value,
                              decltype(std::forward<E>(v))>::type
  {
    return std::forward<E>(v);
  }

  // ndarray case
  template <class T, class pS>
  types::ndarray<T, pS> copy(types::ndarray<T, pS> const &a)
  {
    return a.copy();
  }

  // transposed ndarray case
  template <class T, class pS>
  types::numpy_texpr<types::ndarray<T, pS>>
  copy(types::numpy_texpr<types::ndarray<T, pS>> const &a)
  {
    return a.arg.copy();
  }
}
PYTHONIC_NS_END

#endif
