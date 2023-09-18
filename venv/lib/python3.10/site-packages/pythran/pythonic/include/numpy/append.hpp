#ifndef PYTHONIC_INCLUDE_NUMPY_APPEND_HPP
#define PYTHONIC_INCLUDE_NUMPY_APPEND_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class pS, class F>
  typename std::enable_if<
      !types::is_dtype<F>::value,
      types::ndarray<
          typename __combined<T, typename types::dtype_of<F>::type>::type,
          types::pshape<long>>>::type
  append(types::ndarray<T, pS> const &nto, F const &data);

  template <class T, class pS, class F>
  typename std::enable_if<
      types::is_dtype<F>::value,
      types::ndarray<
          typename __combined<T, typename types::dtype_of<F>::type>::type,
          types::pshape<long>>>::type
  append(types::ndarray<T, pS> const &nto, F const &data);

  template <class T, class F>
  types::ndarray<typename __combined<typename types::dtype_of<T>::type,
                                     typename types::dtype_of<F>::type>::type,
                 types::pshape<long>>
  append(T const &to, F const &data);

  DEFINE_FUNCTOR(pythonic::numpy, append);
}
PYTHONIC_NS_END

#endif
