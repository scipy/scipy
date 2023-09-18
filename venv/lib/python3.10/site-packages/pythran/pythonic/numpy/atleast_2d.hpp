#ifndef PYTHONIC_NUMPY_ATLEAST2D_HPP
#define PYTHONIC_NUMPY_ATLEAST2D_HPP

#include "pythonic/include/numpy/atleast_2d.hpp"

#include "pythonic/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>
  typename std::enable_if<
      types::is_dtype<T>::value,
      types::ndarray<T, types::pshape<std::integral_constant<long, 1>,
                                      std::integral_constant<long, 1>>>>::type
  atleast_2d(T t)
  {
    return {types::pshape<std::integral_constant<long, 1>,
                          std::integral_constant<long, 1>>(),
            t};
  }

  template <class T>
          auto atleast_2d(T const &t) ->
          typename std::enable_if < (!types::is_dtype<T>::value) &&
      T::value<2, types::ndarray<
                      typename T::dtype,
                      types::pshape<std::integral_constant<long, 1>,
                                    typename std::tuple_element<
                                        0, typename T::shape_t>::type>>>::type
  {
    return t.reshape(types::pshape<
        std::integral_constant<long, 1>,
        typename std::tuple_element<0, typename T::shape_t>::type>(
        std::integral_constant<long, 1>(), t.template shape<0>()));
  }

  template <class T>
  auto atleast_2d(T &&t) -> typename std::enable_if<
      (!types::is_dtype<typename std::remove_cv<
          typename std::remove_reference<T>::type>::type>::value) &&
          std::decay<T>::type::value >= 2,
      decltype(std::forward<T>(t))>::type
  {
    return std::forward<T>(t);
  }
}
PYTHONIC_NS_END

#endif
