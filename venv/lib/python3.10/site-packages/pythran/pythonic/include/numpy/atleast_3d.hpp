#ifndef PYTHONIC_INCLUDE_NUMPY_ATLEAST3D_HPP
#define PYTHONIC_INCLUDE_NUMPY_ATLEAST3D_HPP

#include "pythonic/include/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T>
  typename std::enable_if<
      types::is_dtype<T>::value,
      types::ndarray<T, types::pshape<std::integral_constant<long, 1>,
                                      std::integral_constant<long, 1>,
                                      std::integral_constant<long, 1>>>>::type
  atleast_3d(T t);
  template <class T>
  auto atleast_3d(T const &t) -> typename std::enable_if<
      (!types::is_dtype<T>::value) && (T::value == 1),
      types::ndarray<typename T::dtype,
                     types::pshape<std::integral_constant<long, 1>,
                                   typename std::tuple_element<
                                       0, typename T::shape_t>::type,
                                   std::integral_constant<long, 1>>>>::type;

  template <class T>
  auto atleast_3d(T const &t) -> typename std::enable_if<
      (!types::is_dtype<T>::value) && (T::value == 2),
      types::ndarray<
          typename T::dtype,
          types::pshape<
              typename std::tuple_element<0, typename T::shape_t>::type,
              typename std::tuple_element<1, typename T::shape_t>::type,
              std::integral_constant<long, 1>>>>::type;

  template <class T>
  auto atleast_3d(T const &t) ->
      typename std::enable_if<(!types::is_dtype<T>::value) && T::value >= 3,
                              decltype(asarray(t))>::type;

  DEFINE_FUNCTOR(pythonic::numpy, atleast_3d);
}
PYTHONIC_NS_END

#endif
