#ifndef PYTHONIC_INCLUDE_NUMPY_ARRAY_HPP
#define PYTHONIC_INCLUDE_NUMPY_ARRAY_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/nested_container.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T,
            class dtype = types::dtype_t<typename std::decay<T>::type::dtype>>
  typename std::enable_if<
      types::has_size<typename std::decay<T>::type>::value,
      types::ndarray<typename dtype::type,
                     types::array<long, std::decay<T>::type::value>>>::type
  array(T &&iterable, dtype d = dtype());
  template <class T,
            class dtype = types::dtype_t<typename std::decay<T>::type::dtype>>
  typename std::enable_if<
      !types::has_size<typename std::decay<T>::type>::value &&
          !types::is_dtype<typename std::decay<T>::type>::value,
      types::ndarray<typename dtype::type,
                     types::array<long, std::decay<T>::type::value>>>::type
  array(T &&iterable, dtype d = dtype());

  template <class T, class dtype = types::dtype_t<typename types::dtype_of<
                         typename std::decay<T>::type>::type>>
  typename std::enable_if<
      !types::has_size<typename std::decay<T>::type>::value &&
          types::is_dtype<typename std::decay<T>::type>::value,
      typename dtype::type>::type
  array(T &&non_iterable, dtype d = dtype());

  template <class dtype>
  types::ndarray<typename dtype::type,
                 types::pshape<std::integral_constant<long, 0>>>
      array(std::tuple<>, dtype);

  template <class T, class pS>
  types::ndarray<T, pS> array(types::ndarray<T, pS> const &arr);

  template <class T, size_t N, class V, class dtype = types::dtype_of<T>>
  types::ndarray<typename dtype::type,
                 typename types::array_base<T, N, V>::shape_t>
  array(types::array_base<T, N, V> const &, dtype d = dtype());

  template <class T, size_t N, class V, class dtype = types::dtype_of<T>>
  types::ndarray<typename dtype::type,
                 typename types::array_base<T, N, V>::shape_t>
  array(types::array_base<T, N, V> &&, dtype d = dtype());

  DEFINE_FUNCTOR(pythonic::numpy, array);
}
PYTHONIC_NS_END

#endif
