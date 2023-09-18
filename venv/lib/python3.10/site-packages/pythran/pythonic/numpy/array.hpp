#ifndef PYTHONIC_NUMPY_ARRAY_HPP
#define PYTHONIC_NUMPY_ARRAY_HPP

#include "pythonic/include/numpy/array.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/utils/nested_container.hpp"
#include "pythonic/types/ndarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class T, class dtype>
  typename std::enable_if<
      types::has_size<typename std::decay<T>::type>::value,
      types::ndarray<typename dtype::type,
                     types::array<long, std::decay<T>::type::value>>>::type
  array(T &&iterable, dtype d)
  {
    return {std::forward<T>(iterable)};
  }
  template <class T, class dtype>
  typename std::enable_if<
      !types::has_size<typename std::decay<T>::type>::value &&
          !types::is_dtype<typename std::decay<T>::type>::value,
      types::ndarray<typename dtype::type,
                     types::array<long, std::decay<T>::type::value>>>::type
  array(T &&iterable, dtype d)
  {
    types::list<typename std::decay<T>::type::value_type> tmp{iterable.begin(),
                                                              iterable.end()};
    return {tmp};
  }

  template <class T, class dtype>
  typename std::enable_if<
      !types::has_size<typename std::decay<T>::type>::value &&
          types::is_dtype<typename std::decay<T>::type>::value,
      typename dtype::type>::type
  array(T &&non_iterable, dtype d)
  {
    return non_iterable;
  }

  template <class dtype>
  types::ndarray<typename dtype::type,
                 types::pshape<std::integral_constant<long, 0>>>
      array(std::tuple<>, dtype)
  {
    return {types::pshape<std::integral_constant<long, 0>>{},
            types::none_type{}};
  }

  template <class T, class pS>
  types::ndarray<T, pS> array(types::ndarray<T, pS> const &arr)
  {
    return arr.copy();
  }

  template <class T, size_t N, class V, class dtype>
  types::ndarray<typename dtype::type,
                 typename types::array_base<T, N, V>::shape_t>
  array(types::array_base<T, N, V> const &a, dtype)
  {
    return {a};
  }

  template <class T, size_t N, class V, class dtype>
  types::ndarray<typename dtype::type,
                 typename types::array_base<T, N, V>::shape_t>
  array(types::array_base<T, N, V> &&a, dtype)
  {
    return {std::move(a)};
  }
}
PYTHONIC_NS_END

#endif
