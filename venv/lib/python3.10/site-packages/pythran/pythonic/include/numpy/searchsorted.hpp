#ifndef PYTHONIC_INCLUDE_NUMPY_SEARCHSORTED_HPP
#define PYTHONIC_INCLUDE_NUMPY_SEARCHSORTED_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/numpy_conversion.hpp"
#include "pythonic/include/utils/int_.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/str.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class T, class U>
  typename std::enable_if<!types::is_numexpr_arg<T>::value, long>::type
  searchsorted(U const &a, T const &v, types::str const &side = "left");

  template <class E, class T>
  typename std::enable_if<
      types::is_numexpr_arg<E>::value,
      types::ndarray<long, types::array<long, E::value>>>::type
  searchsorted(T const &a, E const &v, types::str const &side = "left");

  DEFINE_FUNCTOR(pythonic::numpy, searchsorted);
}
PYTHONIC_NS_END

#endif
