#ifndef PYTHONIC_INCLUDE_NUMPY_ALL_HPP
#define PYTHONIC_INCLUDE_NUMPY_ALL_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/numpy/multiply.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  typename std::enable_if<types::is_numexpr_arg<E>::value, bool>::type
  all(E const &expr, types::none_type _ = types::none_type());

  template <class E>
  typename std::enable_if<
      std::is_scalar<E>::value || types::is_complex<E>::value, bool>::type
  all(E const &expr, types::none_type _ = types::none_type());

  template <class E>
  auto all(E const &array, long axis) ->
      typename std::enable_if<std::is_scalar<E>::value ||
                                  types::is_complex<E>::value,
                              decltype(all(array))>::type;

  template <class E>
  auto all(E const &array, long axis) ->
      typename std::enable_if<E::value == 1, decltype(all(array))>::type;

  template <class E>
  typename std::enable_if<
      E::value != 1,
      types::ndarray<typename E::dtype, types::array<long, E::value - 1>>>::type
  all(E const &array, long axis);

  DEFINE_FUNCTOR(pythonic::numpy, all);
}
PYTHONIC_NS_END

#endif
