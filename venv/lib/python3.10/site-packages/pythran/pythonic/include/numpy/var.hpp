#ifndef PYTHONIC_INCLUDE_NUMPY_VAR_HPP
#define PYTHONIC_INCLUDE_NUMPY_VAR_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/builtins/None.hpp"
#include "pythonic/include/numpy/add.hpp"
#include "pythonic/include/numpy/mean.hpp"
#include "pythonic/include/numpy/sum.hpp"

#include <algorithm>

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  using var_type = typename std::conditional<
      std::is_integral<typename E::dtype>::value, double,
      decltype(std::real(std::declval<typename E::dtype>()))>::type;

  template <class E>
  auto var(E const &expr, types::none_type axis = builtins::None,
           types::none_type dtype = builtins::None,
           types::none_type out = builtins::None, long ddof = 0)
      -> decltype(var_type<E>(std::real(mean(expr))));

  template <class E>
  auto var(E const &expr, long axis, types::none_type dtype = builtins::None,
           types::none_type out = builtins::None, long ddof = 0) ->
      typename assignable<decltype(var_type<E>() * mean(expr, axis))>::type;

  DEFINE_FUNCTOR(pythonic::numpy, var);
}
PYTHONIC_NS_END

#endif
