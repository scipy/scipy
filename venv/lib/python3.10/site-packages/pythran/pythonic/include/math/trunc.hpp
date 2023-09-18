#ifndef PYTHONIC_INCLUDE_MATH_TRUNC_HPP
#define PYTHONIC_INCLUDE_MATH_TRUNC_HPP

#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace math
{
  template <class T>
  long trunc(T x);

  DEFINE_FUNCTOR(pythonic::math, trunc);
}
PYTHONIC_NS_END

#endif
