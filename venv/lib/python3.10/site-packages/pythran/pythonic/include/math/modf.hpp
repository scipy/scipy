#ifndef PYTHONIC_INCLUDE_MATH_MODF_HPP
#define PYTHONIC_INCLUDE_MATH_MODF_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/tuple.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace math
{
  std::tuple<double, double> modf(double x);
  DEFINE_FUNCTOR(pythonic::math, modf);
}
PYTHONIC_NS_END

#endif
