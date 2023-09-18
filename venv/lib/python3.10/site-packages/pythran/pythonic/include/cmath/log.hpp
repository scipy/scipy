#ifndef PYTHONIC_INCLUDE_CMATH_LOG_HPP
#define PYTHONIC_INCLUDE_CMATH_LOG_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/complex.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  using std::log;
  double log(double x, double base);
  DEFINE_FUNCTOR(pythonic::cmath, log);
}
PYTHONIC_NS_END

#endif
