#ifndef PYTHONIC_MATH_RADIANS_HPP
#define PYTHONIC_MATH_RADIANS_HPP

#include "pythonic/include/math/radians.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/math/pi.hpp"

PYTHONIC_NS_BEGIN

namespace math
{
  template <class T>
  double radians(T x)
  {
    return (x * 2. * pi) / 360.;
  }
}
PYTHONIC_NS_END

#endif
