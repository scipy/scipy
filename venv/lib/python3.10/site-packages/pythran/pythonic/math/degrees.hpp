#ifndef PYTHONIC_MATH_DEGREES_HPP
#define PYTHONIC_MATH_DEGREES_HPP

#include "pythonic/include/math/degrees.hpp"

#include "pythonic/utils/functor.hpp"
#include "pythonic/math/pi.hpp"

PYTHONIC_NS_BEGIN

namespace math
{

  template <class T>
  double degrees(T x)
  {
    return (x * 360.) / (2. * pi);
  }
}
PYTHONIC_NS_END

#endif
