#ifndef PYTHONIC_INCLUDE_CMATH_ASIN_HPP
#define PYTHONIC_INCLUDE_CMATH_ASIN_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/complex.hpp"

#include <cmath>

PYTHONIC_NS_BEGIN

namespace cmath
{
  DEFINE_FUNCTOR_2(asin, std::asin);
}
PYTHONIC_NS_END

#endif
