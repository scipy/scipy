#ifndef PYTHONIC_INCLUDE_RANDOM_UNIFORM_HPP
#define PYTHONIC_INCLUDE_RANDOM_UNIFORM_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/random/random.hpp"

PYTHONIC_NS_BEGIN

namespace random
{
  double uniform(double a, double b);

  DEFINE_FUNCTOR(pythonic::random, uniform);
}
PYTHONIC_NS_END

#endif
