#ifndef PYTHONIC_INCLUDE_RANDOM_RANDINT_HPP
#define PYTHONIC_INCLUDE_RANDOM_RANDINT_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/random/randrange.hpp"

PYTHONIC_NS_BEGIN

namespace random
{

  long randint(long a, long b);

  DEFINE_FUNCTOR(pythonic::random, randint);
}
PYTHONIC_NS_END

#endif
