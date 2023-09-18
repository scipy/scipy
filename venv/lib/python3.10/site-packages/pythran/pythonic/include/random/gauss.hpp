#ifndef PYTHONIC_INCLUDE_RANDOM_GAUSS_HPP
#define PYTHONIC_INCLUDE_RANDOM_GAUSS_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/random/random.hpp"

PYTHONIC_NS_BEGIN

namespace random
{

  double gauss(double mu, double sigma);

  DEFINE_FUNCTOR(pythonic::random, gauss);
}
PYTHONIC_NS_END

#endif
