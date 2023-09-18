#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_RAND_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_RAND_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    template <class... T>
    types::ndarray<double, types::array<long, sizeof...(T)>> rand(T... shape);

    double rand();

    DEFINE_FUNCTOR(pythonic::numpy::random, rand);
  }
}
PYTHONIC_NS_END

#endif
