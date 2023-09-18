#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_LOGNORMAL_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_LOGNORMAL_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/NoneType.hpp"
#include "pythonic/include/types/tuple.hpp"

PYTHONIC_NS_BEGIN
namespace numpy
{
  namespace random
  {
    template <class pS>
    types::ndarray<double, pS> lognormal(double mean, double sigma,
                                         pS const &shape);

    auto lognormal(double mean, double sigma, long size)
        -> decltype(lognormal(mean, sigma, types::array<long, 1>{{size}}));

    double lognormal(double mean = 0.0, double sigma = 1.0,
                     types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, lognormal);
  }
}
PYTHONIC_NS_END

#endif
