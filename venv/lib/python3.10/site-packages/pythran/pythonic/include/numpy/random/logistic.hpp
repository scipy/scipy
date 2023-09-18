#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_LOGISTIC_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_LOGISTIC_HPP

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
    types::ndarray<double, pS> logistic(double loc, double scale,
                                        pS const &shape);

    auto logistic(double loc, double scale, long size)
        -> decltype(logistic(loc, scale, types::array<long, 1>{{size}}));

    double logistic(double loc = 0.0, double scale = 1.0,
                    types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, logistic);
  }
}
PYTHONIC_NS_END

#endif
