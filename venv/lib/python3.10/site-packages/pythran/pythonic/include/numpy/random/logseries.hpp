#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_LOGSERIES_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_LOGSERIES_HPP

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
    types::ndarray<double, pS> logseries(double loc, pS const &shape);

    auto logseries(double loc, long size)
        -> decltype(logseries(loc, types::array<long, 1>{{size}}));

    double logseries(double loc, types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, logseries);
  }
}
PYTHONIC_NS_END

#endif
