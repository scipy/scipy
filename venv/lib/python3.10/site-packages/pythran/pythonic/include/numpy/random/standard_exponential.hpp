#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_STANDARD_EXPONENTIAL_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_STANDARD_EXPONENTIAL_HPP

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
    types::ndarray<double, pS> standard_exponential(pS const &shape);

    auto standard_exponential(long size)
        -> decltype(standard_exponential(types::array<long, 1>{{size}}));

    double standard_exponential(types::none_type d = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, standard_exponential);
  }
}
PYTHONIC_NS_END

#endif
