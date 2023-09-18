#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_STANDARD_GAMMA_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_STANDARD_GAMMA_HPP

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
    types::ndarray<double, pS> standard_gamma(double s, pS const &shape);

    auto standard_gamma(double s, long size)
        -> decltype(standard_gamma(s, types::array<long, 1>{{size}}));

    double standard_gamma(double s, types::none_type d = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, standard_gamma);
  }
}
PYTHONIC_NS_END

#endif
