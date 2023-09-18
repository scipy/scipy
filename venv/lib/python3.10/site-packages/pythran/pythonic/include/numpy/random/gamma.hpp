#ifndef PYTHONIC_INCLUDE_NUMPY_RANDOM_GAMMA_HPP
#define PYTHONIC_INCLUDE_NUMPY_RANDOM_GAMMA_HPP

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
    types::ndarray<double, pS> gamma(double shape, double scale,
                                     pS const &array_shape);

    auto gamma(double shape, double scale, long size)
        -> decltype(gamma(shape, scale, types::array<long, 1>{{size}}));

    double gamma(double shape = 0.0, double scale = 1.0,
                 types::none_type size = {});

    DEFINE_FUNCTOR(pythonic::numpy::random, gamma);
  }
}
PYTHONIC_NS_END

#endif
