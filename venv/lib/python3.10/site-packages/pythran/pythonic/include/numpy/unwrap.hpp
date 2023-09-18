#ifndef PYTHONIC_INCLUDE_NUMPY_UNWRAP_HPP
#define PYTHONIC_INCLUDE_NUMPY_UNWRAP_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/utils/int_.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/numpy/pi.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  types::ndarray<double, typename E::shape_t> unwrap(E const &expr,
                                                     double discont = pi);

  DEFINE_FUNCTOR(pythonic::numpy, unwrap)
}
PYTHONIC_NS_END

#endif
