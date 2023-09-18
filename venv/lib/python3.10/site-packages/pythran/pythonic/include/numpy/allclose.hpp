#ifndef PYTHONIC_INCLUDE_NUMPY_ALLCLOSE_HPP
#define PYTHONIC_INCLUDE_NUMPY_ALLCLOSE_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/numpy/abs.hpp"
#include "pythonic/include/numpy/isfinite.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class U, class V>
  bool allclose(U const &u, V const &v, double rtol = 1e-5, double atol = 1e-8);

  DEFINE_FUNCTOR(pythonic::numpy, allclose);
}
PYTHONIC_NS_END

#endif
