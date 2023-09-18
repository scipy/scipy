#ifndef PYTHONIC_INCLUDE_NUMPY_VDOT_HPP
#define PYTHONIC_INCLUDE_NUMPY_VDOT_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/numpy/dot.hpp"
#include "pythonic/include/numpy/asarray.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class U, class V>
  auto vdot(U const &u, V const &v)
      -> decltype(functor::dot{}(functor::asarray{}(u).flat(),
                                 functor::asarray{}(v).flat()));

  DEFINE_FUNCTOR(pythonic::numpy, vdot);
}
PYTHONIC_NS_END

#endif
