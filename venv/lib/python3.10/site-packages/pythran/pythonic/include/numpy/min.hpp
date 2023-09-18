#ifndef PYTHONIC_INCLUDE_NUMPY_MIN_HPP
#define PYTHONIC_INCLUDE_NUMPY_MIN_HPP

#include "pythonic/include/numpy/reduce.hpp"
#include "pythonic/include/operator_/imin.hpp"
#include "pythonic/include/utils/functor.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{

  template <class... Args>
  auto min(Args &&... args) -> decltype(
      reduce<operator_::functor::imin>(std::forward<Args>(args)...));

  DEFINE_FUNCTOR(pythonic::numpy, min);
}
PYTHONIC_NS_END

#endif
