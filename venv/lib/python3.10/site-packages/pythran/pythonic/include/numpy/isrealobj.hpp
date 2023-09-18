#ifndef PYTHONIC_INCLUDE_NUMPY_ISREALOBJ_HPP
#define PYTHONIC_INCLUDE_NUMPY_ISREALOBJ_HPP

#include "pythonic/include/utils/functor.hpp"
#include "pythonic/include/types/ndarray.hpp"
#include "pythonic/include/types/traits.hpp"

PYTHONIC_NS_BEGIN

namespace numpy
{
  template <class E>
  constexpr bool isrealobj(E const &expr);

  DEFINE_FUNCTOR(pythonic::numpy, isrealobj);
}
PYTHONIC_NS_END

#endif
